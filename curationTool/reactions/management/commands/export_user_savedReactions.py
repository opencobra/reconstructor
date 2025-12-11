"""
Management command to export a user's saved reactions to a CSV file.

Usage:
    python manage.py export_user_savedReactions --userid 4
    python manage.py export_user_savedReactions --userid 4 --output /path/to/output.csv

The CSV contains the following columns:
    - rxn_abbr: Reaction abbreviation (short_name)
    - description: Reaction description
    - vmh_formula: Constructed VMH formula from metabolites
    - new_mets: Comma-separated list of new metabolite abbreviations (not in VMH)
    - references: Formatted references (PMID:..., DOI:..., etc.)
    - GPR: Gene-Protein-Reaction rule
    - confidence_score: Confidence score
    - comments: Reaction comments
    - ext_links: External links formatted as source:id
"""

import csv
import json
import os
from datetime import datetime

from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from reactions.models import User, Reaction, SavedMetabolite
from reactions.utils.search_vmh import search_metabolites_vmh
from reactions.utils.gen_vmh_abbrs import gen_metabolite_abbr
from reactions.utils.add_to_vmh_utils import merge_gene_infos


class Command(BaseCommand):
    help = 'Export a user\'s saved reactions to a CSV file.'

    def add_arguments(self, parser):
        parser.add_argument(
            '--userid',
            type=int,
            required=True,
            help='The ID of the user whose reactions to export.'
        )
        parser.add_argument(
            '--output',
            type=str,
            default=None,
            help='Output file path. Defaults to user_{userid}_reactions_{timestamp}.csv'
        )

    def handle(self, *args, **options):
        user_id = options['userid']
        output_path = options['output']

        # ─────────────────────────────────────────────────────────────────────
        # 1. Fetch user and their saved reactions
        # ─────────────────────────────────────────────────────────────────────
        try:
            user = User.objects.get(pk=user_id)
        except User.DoesNotExist:
            raise CommandError(f"User with ID {user_id} does not exist.")

        reactions = user.saved_reactions.all()
        reaction_count = reactions.count()

        self.stdout.write(
            self.style.SUCCESS(
                f"\nFound {reaction_count} reactions saved for user {user_id}: {user.name or 'Unnamed'}\n"
            )
        )

        if reaction_count == 0:
            self.stdout.write(self.style.WARNING("No reactions to export."))
            return

        # ─────────────────────────────────────────────────────────────────────
        # 2. Collect all non-VMH metabolites that need abbreviation generation
        # ─────────────────────────────────────────────────────────────────────
        new_mets_to_generate = self._collect_new_metabolites(reactions)

        if new_mets_to_generate:
            self.stdout.write(
                f"Found {len(new_mets_to_generate)} new mets, generating abbrs...\n"
            )
            generated_abbrs = self._generate_abbreviations(new_mets_to_generate)
        else:
            self.stdout.write("No new metabolites requiring abbreviation generation.\n")
            generated_abbrs = {}

        # ─────────────────────────────────────────────────────────────────────
        # 3. Build CSV rows for each reaction
        # ─────────────────────────────────────────────────────────────────────
        self.stdout.write("Writing CSV\n")

        csv_rows = []
        for reaction in reactions:
            row = self._build_csv_row(reaction, generated_abbrs)
            csv_rows.append(row)

        # ─────────────────────────────────────────────────────────────────────
        # 4. Write to CSV file
        # ─────────────────────────────────────────────────────────────────────
        if output_path is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = f"user_{user_id}_reactions_{timestamp}.csv"

        self._write_csv(output_path, csv_rows)

        self.stdout.write(
            self.style.SUCCESS(f"\nExported {len(csv_rows)} reactions to: {output_path}\n")
        )

    # =========================================================================
    # Helper Methods
    # =========================================================================

    def _collect_new_metabolites(self, reactions):
        """
        Collect all metabolites that are not in VMH and need abbreviation generation.
        
        Returns:
            dict: Mapping of (met_identifier, met_type, met_name) -> None (to be filled with abbr)
        """
        new_mets = {}

        for reaction in reactions:
            # Process substrates
            substrates = self._safe_json_loads(reaction.substrates, [])
            subs_types = self._safe_json_loads(reaction.substrates_types, [])
            subs_names = self._safe_json_loads(reaction.substrates_names, [])
            subs_found = self._safe_json_loads(reaction.subs_found, [])

            for i, (met, met_type, met_name) in enumerate(zip(substrates, subs_types, subs_names)):
                in_vmh = subs_found[i] if i < len(subs_found) else False
                if not in_vmh and met_type.lower() != 'vmh':
                    key = (met, met_type, met_name)
                    if key not in new_mets:
                        new_mets[key] = None

            # Process products
            products = self._safe_json_loads(reaction.products, [])
            prods_types = self._safe_json_loads(reaction.products_types, [])
            prods_names = self._safe_json_loads(reaction.products_names, [])
            prods_found = self._safe_json_loads(reaction.prod_found, [])

            for i, (met, met_type, met_name) in enumerate(zip(products, prods_types, prods_names)):
                in_vmh = prods_found[i] if i < len(prods_found) else False
                if not in_vmh and met_type.lower() != 'vmh':
                    key = (met, met_type, met_name)
                    if key not in new_mets:
                        new_mets[key] = None

        return new_mets

    def _generate_abbreviations(self, new_mets):
        """
        Generate VMH abbreviations for new metabolites using MATLAB.
        
        Args:
            new_mets: dict with keys (met, met_type, met_name)
            
        Returns:
            dict: Same keys mapped to generated abbreviations
        """
        generated = {}
        total = len(new_mets)
        
        for idx, (met, met_type, met_name) in enumerate(new_mets.keys(), 1):
            self.stdout.write(f"  [{idx}/{total}] Generating abbr for: {met_name[:50]}...")
            try:
                abbr = gen_metabolite_abbr(
                    met, 
                    met_type, 
                    met_name, 
                    search_metabolites_vmh
                )
                generated[(met, met_type, met_name)] = abbr
                self.stdout.write(self.style.SUCCESS(f" -> {abbr}"))
            except Exception as e:
                # Fallback: use a sanitized version of the name
                fallback = self._sanitize_name_as_abbr(met_name)
                generated[(met, met_type, met_name)] = fallback
                self.stdout.write(self.style.WARNING(f" -> {fallback} (fallback, error: {e})"))
        
        return generated

    def _build_csv_row(self, reaction, generated_abbrs):
        """
        Build a single CSV row dict for a reaction.
        """
        # Get metabolite abbreviations and identify new mets
        subs_abbrs, subs_new = self._get_metabolite_abbrs_and_new(
            reaction, 'substrates', generated_abbrs
        )
        prods_abbrs, prods_new = self._get_metabolite_abbrs_and_new(
            reaction, 'products', generated_abbrs
        )

        # Build VMH formula
        vmh_formula = self._build_vmh_formula(reaction, subs_abbrs, prods_abbrs)

        # Combine new mets
        all_new_mets = subs_new + prods_new
        new_mets_str = ', '.join(all_new_mets) if all_new_mets else ''

        # Format references
        references_str = self._format_references(reaction.references)

        # Extract GPR
        gpr = self._extract_gpr(reaction.gene_info)

        # Format comments
        comments_str = self._format_comments(reaction.comments)

        # Format external links
        ext_links_str = self._format_ext_links(reaction.ext_links)

        return {
            'rxn_abbr': reaction.short_name or '',
            'description': reaction.description or '',
            'vmh_formula': vmh_formula,
            'new_mets': new_mets_str,
            'references': references_str,
            'GPR': gpr,
            'confidence_score': reaction.confidence_score or '',
            'comments': comments_str,
            'ext_links': ext_links_str,
        }

    def _get_metabolite_abbrs_and_new(self, reaction, side, generated_abbrs):
        """
        Get abbreviations for metabolites on one side of the reaction.
        Also identify which are new (not in VMH).
        
        Args:
            reaction: Reaction object
            side: 'substrates' or 'products'
            generated_abbrs: dict of pre-generated abbreviations
            
        Returns:
            tuple: (list of abbreviations, list of new met abbreviations)
        """
        if side == 'substrates':
            mets = self._safe_json_loads(reaction.substrates, [])
            types = self._safe_json_loads(reaction.substrates_types, [])
            names = self._safe_json_loads(reaction.substrates_names, [])
            found = self._safe_json_loads(reaction.subs_found, [])
        else:
            mets = self._safe_json_loads(reaction.products, [])
            types = self._safe_json_loads(reaction.products_types, [])
            names = self._safe_json_loads(reaction.products_names, [])
            found = self._safe_json_loads(reaction.prod_found, [])

        abbrs = []
        new_mets = []

        for i, (met, met_type, met_name) in enumerate(zip(mets, types, names)):
            in_vmh = found[i] if i < len(found) else False

            if met_type.lower() == 'vmh':
                # VMH type: the identifier IS the abbreviation
                abbr = met
            elif met_type == 'Saved':
                # Saved metabolite: look up in database
                abbr = self._get_saved_metabolite_abbr(met, met_name, generated_abbrs)
            else:
                # Other types: use generated abbreviation
                key = (met, met_type, met_name)
                abbr = generated_abbrs.get(key, self._sanitize_name_as_abbr(met_name))

            abbrs.append(abbr)

            # Track if this is a new metabolite (not in VMH)
            if not in_vmh:
                new_mets.append(abbr)

        return abbrs, new_mets

    def _get_saved_metabolite_abbr(self, met_id, met_name, generated_abbrs):
        """
        Get abbreviation for a saved metabolite.
        """
        try:
            saved_met = SavedMetabolite.objects.get(id=int(met_id))
            if saved_met.vmh_abbr:
                return saved_met.vmh_abbr
        except (SavedMetabolite.DoesNotExist, ValueError):
            pass
        
        # Fallback to generated or sanitized name
        for key, abbr in generated_abbrs.items():
            if key[2] == met_name:  # Match by name
                return abbr
        
        return self._sanitize_name_as_abbr(met_name)

    def _build_vmh_formula(self, reaction, subs_abbrs, prods_abbrs):
        """
        Construct the VMH formula string.
        
        Format: stoich0 abbr0[comp0] + stoich1 abbr1[comp1] + ... --> ...
        """
        subs_stoich = self._safe_json_loads(reaction.subs_sch, [])
        prods_stoich = self._safe_json_loads(reaction.prods_sch, [])
        subs_comps = self._safe_json_loads(reaction.subs_comps, [])
        prods_comps = self._safe_json_loads(reaction.prods_comps, [])

        # Build substrate side
        subs_parts = []
        for i, abbr in enumerate(subs_abbrs):
            stoich = subs_stoich[i] if i < len(subs_stoich) else 1
            comp = subs_comps[i] if i < len(subs_comps) else 'c'
            subs_parts.append(self._format_metabolite_term(stoich, abbr, comp))

        # Build product side
        prods_parts = []
        for i, abbr in enumerate(prods_abbrs):
            stoich = prods_stoich[i] if i < len(prods_stoich) else 1
            comp = prods_comps[i] if i < len(prods_comps) else 'c'
            prods_parts.append(self._format_metabolite_term(stoich, abbr, comp))

        # Determine arrow based on direction
        direction = reaction.direction or 'forward'
        arrow = '-->' if direction.lower() == 'forward' else '<-->'

        # Combine
        lhs = ' + '.join(subs_parts)
        rhs = ' + '.join(prods_parts)

        return f"{lhs} {arrow} {rhs}"

    def _format_metabolite_term(self, stoich, abbr, comp):
        """
        Format a single metabolite term: stoich abbr[comp]
        """
        # Convert stoich to clean string (remove trailing .0 for integers)
        try:
            stoich_val = float(stoich)
            if stoich_val == int(stoich_val):
                stoich_str = str(int(stoich_val))
            else:
                stoich_str = str(stoich_val)
        except (ValueError, TypeError):
            stoich_str = str(stoich) if stoich else '1'

        return f"{stoich_str} {abbr}[{comp}]"

    def _format_references(self, references):
        """
        Format references as: PMID:12345, DOI:10.xxxx, ...
        """
        if not references:
            return ''

        formatted = []
        for ref in references:
            if isinstance(ref, dict):
                ref_type = ref.get('ref_type', '')
                info = ref.get('info', '')
                if info:
                    formatted.append(f"{ref_type}:{info}" if ref_type else info)
            elif isinstance(ref, str):
                formatted.append(ref)

        return ', '.join(formatted)

    def _extract_gpr(self, gene_info):
        """
        Extract GPR rule from gene_info using the same logic as prepare_vmh_update_json_files.
        """
        if not gene_info:
            return ''

        try:
            return merge_gene_infos(gene_info)
        except Exception:
            # Fallback: try to extract manually
            gprs = []
            for item in gene_info:
                if isinstance(item, dict) and 'info' in item:
                    info = item['info']
                    # Parse: take part before semicolon, remove "GPR: " prefix
                    first_part = info.split(';')[0].strip()
                    if first_part.startswith("GPR: "):
                        first_part = first_part[5:]
                    if first_part:
                        gprs.append(first_part)
            
            if len(gprs) == 1:
                return gprs[0]
            elif len(gprs) > 1:
                return ' or '.join(gprs)
            return ''

    def _format_comments(self, comments):
        """
        Format comments as a comma-separated string.
        """
        if not comments:
            return ''

        formatted = []
        for comment in comments:
            if isinstance(comment, dict):
                info = comment.get('info', '')
                if info:
                    formatted.append(info)
            elif isinstance(comment, str):
                formatted.append(comment)

        return ', '.join(formatted)

    def _format_ext_links(self, ext_links):
        """
        Format external links as: source:id, source:id, ...
        """
        if not ext_links:
            return ''

        formatted = []
        for link in ext_links:
            if isinstance(link, dict):
                link_type = link.get('ext_link_type', '')
                info = link.get('info', '')
                if info:
                    formatted.append(f"{link_type}:{info}" if link_type else info)
            elif isinstance(link, str):
                formatted.append(link)

        return ', '.join(formatted)

    def _write_csv(self, output_path, rows):
        """
        Write rows to a CSV file.
        """
        if not rows:
            return

        fieldnames = [
            'rxn_abbr',
            'description', 
            'vmh_formula',
            'new_mets',
            'references',
            'GPR',
            'confidence_score',
            'comments',
            'ext_links',
        ]

        with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    def _safe_json_loads(self, value, default=None):
        """
        Safely load JSON from a string, returning default on failure.
        """
        if default is None:
            default = []
        
        if not value:
            return default
        
        if isinstance(value, (list, dict)):
            return value
        
        try:
            return json.loads(value)
        except (json.JSONDecodeError, TypeError):
            return default

    def _sanitize_name_as_abbr(self, name):
        """
        Create a fallback abbreviation from a metabolite name.
        """
        if not name:
            return 'unknown'
        
        # Take first 20 chars, replace spaces with underscores, remove special chars
        sanitized = name[:20].replace(' ', '_')
        sanitized = ''.join(c for c in sanitized if c.isalnum() or c == '_')
        return sanitized or 'unknown'
