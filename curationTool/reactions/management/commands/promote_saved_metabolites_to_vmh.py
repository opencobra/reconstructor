import json
from collections import defaultdict
from typing import Dict, Optional, Set

import urllib3
from django.core.management.base import BaseCommand
from tqdm import tqdm

from reactions.models import Reaction, SavedMetabolite
from reactions.utils.vmh_api import search_metabolites_new, vmh_metabolite_url

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


class Command(BaseCommand):
    help = (
        "Promote Saved metabolites that now exist in VMH (matched by EXACT "
        "fullName via the VMH API) by converting reaction entries from "
        "type 'Saved' to type 'VMH'."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Preview changes without writing to the database.",
        )
        parser.add_argument(
            "--keep-promoted",
            action="store_true",
            help="Do not delete promoted SavedMetabolite rows after conversion.",
        )

    @staticmethod
    def _safe_json_list(raw_value, default=None):
        if default is None:
            default = []
        if raw_value is None:
            return list(default)
        if isinstance(raw_value, list):
            return raw_value
        try:
            parsed = json.loads(raw_value)
            return parsed if isinstance(parsed, list) else list(default)
        except Exception:
            return list(default)

    @staticmethod
    def _exact_name_row(rows: list, name: str) -> Optional[Dict]:
        target = (name or "").strip()
        if not target:
            return None
        for row in rows:
            if (row.get("fullName") or "").strip() == target:
                return row
        return None

    def _resolve_saved_metabolite(self, saved_met: SavedMetabolite) -> Optional[Dict[str, str]]:
        name = (saved_met.name or "").strip()
        if not name:
            return None

        try:
            rows = search_metabolites_new({"fullName": name})
        except Exception as exc:
            self.stderr.write(f"VMH search failed for '{name}' (id={saved_met.id}): {exc}")
            return None

        row = self._exact_name_row(rows, name)
        if not row:
            return None

        abbr = row.get("abbreviation") or ""
        if not abbr:
            return None

        return {
            "abbr": abbr,
            "name": row.get("fullName") or name,
            "url": vmh_metabolite_url(abbr),
        }

    def _remaining_saved_ids(self) -> Set[str]:
        remaining: Set[str] = set()
        reactions = Reaction.objects.only("substrates", "products", "substrates_types", "products_types")
        for reaction in reactions.iterator():
            substrates = self._safe_json_list(reaction.substrates)
            products = self._safe_json_list(reaction.products)
            subs_types = self._safe_json_list(reaction.substrates_types)
            prod_types = self._safe_json_list(reaction.products_types)

            for met_id, met_type in zip(substrates, subs_types):
                if met_type == "Saved":
                    remaining.add(str(met_id))
            for met_id, met_type in zip(products, prod_types):
                if met_type == "Saved":
                    remaining.add(str(met_id))
        return remaining

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        keep_promoted = options["keep_promoted"]

        saved_lookup = {
            str(saved_met.id): saved_met
            for saved_met in SavedMetabolite.objects.select_related("owner").only(
                "id", "name", "vmh_abbr", "owner__id", "owner__name"
            )
        }
        promotion_cache: Dict[str, Optional[Dict[str, str]]] = {}
        promoted_saved_ids: Set[str] = set()
        user_matched: Dict[str, Set[str]] = defaultdict(set)
        user_unmatched: Dict[str, list] = defaultdict(list)

        reactions_updated = 0
        replacements_count = 0
        vmh_matches = 0

        reactions = Reaction.objects.only(
            "id",
            "substrates",
            "products",
            "substrates_types",
            "products_types",
            "substrates_names",
            "products_names",
            "subs_found",
            "prod_found",
            "subs_miriams",
            "prod_miriams",
        )

        for reaction in tqdm(reactions.iterator(), total=reactions.count(), desc="Processing reactions"):
            substrates = self._safe_json_list(reaction.substrates)
            products = self._safe_json_list(reaction.products)
            subs_types = self._safe_json_list(reaction.substrates_types)
            prod_types = self._safe_json_list(reaction.products_types)
            subs_names = self._safe_json_list(reaction.substrates_names)
            prod_names = self._safe_json_list(reaction.products_names)
            subs_found = self._safe_json_list(reaction.subs_found, [False] * len(substrates))
            prod_found = self._safe_json_list(reaction.prod_found, [False] * len(products))
            subs_miriams = self._safe_json_list(reaction.subs_miriams, [""] * len(substrates))
            prod_miriams = self._safe_json_list(reaction.prod_miriams, [""] * len(products))

            while len(subs_found) < len(substrates):
                subs_found.append(False)
            while len(prod_found) < len(products):
                prod_found.append(False)
            while len(subs_miriams) < len(substrates):
                subs_miriams.append("")
            while len(prod_miriams) < len(products):
                prod_miriams.append("")
            while len(subs_names) < len(substrates):
                subs_names.append("")
            while len(prod_names) < len(products):
                prod_names.append("")

            reaction_changed = False

            def maybe_promote(met_id, met_type):
                nonlocal vmh_matches
                if met_type != "Saved":
                    return None
                saved_id = str(met_id)
                if saved_id in promotion_cache:
                    return promotion_cache[saved_id]

                saved_met = saved_lookup.get(saved_id)
                if not saved_met:
                    promotion_cache[saved_id] = None
                    return None

                owner = saved_met.owner
                user_label = f"{owner.name} (id={owner.id})" if owner else "Unknown"

                promoted = self._resolve_saved_metabolite(saved_met)
                promotion_cache[saved_id] = promoted
                if promoted:
                    vmh_matches += 1
                    user_matched[user_label].add(saved_id)
                    if not dry_run and saved_met.vmh_abbr != promoted["abbr"]:
                        saved_met.vmh_abbr = promoted["abbr"]
                        saved_met.save(update_fields=["vmh_abbr"])
                else:
                    user_unmatched[user_label].append((saved_met.name or "", saved_id))
                return promoted

            for idx, (met_id, met_type) in enumerate(zip(substrates, subs_types)):
                promoted = maybe_promote(met_id, met_type)
                if not promoted:
                    continue
                substrates[idx] = promoted["abbr"]
                subs_types[idx] = "VMH"
                subs_found[idx] = True
                subs_miriams[idx] = promoted["url"]
                if not subs_names[idx]:
                    subs_names[idx] = promoted["name"]
                reaction_changed = True
                replacements_count += 1
                promoted_saved_ids.add(str(met_id))

            for idx, (met_id, met_type) in enumerate(zip(products, prod_types)):
                promoted = maybe_promote(met_id, met_type)
                if not promoted:
                    continue
                products[idx] = promoted["abbr"]
                prod_types[idx] = "VMH"
                prod_found[idx] = True
                prod_miriams[idx] = promoted["url"]
                if not prod_names[idx]:
                    prod_names[idx] = promoted["name"]
                reaction_changed = True
                replacements_count += 1
                promoted_saved_ids.add(str(met_id))

            if reaction_changed:
                reactions_updated += 1
                if not dry_run:
                    reaction.substrates = json.dumps(substrates)
                    reaction.products = json.dumps(products)
                    reaction.substrates_types = json.dumps(subs_types)
                    reaction.products_types = json.dumps(prod_types)
                    reaction.substrates_names = json.dumps(subs_names)
                    reaction.products_names = json.dumps(prod_names)
                    reaction.subs_found = json.dumps(subs_found)
                    reaction.prod_found = json.dumps(prod_found)
                    reaction.subs_miriams = json.dumps(subs_miriams)
                    reaction.prod_miriams = json.dumps(prod_miriams)
                    reaction.save(
                        update_fields=[
                            "substrates",
                            "products",
                            "substrates_types",
                            "products_types",
                            "substrates_names",
                            "products_names",
                            "subs_found",
                            "prod_found",
                            "subs_miriams",
                            "prod_miriams",
                        ]
                    )

        deleted_count = 0
        if promoted_saved_ids and not keep_promoted:
            remaining = self._remaining_saved_ids()
            deletable_ids = [int(saved_id) for saved_id in promoted_saved_ids if saved_id not in remaining]
            if deletable_ids and not dry_run:
                deleted_count, _ = SavedMetabolite.objects.filter(id__in=deletable_ids).delete()
            elif deletable_ids and dry_run:
                deleted_count = len(deletable_ids)

        self.stdout.write(self.style.SUCCESS(
            f"\nPromotion complete (dry_run={dry_run}): "
            f"vmh_matches={vmh_matches}, reactions_updated={reactions_updated}, "
            f"replacements={replacements_count}, deleted_saved_metabolites={deleted_count}"
        ))
        self.stdout.write("\n--- Per-user summary ---")
        all_users = set(user_matched) | set(user_unmatched)
        for user_label in sorted(all_users):
            matched_count = len(user_matched.get(user_label, set()))
            unmatched = user_unmatched.get(user_label, [])
            self.stdout.write(f"\n  {user_label}: {matched_count} matched, {len(unmatched)} unmatched")
            if unmatched:
                self.stdout.write("    Unmatched saved metabolites:")
                for name, saved_id in unmatched:
                    self.stdout.write(f"      - id={saved_id}  {name}")
