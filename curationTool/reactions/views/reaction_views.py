"""
This module provides Django views for handling reactions and related data.
It includes functionalities for:
- Inputting and retrieving reaction details
- Saving and managing user reactions
- Handling reaction-related metadata (e.g., gene info, references, comments)
- Checking reaction uniqueness and cloning reactions
- Updating confidence scores and other attributes

Dependencies:
- Django
- RDKit
- JSON handling

"""
import json
import re
import os

from rdkit import RDLogger, Chem
from rdkit.Chem import MolToInchiKey, AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula # pylint: disable=import-error
from rdkit.Chem.Descriptors import MolWt # pylint: disable=import-error
from django.shortcuts import render, redirect
from django.urls import reverse
from django.http import JsonResponse
from django.core import serializers
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_POST
from django.conf import settings
from django.views.decorators.http import require_GET
from django.shortcuts import get_object_or_404
from django.db.models import Q

from reactions.views.user_views import validate_user_ID
from reactions.utils.utils import safe_json_loads
from reactions.reaction_info import get_reaction_info
from reactions.utils.process_strings import construct_reaction_rxnfile
from reactions.utils.get_mol_info import get_mol_info
from reactions.utils.search_vmh import search_metabolites_vmh, check_reaction_vmh
from reactions.utils.vmh_api import vmh_public_base_url
from reactions.utils.to_mol import any_to_mol
from reactions.utils.utils import get_fields
from reactions.utils.RDT import RDT
from reactions.utils.to_mol import smiles_with_explicit_hydrogens
from reactions.forms import ReactionForm
from reactions.models import User, CreatedReaction, Reaction, ReactionsAddedVMH, SavedMetabolite
# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')


def _normalize_direction(direction):
    """Normalize direction values coming from different clients/sources."""
    value = (direction or '').strip().lower()
    if value in {'forward', '->', '=>', '='}:
        return 'forward'
    if value in {'bidirectional', 'reversible', '<=>', '<->', '<-->', '↔', '⇌'}:
        return 'bidirectional'
    if value in {'reverse', 'backward', '<=', '<-'}:
        # The app currently supports forward/bidirectional only.
        return 'bidirectional'
    return 'forward'


def input_reaction(request):
    """
    Process and store a reaction submitted by a user.

    Process:
        - Parses and validates user input.
        - Converts metabolite data to molecular representations.
        - Generates a reaction signature for uniqueness.
        - Stores reaction data and metadata.

    Parameters:
        request (HttpRequest): The HTTP request containing reaction data.

    Returns:
        JsonResponse:
            - Success: Contains processed reaction details.
            - Error: If user input is invalid or processing fails.
    """
    if request.method != 'POST':
        return render(request, 'reactions/home_page.html', {'form': ReactionForm()})

    # Action to perform (either 'create' or 'edit')
    action = request.POST.get('action')
    form = ReactionForm(request.POST, request.FILES)
    user_id = request.POST.get('userID')
    user = User.objects.get(pk=user_id) if (user_id and user_id != 'null') else None

    if action == 'edit':
        reaction_id = request.POST.get('reaction_id')
        try:
            reaction = Reaction.objects.get(id=reaction_id)
        except Reaction.DoesNotExist:
            return JsonResponse(
                {'message': 'Reaction not found', 'status': 'error'}
            )
    else:
        # Create a new Reaction object if action is not 'edit'
        reaction = Reaction()

    # Get multiple substrates, products, and their stoichiometry as lists
    substrates_list = request.POST.getlist('substrates')
    products_list = request.POST.getlist('products')
    names_dict = request.POST.get('nameData')
    organs = request.POST.get('organs')

    names_dict = json.loads(names_dict)
    substrates_names = []
    products_names = []

    for key, value in names_dict.items():
        if 'substrate' in key:
            substrates_names.append(value)
        elif 'product' in key:
            products_names.append(value)
        else:
            raise ValueError(f"Invalid key: {key}")

    # Stoichiometry for substrates
    subs_sch = request.POST.getlist('subs_sch')
    prod_sch = request.POST.getlist(
        'prod_sch')  # Stoichiometry for products
    subs_comp = request.POST.getlist(
        'subs_comps')  # Compartments for substrates
    prod_comp = request.POST.getlist('prod_comps')
    substrates_types = request.POST.getlist('substrates_type')
    products_types = request.POST.getlist('products_type')
    if ('Saved' in substrates_types or 'Saved' in products_types) and not user:
        return JsonResponse(
            {'status': 'error', 'message': 'Cannot use saved metabolites without signing in.'})
    direction = _normalize_direction(request.POST.get('direction'))
    subs_sch = [int(s) for s in subs_sch]
    prod_sch = [int(s) for s in prod_sch]

    subs_mols, subs_errors, _ = any_to_mol(
        substrates_list, substrates_types, request, side='substrates')
    prod_mols, prod_errors, _ = any_to_mol(
        products_list, products_types, request, side='products')
    subsystem = request.POST.get('subsystem')
    all_errors = subs_errors + prod_errors
    if any(elem is not None for elem in all_errors):
        error_message = "\n".join(
            [error for error in all_errors if error is not None])
        if request.headers.get('X-Requested-With') == 'XMLHttpRequest':
            return JsonResponse(
                {'status': 'error', 'message': error_message})
        # Return error message in context for non-AJAX requests
        context = {'form': form, 'error_message': error_message}
        return render(request, 'reactions/home_page.html', context)
    mol_data = get_mol_info(subs_mols + prod_mols)
    smiles = mol_data['smiles']
    inchis = mol_data['inchis']
    inchi_keys = mol_data['inchi_keys']
    mol_weights = mol_data['mol_weights']
    metabolite_formulas = mol_data['metabolite_formulas']
    metabolite_charges = mol_data['metabolite_charges']
    metabolite_mol_file_strings = mol_data['metabolite_mol_file_strings']
    stereo_counts = mol_data['stereo_counts']
    stereo_locations_list = mol_data['stereo_locations_list']
    metabolite_names = substrates_names + products_names
    reaction_rxn_file = construct_reaction_rxnfile(
        subs_mols, subs_sch, prod_mols, prod_sch, substrates_names, products_names)
    subs_inchi_keys = [
        (MolToInchiKey(mol), int(stoich))
        for mol, stoich in zip(subs_mols, subs_sch)
    ]
    prods_inchi_keys = [
        (MolToInchiKey(mol), int(stoich))
        for mol, stoich in zip(prod_mols, prod_sch)
    ]
    # Sort to make it order-independent
    subs_inchi_keys.sort()
    prods_inchi_keys.sort()

    # Create a unique reaction signature (store as JSON string)
    reaction_signature = json.dumps({
        "substrates": subs_inchi_keys,
        "products": prods_inchi_keys
    }, sort_keys=True)

    reaction.save()

    # Skip atom mapping if any product fields are empty or only one
    skip_atom_mapping = request.POST.get('skipAtomMapping') == 'true' or (
        len(substrates_list) == 1 and len(products_list) == 0)
    if skip_atom_mapping:
        response_data = {'visualizations': [
            '/images/atom_mapping_skip.png']}
        reaction_info = get_reaction_info(reaction_rxn_file, direction)
    else:
        # Use relative paths for RDT, but ensure directories exist with absolute paths
        png_rel_path = f'media/images/visual{reaction.id}.png'
        rxn_rel_path = f'media/rxn_files/rxn{reaction.id}.rxn'
        os.makedirs(os.path.join(settings.MEDIA_ROOT, 'images'), exist_ok=True)
        os.makedirs(os.path.join(settings.MEDIA_ROOT, 'rxn_files'), exist_ok=True)
        
        response_data = RDT(
            reaction_rxn_file,
            destination_path_png=png_rel_path,
            destination_path_rxn=rxn_rel_path)
        reaction_info = get_reaction_info(
            rxn_rel_path, direction
        )
    balanced_count = reaction_info['balanced_count']
    subs_atoms = reaction_info['subs_atoms']
    prods_atoms = reaction_info['prods_atoms']
    balanced_charge = reaction_info['balanced_charge']
    subs_charge = reaction_info['subs_charge']
    prods_charge = reaction_info['prods_charge']
    molc_formula = reaction_info['molc_formula']
    symb_to_name = reaction_info['symb_to_name']

    vmh_found = check_reaction_vmh(
        substrates_list,
        products_list,
        subs_sch,
        prod_sch,
        substrates_types,
        products_types,
        subs_mols,
        prod_mols,
        direction,
        subsystem,
        subs_comp,
        prod_comp)

    if 'error' in response_data:
        context = {'form': form, 'error_message': response_data['error']}
        return render(request, 'reactions/home_page.html', context)

    subs_found, subs_miriams = search_metabolites_vmh(
        substrates_list, substrates_types, request, side='substrates')
    prod_found, prod_miriams = search_metabolites_vmh(
        products_list, products_types, request, side='products')
    substrates_list = get_fields(
        request,
        substrates_list,
        substrates_types,
        settings.MEDIA_ROOT,
        settings.MEDIA_URL,
        side='substrates')
    products_list = get_fields(
        request,
        products_list,
        products_types,
        settings.MEDIA_ROOT,
        settings.MEDIA_URL,
        side='products')
    # Save non-VMH metabolites
    if user:
        def save_metabolites(metabolites, found_list, names, types, user,substrates=True):
            for i, (mol, found) in enumerate(zip(metabolites, found_list)):
                if not found and user:
                    mol_smiles = Chem.MolToSmiles(mol, allHsExplicit=True)
                    mol_smiles = smiles_with_explicit_hydrogens(mol_smiles)
                    mol = Chem.MolFromSmiles(mol_smiles)
                    mol = Chem.AddHs(mol, explicitOnly=True)
                    inchi_key = MolToInchiKey(mol)
                    smiles = Chem.MolToSmiles(mol)
                    inchi = Chem.MolToInchi(mol)
                    mol_w = MolWt(mol)
                    mol_formula = CalcMolFormula(mol)
                    try:
                        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                            raise ValueError("Embedding failed for molecule.")
                        AllChem.UFFOptimizeMolecule(mol)
                        mol_file = Chem.MolToMolBlock(mol)
                    except Exception as e:
                        mol_file = Chem.MolToMolBlock(mol)
                    # Check if the metabolite already exists for the user (owner or shared)
                    if not SavedMetabolite.objects.filter(
                        Q(owner=user) | Q(shared_with=user),
                        inchi_key=inchi_key
                    ).exists():
                        if SavedMetabolite.objects.filter(
                            Q(owner=user) | Q(shared_with=user),
                            name=names[i].strip()
                        ).exists():
                            return JsonResponse(
                                {
                                    'status': 'error',
                                    'message': (
                                        f"Metabolite with name {names[i]} "
                                        "already exists in your saved metabolites. "
                                        "Choose a different name or "
                                        "use the existing metabolite."
                                    )
                                }
                            )
                        # Create a new metabolite with the user as the owner
                        obj = SavedMetabolite.objects.create(
                            owner=user,
                            name=names[i].strip(),
                            inchi_key=inchi_key,
                            inchi=inchi,
                            smiles=smiles,
                            mol_w=mol_w,
                            mol_formula=mol_formula,
                            mol_file=mol_file,
                            source_type=types[i],
                            original_identifier=types[i]
                        )
                        obj.save()
                    else:
                        # Get the existing metabolite
                        obj = SavedMetabolite.objects.filter(
                        Q(owner=user) | Q(shared_with=user),
                        inchi_key=inchi_key).first()

                    list_to_update = substrates_list if substrates else products_list
                    type_list = substrates_types if substrates else products_types
                    type_list[i] = 'Saved'
                    list_to_update[i] = obj.pk
            return None
        response = save_metabolites(
            subs_mols, subs_found, substrates_names, substrates_types, user, substrates=True
        )
        if response:
            return response

        response = save_metabolites(
            prod_mols, prod_found, products_names, products_types, user, substrates=False
        )
        if response:
            return response


    # Assign the values directly to the reaction instance
    reaction.Organs = json.dumps(organs)
    reaction.subs_sch = json.dumps(subs_sch)
    reaction.prods_sch = json.dumps(prod_sch)
    reaction.substrates_types = json.dumps(substrates_types)
    reaction.products_types = json.dumps(products_types)
    reaction.subs_comps = json.dumps(subs_comp)
    reaction.prods_comps = json.dumps(prod_comp)
    reaction.substrates_names = json.dumps(substrates_names)
    reaction.products_names = json.dumps(products_names)
    reaction.substrates = json.dumps(substrates_list)
    reaction.products = json.dumps(products_list)
    reaction.direction = direction
    reaction.subsystem = subsystem
    reaction.visualization = json.dumps(response_data['visualizations'])
    reaction.molc_formula = json.dumps([molc_formula])
    reaction.balanced_count = json.dumps([balanced_count])
    reaction.balanced_charge = json.dumps([balanced_charge])
    reaction.subs_atoms = json.dumps([subs_atoms])
    reaction.prods_atoms = json.dumps([prods_atoms])
    reaction.subs_charge = json.dumps([subs_charge])
    reaction.prods_charge = json.dumps([prods_charge])
    reaction.symb_to_name = json.dumps([symb_to_name])
    reaction.subs_found = json.dumps(subs_found)
    reaction.subs_miriams = json.dumps(subs_miriams)
    reaction.prod_found = json.dumps(prod_found)
    reaction.prod_miriams = json.dumps(prod_miriams)
    reaction.vmh_found = vmh_found['found']
    reaction.vmh_found_similar = vmh_found['similar']
    reaction.vmh_url = json.dumps(
        vmh_found['url']) if vmh_found['found'] else None
    reaction.vmh_formula = json.dumps(
        vmh_found['formula']) if vmh_found['found'] else None
    reaction.metabolite_names = json.dumps(metabolite_names)
    reaction.metabolite_formulas = json.dumps(metabolite_formulas)
    reaction.metabolite_smiles = json.dumps(smiles)
    reaction.metabolite_inchis = json.dumps(inchis)
    reaction.metabolite_inchi_keys = json.dumps(inchi_keys)
    reaction.metabolite_mol_weights = json.dumps(mol_weights)
    reaction.metabolite_charges = json.dumps(metabolite_charges)
    reaction.metabolite_mol_file_strings = json.dumps(
        metabolite_mol_file_strings)
    reaction.stereo_counts = json.dumps(stereo_counts)
    reaction.stereo_locations_list = json.dumps(stereo_locations_list)
    reaction.reaction_signature = reaction_signature
    reaction.save()
    if request.headers.get('X-Requested-With') == 'XMLHttpRequest':
        data = {
            'visualization': json.loads(reaction.visualization),
            'molc_formula': json.loads(reaction.molc_formula),
            'balanced_count': json.loads(reaction.balanced_count),
            'balanced_charge': json.loads(reaction.balanced_charge),
            'subs_atoms': json.loads(reaction.subs_atoms),
            'prods_atoms': json.loads(reaction.prods_atoms),
            'subs_charge': json.loads(reaction.subs_charge),
            'prods_charge': json.loads(reaction.prods_charge),
            'symb_to_name': json.loads(reaction.symb_to_name),
            'subs_found': json.loads(reaction.subs_found),
            'subs_miriams': json.loads(reaction.subs_miriams),
            'prod_found': json.loads(reaction.prod_found),
            'prod_miriams': json.loads(reaction.prod_miriams),
            'reaction_id': reaction.id,
            'metabolite_names': json.loads(reaction.metabolite_names),
            'metabolite_formulas': json.loads(reaction.metabolite_formulas),
            'metabolite_charges': json.loads(reaction.metabolite_charges),
            'metabolite_mol_file_strings': json.loads(reaction.metabolite_mol_file_strings),
            'sterio_counts': json.loads(reaction.stereo_counts),
            'sterio_locations_list': json.loads(reaction.stereo_locations_list),
        }
        if vmh_found['found']:
            data['vmh_found'] = vmh_found['found']
            data['vmh_found_similar'] = vmh_found['similar']
            data['vmh_url'] = vmh_found['url']
            data['vmh_formula'] = vmh_found['formula']
        data['status'] = 'success'
        return JsonResponse(data)


def get_reaction(request, reaction_id):
    """
    Retrieve details of a reaction by its ID.

    Process:
        - Fetches the reaction from the database.
        - Converts stored molecular data if necessary.
        - Formats the response with reaction metadata.

    Parameters:
        request (HttpRequest): The HTTP request object.
        reaction_id (int): The ID of the reaction.

    Returns:
        JsonResponse:
            - Success: Contains reaction details.
            - Error: If the reaction ID is invalid or not found.
    """
    try:
        reaction = Reaction.objects.get(pk=reaction_id)
        if not reaction.metabolite_smiles:
            metabolite_smiles = []
            metabolite_inchis = []
            metabolite_inchi_keys = []
            metabolite_mol_weights = []
            for mol_file in json.loads(reaction.metabolite_mol_file_strings):
                mol = Chem.MolFromMolBlock(mol_file)
                mol = Chem.AddHs(mol, explicitOnly=True)
                smiles = Chem.MolToSmiles(mol)
                inchi = Chem.MolToInchi(mol)
                inchi_key = MolToInchiKey(mol)
                metabolite_smiles.append(smiles)
                metabolite_inchis.append(inchi)
                metabolite_inchi_keys.append(inchi_key)
                metabolite_mol_weights.append(MolWt(mol))
            reaction.metabolite_smiles = json.dumps(metabolite_smiles)
            reaction.metabolite_inchis = json.dumps(metabolite_inchis)
            reaction.metabolite_inchi_keys = json.dumps(metabolite_inchi_keys)
            reaction.metabolite_mol_weights = json.dumps(metabolite_mol_weights)
            reaction.save()

        reaction_data = {
            'Organs': reaction.Organs,
            'reaction_id': reaction.id,
            'short_name': reaction.short_name,
            'description': reaction.description,
            'substrates': safe_json_loads(reaction.substrates),
            'products': safe_json_loads(reaction.products),
            'substrates_names': safe_json_loads(reaction.substrates_names),
            'products_names': safe_json_loads(reaction.products_names),
            'direction': _normalize_direction(reaction.direction),
            'subsystem': reaction.subsystem,
            'subs_comps': safe_json_loads(reaction.subs_comps),
            'prods_comps': safe_json_loads(reaction.prods_comps),
            'visualization': safe_json_loads(reaction.visualization),
            'rxn_formula': safe_json_loads(reaction.rxn_formula),
            'molc_formula': safe_json_loads(reaction.molc_formula),
            'balanced_count': safe_json_loads(reaction.balanced_count),
            'balanced_charge': safe_json_loads(reaction.balanced_charge),
            'subs_sch': safe_json_loads(reaction.subs_sch),
            'prod_sch': safe_json_loads(reaction.prods_sch),
            'subs_types': safe_json_loads(reaction.substrates_types),
            'prods_types': safe_json_loads(reaction.products_types),
            'subs_atoms': safe_json_loads(reaction.subs_atoms),
            'prods_atoms': safe_json_loads(reaction.prods_atoms),
            'subs_charge': safe_json_loads(reaction.subs_charge),
            'prods_charge': safe_json_loads(reaction.prods_charge),
            'symb_to_name': safe_json_loads(reaction.symb_to_name),
            'subs_found': safe_json_loads(reaction.subs_found),
            'subs_miriams': safe_json_loads(reaction.subs_miriams),
            'prod_found': safe_json_loads(reaction.prod_found),
            'prod_miriams': safe_json_loads(reaction.prod_miriams),
            'vmh_found': reaction.vmh_found,
            'vmh_found_similar': reaction.vmh_found_similar,
            'vmh_url': reaction.vmh_url,
            'vmh_formula': reaction.vmh_formula,
            'metabolite_names': safe_json_loads(reaction.metabolite_names),
            'metabolite_formulas': safe_json_loads(reaction.metabolite_formulas),
            'metabolite_charges': safe_json_loads(reaction.metabolite_charges),
            'metabolite_mol_file_strings': safe_json_loads(reaction.metabolite_mol_file_strings),
            'stereo_counts': safe_json_loads(reaction.stereo_counts),
            'stereo_locations_list': safe_json_loads(reaction.stereo_locations_list),
            'metabolite_smiles': safe_json_loads(reaction.metabolite_smiles),
            'metabolite_inchis': safe_json_loads(reaction.metabolite_inchis),
            'metabolite_inchi_keys': safe_json_loads(reaction.metabolite_inchi_keys),
            'metabolite_mol_weights': safe_json_loads(reaction.metabolite_mol_weights),
        }
        return JsonResponse(reaction_data)

    except Reaction.DoesNotExist:
        return JsonResponse({'error': 'Reaction not found'}, status=404)

    except Exception as e:
        return JsonResponse(
            {'error': f'An unexpected error occurred ({str(e)})'},
            status=500
        )


@csrf_exempt
def add_info_to_reaction(request):
    """
    Add metadata (e.g., gene info, external links, references, comments) to a reaction.

    Process:
        - Validates user input.
        - Associates the provided metadata with the reaction.
        - Saves the updated reaction information.

    Parameters:
        request (HttpRequest): The HTTP request containing metadata and reaction ID.

    Returns:
        JsonResponse:
            - Success: Confirmation of metadata addition.
            - Error: If user input is missing or the reaction is not found.
    """
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('userID')
            info_type = data.get('infoType')
            info_text = data.get('infoText')
            ext_link_type = data.get('extLinkType', '')
            ref_type = data.get('refType', '')
            reaction_id = data.get('reactionId')

            if not all([user_id, info_type, info_text]) or (
                    info_type != 'Gene Info' and not reaction_id):
                return JsonResponse({'status': 'error',
                                     'message': 'All fields are required.',
                                     'info_type': info_type},
                                    status=400)

            user = User.objects.get(pk=user_id)
            username = user.name

            if info_type == 'Gene Info' and not reaction_id:
                # Store gene info in session if no reaction_id is provided
                if 'gene_info' not in request.session:
                    request.session['gene_info'] = []

                info_data = {'info': info_text, 'user_name': username}
                request.session['gene_info'].append(info_data)
                request.session.modified = True

                return JsonResponse({'status': 'success',
                                     'message': 'Gene information added to session.',
                                     'info_type': info_type})

            # Fetch the reaction if reaction_id is provided
            reaction = Reaction.objects.get(id=reaction_id)

            if info_type == 'Reference':
                if reaction.references is None:
                    reaction.references = []
                info_data_template = {
                    'user_name': username, 'ref_type': ref_type}

                # Determine the delimiter based on the reference type
                if 'PMID' in ref_type:
                    ref_list = [ref.strip() for ref in info_text.split(';')]
                    existing_refs = {
                        ref['info'] for ref in reaction.references if 'PMID' in ref['info']}
                elif 'DOI' in ref_type:
                    ref_list = [ref.strip() for ref in info_text.split(',')]
                    existing_refs = {
                        ref['info'] for ref in reaction.references if 'DOI' in ref['info']}
                else:
                    ref_list = [info_text.strip()]
                    existing_refs = {ref['info']
                                     for ref in reaction.references}

                for ref in ref_list:
                    if ref not in existing_refs:
                        info_data = info_data_template.copy()
                        info_data['info'] = ref
                        reaction.references.append(info_data)
                        existing_refs.add(ref)

            elif info_type == 'External Link':
                if reaction.ext_links is None:
                    reaction.ext_links = []
                info_data = {
                    'info': info_text,
                    'user_name': username,
                    'ext_link_type': ext_link_type}
                reaction.ext_links.append(info_data)

            elif info_type == 'Gene Info':
                if reaction.gene_info is None:
                    reaction.gene_info = []
                info_data = {'info': info_text, 'user_name': username}
                reaction.gene_info.append(info_data)

            elif info_type == 'Comment':
                if reaction.comments is None:
                    reaction.comments = []
                info_data = {'info': info_text, 'user_name': username}
                reaction.comments.append(info_data)

            reaction.save()
            return JsonResponse({'status': 'success',
                                 'message': 'Information added successfully.',
                                 'reaction_id': reaction_id,
                                 'info_type': info_type})

        except User.DoesNotExist:
            return JsonResponse({'status': 'error',
                                 'message': 'Invalid user key.',
                                 'info_type': info_type},
                                status=404)
        except Reaction.DoesNotExist:
            return JsonResponse({'status': 'error',
                                 'message': 'Invalid reaction ID.',
                                 'info_type': info_type},
                                status=404)
        except json.JSONDecodeError:
            return JsonResponse({'status': 'error',
                                 'message': 'Invalid JSON.',
                                 'info_type': info_type},
                                status=400)
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(
                e), 'info_type': info_type}, status=500)
    else:
        return JsonResponse({'status': 'error',
                             'message': 'Invalid request method.',
                             'info_type': info_type},
                            status=405)


def update_gene_info(request):
    """
    Update gene-related information for a reaction.

    Process:
        - Validates user input.
        - Updates the specific gene-related metadata.
        - Saves the modified reaction information.

    Parameters:
        request (HttpRequest): The HTTP request containing gene details and update data.

    Returns:
        JsonResponse:
            - Success: Confirmation of the gene info update.
            - Error: If gene information is not found or update fails.
    """
    try:
        data = json.loads(request.body)
        user_id = data.get('userID')
        reaction_id = data.get('reactionID')
        gene = data.get('gene')
        field_type = data.get('fieldType')
        updated_value = data.get('updatedValue')

        # Debugging statements to trace values

        # Check if any required field is missing or updated_value is empty
        if not all([user_id, reaction_id, gene, field_type]
                   ) or updated_value is None or updated_value.strip() == "":
            return JsonResponse(
                {
                    'status': 'error',
                    'message': 'All fields are required and updated value cannot be empty.'},
                status=400)

        reaction = Reaction.objects.get(id=reaction_id)

        if reaction.gene_info is None:
            return JsonResponse(
                {'status': 'error', 'message': 'No gene information found.'}, status=404)

        gene_info_updated = False
        for gene_info in reaction.gene_info:
            if gene in gene_info['info']:
                if field_type == 'Organs':
                    organ_match = re.search(
                        r'ORGAN\(([^)]+)\)', gene_info['info'])
                    if organ_match:
                        old_organs = organ_match.group(1)
                        new_info = gene_info['info'].replace(
                            f'ORGAN({old_organs})', f'ORGAN({updated_value})')
                        gene_info['info'] = new_info
                        gene_info_updated = True
                elif field_type == 'SubcellularLocations':
                    subcellular_match = re.search(
                        r'SUBCELLULAR\(([^)]+)\)', gene_info['info'])
                    if subcellular_match:
                        old_subcellular = subcellular_match.group(1)
                        new_info = gene_info['info'].replace(
                            f'SUBCELLULAR({old_subcellular})', f'SUBCELLULAR({updated_value})')
                        gene_info['info'] = new_info
                        gene_info_updated = True

        if gene_info_updated:
            reaction.save()
            return JsonResponse(
                {
                    "status": "success",
                    "message": "Gene information updated successfully.",
                }
            )

        return JsonResponse(
            {
                "status": "error",
                "message": "Gene information not found or not updated.",
            },
            status=404,
        )

    except User.DoesNotExist:
        return JsonResponse(
            {
                "status": "error",
                "message": "Invalid user key.",
            },
            status=404,
        )
    except Reaction.DoesNotExist:
        return JsonResponse(
            {
                "status": "error",
                "message": "Invalid reaction ID.",
            },
            status=404,
        )
    except json.JSONDecodeError:
        return JsonResponse(
            {
                "status": "error",
                "message": "Invalid JSON.",
            },
            status=400,
        )
    except Exception as e:
        return JsonResponse(
            {
                "status": "error",
                "message": str(e),
            },
            status=500,
        )

@csrf_exempt
def delete_reaction_info(request):
    """
    Delete specific reaction metadata (e.g., references, comments, external links, gene info).

    Process:
        - Validates reaction ID and metadata entry.
        - Identifies the metadata category (e.g., references, external links, comments).
        - Removes the specified metadata from the reaction.
        - Saves the updated reaction data.

    Parameters:
        request (HttpRequest): The HTTP request containing `reaction_id`, `tab_id`, and `item_to_delete`.

    Returns:
        JsonResponse:
            - Success: Confirmation of metadata deletion.
            - Error: If the reaction or metadata entry is not found.
    """
    if request.method != "POST":
        return JsonResponse(
            {"error": True, "message": "Invalid request method."}, status=400
        )

    try:
        req_body = json.loads(request.body)
        reaction_id = req_body.get("reaction_id")
        tab_id = req_body.get("tab_id")
        item_to_delete = req_body.get("item_to_delete")

        # Fetch the reaction object
        react_obj = Reaction.objects.get(pk=reaction_id)

        # Handle deletion based on the tab_id
        if tab_id == "refs-content":
            for ref in react_obj.references:
                if (
                    ref["info"] == item_to_delete["info"]
                    and ref["ref_type"] == item_to_delete["ref_type"]
                ):
                    react_obj.references.remove(ref)
                    react_obj.save()
                    return JsonResponse(
                        {
                            "status": "success",
                            "message": "Reference deleted successfully.",
                        }
                    )

        elif tab_id == "ext-links-content":
            for ext_link in react_obj.ext_links:
                if (
                    ext_link["info"] == item_to_delete["info"]
                    and ext_link["ext_link_type"] == item_to_delete["ext_link_type"]
                ):
                    react_obj.ext_links.remove(ext_link)
                    react_obj.save()
                    return JsonResponse(
                        {
                            "status": "success",
                            "message": "External link deleted successfully.",
                        }
                    )

        elif tab_id == "comments-content":
            for comment in react_obj.comments:
                if comment["info"] == item_to_delete["info"]:
                    react_obj.comments.remove(comment)
                    react_obj.save()
                    return JsonResponse(
                        {
                            "status": "success",
                            "message": "Comment deleted successfully.",
                        }
                    )

        elif tab_id == "gene-info-content":
            for gene_info in react_obj.gene_info:
                if gene_info["info"] == item_to_delete["info"]:
                    react_obj.gene_info.remove(gene_info)
                    react_obj.save()
                    return JsonResponse(
                        {
                            "status": "success",
                            "message": "Gene info deleted successfully.",
                        }
                    )

        # If no valid tab_id matched
        return JsonResponse(
            {"error": True, "message": "Invalid tab_id or item not found."}, status=400
        )

    except Reaction.DoesNotExist:
        return JsonResponse(
            {"error": True, "message": "Reaction not found."}, status=404
        )

    except Exception as e:
        return JsonResponse(
            {"error": True, "message": str(e)}, status=500
        )

def get_reaction_details(request):
    """
    Retrieve detailed metadata about a reaction.

    Process:
        - Fetches reaction metadata from the database.
        - Includes references, external links, gene info, and comments.
        - Formats and returns the metadata.

    Parameters:
        request (HttpRequest): The HTTP request containing the reaction ID.

    Returns:
        JsonResponse:
            - Success: Contains reaction metadata details.
            - Error: If the reaction is not found.
    """
    try:
        reaction_id = json.loads(request.body)
        reaction = Reaction.objects.get(pk=reaction_id)

        if reaction.references is None:
            reaction.references = []
        if reaction.ext_links is None:
            reaction.ext_links = []
        if reaction.gene_info is None:
            reaction.gene_info = []
        if reaction.comments is None:
            reaction.comments = []

        return JsonResponse({
            'references': reaction.references,
            'external_links': reaction.ext_links,
            'gene_info': reaction.gene_info,
            'comments': reaction.comments,
        })
    except Reaction.DoesNotExist:
        return JsonResponse({'error': 'Reaction not found'}, status=404)


def delete_reaction(request):
    """
    Remove a reaction from a user's saved reactions list.

    Process:
        - Validates user authentication.
        - Ensures the reaction exists in the user's saved list.
        - Removes the reaction from the user's saved reactions.
        - Redirects the user to the saved reactions page.

    Parameters:
        request (HttpRequest): The HTTP request containing `reaction_id` and `userID`.

    Returns:
        HttpResponse:
            - Success: Redirects to the saved reactions page.
            - Error: If the reaction does not exist in the user’s saved list.
    """
    if request.method == 'POST':
        user_id = request.POST.get('userID')
        reaction_id = request.POST.get('reaction_id')
        user = validate_user_ID(user_id)
        if user:
            reaction = Reaction.objects.get(pk=reaction_id)
            user.saved_reactions.remove(reaction)
            saved_reactions_url = reverse('saved_reactions')
            return redirect(saved_reactions_url)
    return render(request, '404.html')


def saved_reactions(request, modal=False):
    """
    Fetch and display all saved reactions for a user.

    Process:
        - Validates user authentication.
        - Fetches reactions saved by the user.
        - Parses and formats reaction details.
        - Renders a template to display the saved reactions.

    Parameters:
        request (HttpRequest): The HTTP request containing the user ID.
        modal (bool, optional): Whether to return a modal view or a full-page response.

    Returns:
        HttpResponse:
            - Success: Rendered template with saved reaction details.
            - Error: If the user is not found.
    """
    userID = request.session.get('userID')
    user = validate_user_ID(userID)
    if user:
        reactions = user.saved_reactions.all().order_by('id')
        reactions_json = serializers.serialize('json', reactions)
        combined_reactions_details = []

        for reaction in reactions:
            # Parse JSON fields
            subs_sch = json.loads(reaction.subs_sch)
            prods_sch = json.loads(reaction.prods_sch)
            subs_comps = json.loads(reaction.subs_comps)
            prods_comps = json.loads(reaction.prods_comps)
            short_name = reaction.short_name
            description = reaction.description  # Get the reaction description

            # Construct details strings
            subs_details = " + ".join(
                f"{float(sch)} {json.loads(reaction.substrates_names)[idx]} [{comp}]" for idx,
                (sch,
                 comp) in enumerate(
                    zip(
                        subs_sch,
                        subs_comps)))
            prods_details = " + ".join(
                f"{float(sch)} {json.loads(reaction.products_names)[idx]} [{comp}]" for idx,
                (sch,
                 comp) in enumerate(
                    zip(
                        prods_sch,
                        prods_comps)))

            # Check if gene_info is None
            if reaction.gene_info:
                if isinstance(reaction.gene_info, str):
                    gene_info_data = json.loads(reaction.gene_info)
                else:
                    gene_info_data = reaction.gene_info

                gene_info_list = []
                for gene in gene_info_data:
                    if 'info' in gene:
                        info = gene['info']
                        gpr_start = info.find('GPR: ')
                        if gpr_start != -1:
                            gpr_end = info.find(';', gpr_start)
                            if gpr_end != -1:
                                gpr_info = info[gpr_start + 5:gpr_end]
                                gene_info_list.append(gpr_info.strip())
                        else:
                            gene_info_list.append(info.strip())
            else:
                gene_info_list = []

            # Get associated flags and their colors
            flags = reaction.flags.all()
            flag_details = [{"name": flag.name_flag,
                             "id": flag.pk,
                             "color": flag.color} for flag in flags]

            combined_reactions_details.append({
                'reaction': reaction,
                'details': {
                    'short_name': short_name,
                    'description': description,  # Include description in details
                    'subs_details': subs_details,
                    'prods_details': prods_details,
                    'molc_formula': reaction.molc_formula,
                    'balanced_count': json.loads(reaction.balanced_count)[0] if reaction.balanced_count else None,
                    'balanced_charge': json.loads(reaction.balanced_charge)[0] if reaction.balanced_charge else None,
                    'subsystem': reaction.subsystem,
                    'direction': _normalize_direction(reaction.direction),
                    'gene_info': gene_info_list,
                    'flags': flag_details,  # Include flag details with name and color
                    'confidence_score': reaction.confidence_score,
                }
            })

        context = {
            'reactions': reactions,
            'reactions_json': reactions_json,
            'userID': userID,
            'user_name': user.name,
            'vmh_base_url': vmh_public_base_url(),
            'combined_reactions_details': combined_reactions_details,
        }

        if modal:
            return render(
                request,
                'reactions/saved_reactions_modal.html',
                context)
        
        return render(request, 'reactions/saved_reactions.html', context)
    
    return render(request, 'reactions/error.html',
                    {'error_message': 'Invalid key'})

def identical_reaction(request):
    """
    Check if the user has already saved an identical reaction.

    Process:
        - Converts metabolites into molecular representations.
        - Generates a reaction signature for comparison.
        - Checks the user’s saved reactions for a matching signature.

    Parameters:
        request (HttpRequest): The HTTP request containing reaction data.

    Returns:
        JsonResponse:
            - Success: Indicates if a matching reaction exists.
            - Error: If metabolite processing fails.
    """
    user_id = request.POST.get('userID')
    if not user_id or user_id == 'null':
        return JsonResponse({'exists': False, 'status': 'success'})

    # Ensure the user exists
    try:
        user = User.objects.get(pk=user_id)
    except User.DoesNotExist:
        return JsonResponse({'exists': False, 'status': 'success'})

    # Gather substrates/products and their stoichiometries from request
    substrates_list = request.POST.getlist('substrates')
    subs_sch_list = request.POST.getlist('subs_sch')
    products_list = request.POST.getlist('products')
    prods_sch_list = request.POST.getlist('prod_sch')

    substrates_types = request.POST.getlist('substrates_type')
    products_types = request.POST.getlist('products_type')

    # Convert stoichiometries to integers
    try:
        subs_sch = [int(s) for s in subs_sch_list]
        prods_sch = [int(s) for s in prods_sch_list]
    except ValueError:
        return JsonResponse({'error': 'Invalid stoichiometry values', 'status': 'error'})

    # Convert metabolites to RDKit Mol objects
    subs_mols, subs_errors, _ = any_to_mol(
        substrates_list, substrates_types, request, side="substrates"
    )
    prod_mols, prod_errors, _ = any_to_mol(
        products_list, products_types, request, side="products"
    )

    if any(subs_errors) or any(prod_errors):
        return JsonResponse({'error': 'Failed to process some metabolites', 'status': 'error'})
    # Generate InChIKey-based reaction signature
    subs_inchi_keys = [
        (MolToInchiKey(mol), int(stoich))
        for mol, stoich in zip(subs_mols, subs_sch)
        if mol
    ]
    prods_inchi_keys = [
        (MolToInchiKey(mol), int(stoich))
        for mol, stoich in zip(prod_mols, prods_sch)
        if mol
    ]
    # Sort InChIKeys and create reaction signature
    subs_inchi_keys.sort()
    prods_inchi_keys.sort()
    reaction_signature = json.dumps({
        "substrates": subs_inchi_keys,
        "products": prods_inchi_keys
    }, sort_keys=True)

    # Look for existing reactions with the same signature (only for this user)
    matching_reactions = user.saved_reactions.filter(reaction_signature=reaction_signature)

    if matching_reactions.exists():
        matches_data = [
            {
                "reaction_id": r.id,
                "reaction_name": r.short_name or f"Reaction {r.id}",
            }
            for r in matching_reactions
        ]
        return JsonResponse(
            {
                "exists": True,
                "matches": matches_data,
                "status": "success",
            }
        )

    return JsonResponse({'exists': False, 'status': 'success'})


@require_POST
def clone_reaction_view(request):
    """
    Create a duplicate of a reaction for a user.

    Process:
        - Fetches the original reaction from the database.
        - Creates a new reaction entry with a unique name.
        - Saves the cloned reaction under the user’s account.

    Parameters:
        request (HttpRequest): The HTTP request containing `reaction_id` and `name`.

    Returns:
        JsonResponse:
            - Success: Confirmation of reaction cloning.
            - Error: If cloning fails.
    """
    try:
        reaction_id = request.POST.get('reaction_id')
        user_id = request.POST.get('userID')
        reaction_name = request.POST.get('name')
        cloned_reaction = Reaction.objects.get(pk=reaction_id)

        cloned_reaction.pk = None  # Set the primary key to None to create a new instance
        cloned_reaction.short_name = reaction_name
        # Save the cloned reaction object to generate a new ID
        cloned_reaction.save()
        user = User.objects.get(pk=user_id)
        user.saved_reactions.add(cloned_reaction)
        user.save()
        return JsonResponse({'status': 'success',
                             'message': 'Reaction cloned successfully'})
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)})


@require_GET
def get_user_reactions_and_vmh(request):
    """
    Retrieve reaction statistics for all users.

    Process:
        - Fetches all users from the database.
        - Counts saved reactions, created reactions, and reactions added to VMH.
        - Formats and returns the user reaction statistics.

    Parameters:
        request (HttpRequest): The HTTP request object.

    Returns:
        JsonResponse:
            - Success: Contains user reaction statistics.
            - Error: If data retrieval fails.
    """
    try:
        users_data = []
        users = User.objects.all()

        for user in users:
            saved_reactions_count = user.saved_reactions.count()
            reactions_added_vmh_count = ReactionsAddedVMH.objects.filter(
                user=user).count()
            created_reactions_count = CreatedReaction.objects.filter(
                user=user).count()

            user_data = {
                'full_name': user.full_name,
                'saved': saved_reactions_count,
                'added': reactions_added_vmh_count,
                'created': created_reactions_count
            }
            users_data.append(user_data)

        return JsonResponse(users_data, safe=False)

    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)


@csrf_exempt
def create_reaction(request):
    """
    Create a new reaction entry and associate it with a user.

    Process:
        - Validates user authentication.
        - Fetches the specified reaction from the database.
        - Links the reaction to the user.

    Parameters:
        request (HttpRequest): The HTTP request containing `user_id` and `reaction_id`.

    Returns:
        JsonResponse:
            - Success: Confirmation of reaction creation.
            - Error: If the user or reaction does not exist.
    """
    if request.method == 'POST':
        data = json.loads(request.body)
        user_id = data.get('user_id')
        reaction_id = data.get('reaction_id')
        # Validate the user ID using the provided function
        user = validate_user_ID(user_id)

        if not user:
            return JsonResponse(
                {'success': False, 'error': 'User does not exist'}, status=400)

        # Fetch the reaction based on reaction_id
        reaction = get_object_or_404(Reaction, id=reaction_id)

        created_reaction = CreatedReaction.objects.create(
            user=user, reaction=reaction)
        return JsonResponse({'success': True,
                             'created_reaction_id': created_reaction.id})

    return JsonResponse(
        {'success': False, 'error': 'Invalid request method'}, status=400)


@csrf_exempt
def reaction_name_exists(request):
    """
    Check if a reaction name is already used by a user.

    Process:
        - Fetches the user’s saved reactions.
        - Checks if the provided reaction name already exists.

    Parameters:
        request (HttpRequest): The HTTP request containing `user_id` and `short_name`.

    Returns:
        JsonResponse:
            - Success: Indicates if the reaction name is already used.
            - Error: If the user is not found.
    """
    if request.method != 'POST':
        return JsonResponse({'error': 'Invalid request method'}, status=405)

    try:
        data = json.loads(request.body)
        user_id = data.get('user_id')
        reaction_name = data.get('short_name')
        user = User.objects.get(pk=user_id)
        saved_reaction_names = list(
            user.saved_reactions.all().values_list(
                'short_name', flat=True))
        is_name_saved = reaction_name in saved_reaction_names
        return JsonResponse({'is_name_saved': is_name_saved})
    except User.DoesNotExist:
        return JsonResponse({'error': 'User does not exist'}, status=404)
    except (KeyError, ValueError):
        return JsonResponse({'error': 'Invalid data'}, status=400)

def get_available_reactions(request):
    """
    Retrieve a list of available reactions for a user.

    Process:
        - Fetches the user's saved reaction IDs.
        - Retrieves all reactions excluding the user’s saved reactions.
        - Returns the list of available reaction IDs.

    Parameters:
        request (HttpRequest): The HTTP request containing `user_id`.

    Returns:
        JsonResponse:
            - Success: List of available reaction IDs.
            - Error: If the user is not found.
    """
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('user_id')

            user = User.objects.get(pk=user_id)
            saved_reaction_ids = list(
                user.saved_reactions.all().values_list(
                    'id', flat=True))
            available_reactions = Reaction.objects.exclude(
                id__in=saved_reaction_ids)
            available_reaction_ids = list(
                available_reactions.values_list(
                    'id', flat=True))
            last_index = saved_reaction_ids[-1] if saved_reaction_ids else None

            return JsonResponse(
                {'available_reaction_ids': available_reaction_ids, 'last_index': last_index})
        except User.DoesNotExist:
            return JsonResponse({'error': 'User does not exist'}, status=404)
        except (KeyError, ValueError):
            return JsonResponse({'error': 'Invalid data'}, status=400)
    return JsonResponse({'error': 'Invalid request method'}, status=405)


@csrf_exempt
def already_saved(request):
    """
    Check if a reaction is already saved by a user.

    Process:
        - Fetches the user’s saved reaction IDs.
        - Checks if the provided reaction ID is in the user’s saved list.

    Parameters:
        request (HttpRequest): The HTTP request containing `user_id` and `reaction_id`.

    Returns:
        JsonResponse:
            - Success: Indicates whether the reaction is saved.
            - Error: If the user is not found.
    """
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('user_id')
            # Convert reaction_id to integer
            reaction_id = int(data.get('reaction_id'))

            user = User.objects.get(pk=user_id)
            saved_reaction_ids = list(
                user.saved_reactions.all().values_list(
                    'id', flat=True))

            is_reaction_saved = reaction_id in saved_reaction_ids

            return JsonResponse({'is_reaction_saved': is_reaction_saved})
        except User.DoesNotExist:
            return JsonResponse({'error': 'User does not exist'}, status=404)
        except (KeyError, ValueError):
            return JsonResponse({'error': 'Invalid data'}, status=400)

    return JsonResponse({'error': 'Invalid request method'}, status=405)


@csrf_exempt
def edit_reaction_info(request):
    """
    Edit the short name and/or description of a reaction.

    Process:
        - Fetches the reaction from the database.
        - Updates the reaction name and/or description.
        - Ensures that the new name is unique among the user’s saved reactions.

    Parameters:
        request (HttpRequest): The HTTP request containing `reaction_id`, `new_name`, and `new_description`.

    Returns:
        JsonResponse:
            - Success: Confirmation of reaction update.
            - Error: If the reaction does not exist or the new name is already used.
    """
    if request.method != "POST":
        return JsonResponse(
            {"error": "Invalid request method."}, status=405
        )

    try:
        data = json.loads(request.body)
        reaction_id = int(data.get("reaction_id"))
        new_name = data.get("new_name")
        new_description = data.get("new_description")
        user_id = data.get("user_id")

        reaction = Reaction.objects.get(pk=reaction_id)

        if new_name is not None:
            new_name = new_name.strip()
            if new_name != reaction.short_name:
                user_reactions = User.objects.get(pk=user_id).saved_reactions.exclude(pk=reaction_id)

                # Only check for duplicates if name is actually changing
                if user_reactions.filter(short_name=new_name).exists():
                    return JsonResponse(
                        {
                            "error": f'Reaction name "{new_name}" already exists.',
                            "original_name": reaction.short_name,
                        },
                        status=400,
                    )

                reaction.short_name = new_name

        if new_description is not None:
            reaction.description = new_description

        reaction.save()

        return JsonResponse({"success": True})

    except Reaction.DoesNotExist:
        return JsonResponse(
            {"error": "Reaction does not exist."}, status=404
        )

    except (KeyError, ValueError, json.JSONDecodeError):
        return JsonResponse(
            {
                "error": "Invalid data.",
                "original_name": reaction.short_name,
            },
            status=400,
        )

@csrf_exempt
def update_confidence_score(request):
    """
    Update the confidence score of one or more reactions.

    Process:
        - Validates reaction IDs and confidence score values.
        - Updates the confidence score for the specified reactions in the database.

    Parameters:
        request (HttpRequest): The HTTP request containing `reaction_ids` and `confidence_score`.

    Returns:
        JsonResponse:
            - Success: Confirmation of confidence score update.
            - Error: If the input data is invalid.
    """
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            reaction_ids = data.get("reaction_ids", [])
            confidence_score = data.get("confidence_score")

            if not reaction_ids or confidence_score not in ["1", "2", "3", "4"]:
                return JsonResponse({"status": "error", "message": "Invalid data."})

            Reaction.objects.filter(id__in=reaction_ids).update(confidence_score=confidence_score)

            return JsonResponse({"status": "success"})
        except Exception as e:
            return JsonResponse({"status": "error", "message": str(e)})

    return JsonResponse({"status": "error", "message": "Invalid request"})
