"""
This module provides Django views for handling the addition of reactions 
and metabolites to the Virtual Metabolic Human (VMH) database.

Functionalities:
- Retrieving metabolite abbreviations.
- Fetching and updating subsystems.
- Preparing reactions for submission to VMH.
- Validating and generating metabolite and reaction details.
- Sending reaction and metabolite data to VMH via MATLAB.

Dependencies:
- Django
- JSON handling
- Requests for external API calls
- MATLAB integration for reaction/metabolite processing

"""
import json
import os
import requests

from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt

from django.shortcuts import redirect, render, get_object_or_404
from django.views.decorators.http import require_POST
from django.core import serializers 

from django.conf import settings

from reactions.models import (
    User,
    Reaction,
    MetabolitesAddedVMH,
    ReactionsAddedVMH,
    Subsystem,
    SavedMetabolite,
    Workspace
)
from reactions.reaction_info import construct_vmh_formula
from reactions.utils.search_vmh import search_metabolites_vmh, is_name_in_vmh
from reactions.utils.to_smiles import any_to_smiles
from reactions.utils.utils import capitalize_first_letter
from reactions.utils.gen_vmh_abbrs import gen_metabolite_abbr
from reactions.utils.search_vmh import get_from_vmh
from reactions.utils.add_to_vmh_utils import (
    validate_metabolite_existence,
    validate_reaction_existence,
    validate_reaction_fields,
    validate_reaction_objects,
    check_reaction_vmh,
    gather_reaction_details,
    rxn_prepare_json_paths_and_variables,
    met_prepare_json_paths_and_variables,
    add_reaction_matlab,
    add_metabolites_matlab,
    smiles_to_inchikeys,
    smiles_to_charged_formula,
    get_nonfound_metabolites
)
from reactions.utils.MatlabSessionManager import MatlabSessionManager
from reactions.utils.utils import get_external_ids, get_mol_weights, reactions_to_json

def get_metabolite_abbrs(reaction_objs, attr_key, attr_type_key, attr_name_key):
    """
    Retrieve abbreviations for metabolites in a reaction.

    Process:
        - Iterates through reaction metabolites.
        - Retrieves stored abbreviations for saved metabolites.
        - Generates abbreviations for other metabolite types if needed.

    Parameters:
        reaction_objs (list): 
            List of Reaction objects.
        attr_key (str): 
            Attribute name for the metabolite list 
            (e.g., 'substrates' or 'products').
        attr_type_key (str): 
            Attribute name for the metabolite type list 
            (e.g., 'substrates_types' or 'products_types').
        attr_name_key (str): 
            Attribute name for the metabolite name list 
            (e.g., 'substrates_names' or 'products_names').

    Returns:
        list: 
            A list of lists containing abbreviations for each reaction's metabolites.
    """
    abbr_list = []

    for reaction in reaction_objs:
        reaction_abbrs = []
        metabolites = json.loads(getattr(reaction, attr_key))
        metabolite_types = json.loads(getattr(reaction, attr_type_key))
        metabolite_names = json.loads(getattr(reaction, attr_name_key))

        for met, met_type, met_name in zip(metabolites, metabolite_types, metabolite_names):
            if met_type == 'Saved':
                saved_metabolite = SavedMetabolite.objects.get(id=int(met))
                abbr = saved_metabolite.vmh_abbr
                if not abbr:  # Generate if abbreviation doesn't exist
                    abbr = gen_metabolite_abbr(met, met_type, met_name, search_metabolites_vmh)
            else:
                abbr = gen_metabolite_abbr(met, met_type, met_name, search_metabolites_vmh)
            reaction_abbrs.append(abbr)
        abbr_list.append(reaction_abbrs)

    return abbr_list

def get_vmh_subsystems():
    """
    Retrieve subsystem names from the Virtual Metabolic Human (VMH) database.

    Process:
        - Fetches subsystem names from VMH via an API request.
        - Iterates through paginated results to collect all subsystems.

    Returns:
        list:
            A list of subsystem names from VMH.
    """
    base_url = settings.OLD_VMH_BASE_URL
    endpoint = f"{base_url}_api/subsystems/"
    subsystems = []

    # Fetch subsystems from VMH
    while True:
        response = requests.get(endpoint, verify=False,timeout=10)
        data = response.json()['results']
        subsystems.extend([subsystem['name'] for subsystem in data])
        endpoint = response.json().get('next')
        if not endpoint:
            break

    return subsystems


def get_subsystems(request):
    """
    Retrieve subsystem names from VMH and merge them with stored subsystems.

    Process:
        - Fetches subsystems from VMH.
        - Retrieves stored subsystems from the local database.
        - Combines both sets of subsystem names.

    Parameters:
        request (HttpRequest): 
            The HTTP request object.

    Returns:
        JsonResponse:
            - Success: A list of subsystem names.
            - Error: If the request fails.
    """
    try:
        subsystems = get_vmh_subsystems()

        # Fetch subsystems from the database
        stored_subsystems = Subsystem.objects.values_list('name', flat=True)

        # Merge the two lists
        combined_subsystems = set(subsystems).union(set(stored_subsystems))
        combined_subsystems = list(combined_subsystems)

        return JsonResponse({'subsystem_list': combined_subsystems})

    except Exception as e:
        return JsonResponse({'error': True, 'message': str(e)}, status=500)


@csrf_exempt
def update_subsystems(request):
    """
    Update the local database with new subsystems.

    Process:
        - Parses the list of new subsystems from the request.
        - Adds each subsystem to the database if it does not already exist.

    Parameters:
        request (HttpRequest): 
            The HTTP request containing `subsystems` (list of subsystem names).

    Returns:
        JsonResponse:
            - Success: Confirmation of successful update.
            - Error: If an exception occurs or the request method is invalid.
    """
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            subsystems = data.get('subsystems', [])
            # Add new subsystems to the database
            for subsystem in subsystems:
                Subsystem.objects.get_or_create(name=subsystem)

            return JsonResponse({'message': 'Subsystems updated successfully'})
        except Exception as e:
            return JsonResponse({'error': True, 'message': str(e)}, status=500)
    return JsonResponse(
        {'error': True, 'message': 'Invalid request method'}, status=400)


@csrf_exempt
def prepare_add_to_vmh(request):
    """
    Prepare reactions for submission to VMH.

    Process:
        - Validates and fetches reactions by their IDs.
        - Checks if reactions are already in VMH.
        - Retrieves metabolite abbreviations.
        - Identifies metabolites that need new names in VMH.
        - Returns structured data for reaction submission.

    Parameters:
        request (HttpRequest): 
            The HTTP request containing `reactionIds` (list of reaction IDs).

    Returns:
        JsonResponse:
            - Success: Data required for reaction submission.
            - Error: If reactions are missing, already in VMH, or an error occurs.
    """
    if request.method != 'POST':
        return JsonResponse({'status': 'error',
                             'message': 'Invalid request method. Use POST instead.'},
                            status=400)
    try:
        request_data = json.loads(request.body)
        reaction_ids = request_data['reactionIds']
    except json.JSONDecodeError:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid JSON data.'}, status=400)
    except KeyError:
        return JsonResponse({'status': 'error',
                             'message': 'Missing reactionIds in the request data.'},
                            status=400)
    except Exception as e:
        return JsonResponse(
            {'status': 'error', 'message': 'An unexpected error occurred.'}, status=500)

    try:
        reaction_objs = [
            Reaction.objects.get(
                pk=int(reaction_id)) for reaction_id in reaction_ids]

        in_vmh = [
            reaction.vmh_found and not reaction.vmh_found_similar for reaction in reaction_objs]
        if True in in_vmh:
            reaction_objs_in_vmh = [
                reaction for reaction, found in zip(
                    reaction_objs, in_vmh) if found]
            names = [reaction.short_name for reaction in reaction_objs_in_vmh]
            return JsonResponse(
                {'status': 'error',
                 'message': f'The following reactions are already in VMH: {", ".join(names)}'}
            )
        subs_in_vmh = [
            json.loads(reaction.subs_found) if reaction.subs_found else []
            for reaction in reaction_objs
        ]
        prods_in_vmh = [
            json.loads(reaction.prod_found) if reaction.prod_found else []
            for reaction in reaction_objs
        ]
        subs_abbr = get_metabolite_abbrs(
            reaction_objs,
            'substrates',
            'substrates_types',
            'substrates_names'
        )

        prods_abbr = get_metabolite_abbrs(
            reaction_objs,
            'products',
            'products_types',
            'products_names'
        )
        subs_need_new_names = [[] for _ in reaction_ids]
        prods_need_new_names = [[] for _ in reaction_ids]

        for idx, in_vmh_list in enumerate(subs_in_vmh):
            for j, sub_in_vmh in enumerate(in_vmh_list):
                if not sub_in_vmh:
                    need_new_name = is_name_in_vmh(json.loads(
                        reaction_objs[idx].substrates_names)[j])
                    subs_need_new_names[idx].append(need_new_name)
                else:
                    subs_need_new_names[idx].append(False)
        for idx, in_vmh_list in enumerate(prods_in_vmh):
            for j, prod_in_vmh in enumerate(in_vmh_list):
                if prod_in_vmh:
                    prods_need_new_names[idx].append(False)
                else:
                    need_new_name = is_name_in_vmh(json.loads(
                        reaction_objs[idx].products_names)[j])
                    prods_need_new_names[idx].append(need_new_name)

        reaction_abbrs = ['' for _ in reaction_ids]

        return JsonResponse({
            'status': 'success',
            'subs_in_vmh': subs_in_vmh,
            'prods_in_vmh': prods_in_vmh,
            'subs_abbr': subs_abbr,
            'prods_abbr': prods_abbr,
            'subs_need_new_names': subs_need_new_names,
            'prods_need_new_names': prods_need_new_names,
            'reaction_abbrs': reaction_abbrs,
        })
    except Exception as e:
        return JsonResponse({'status': 'error',
                             'message': 'An error occurred while processing reactions.'},
                            status=500)


def bypass_search_func(metabolites, types, *args, **kwargs):
    """
    Dummy function to bypass metabolite search.

    Process:
        - Always returns False for metabolite found status.
        - Returns None for metabolite abbreviations.

    Parameters:
        metabolites (list): List of metabolite identifiers.
        types (list): List of metabolite types.

    Returns:
        tuple:
            - A list of False values (indicating metabolites are not found).
            - A list of None values (indicating no abbreviations).
    """
    return [False], [None]


@csrf_exempt
def create_formula_abbr(request):
    """
    Generate a metabolite abbreviation.

    Process:
        - Parses metabolite details from the request.
        - Generates a new abbreviation using the `gen_metabolite_abbr` function.

    Parameters:
        request (HttpRequest): 
            The HTTP request containing `metabolite`, `mtype`, and `metabolite_name`.

    Returns:
        JsonResponse:
            - Success: The generated abbreviation.
            - Error: If input data is missing or the request method is invalid.
    """
    if request.method == 'POST':
        # Parse JSON data from the request body
        try:
            data = json.loads(request.body)
            metabolite = data.get('metabolite')
            mtype = data.get('mtype')
            metabolite_name = data.get('metabolite_name')
        except json.JSONDecodeError:
            return JsonResponse({'error': 'Invalid JSON'}, status=400)

        # Ensure all required fields are present
        if not all([mtype, metabolite_name]):
            return JsonResponse({'error': 'Missing data'}, status=400)

        # Use the bypass function to force the else clause
        abbr = gen_metabolite_abbr(
            metabolite,
            mtype,
            metabolite_name,
            bypass_search_func)

        # Return the abbreviation as JSON
        return JsonResponse({'abbr': abbr})

    return JsonResponse({'error': 'Invalid request method'}, status=400)

def add_to_vmh(request):
    """
    Submit reactions and metabolites to VMH.

    Process:
        - Validates the user and their permissions.
        - Checks for missing reaction details (e.g., name, abbreviation, confidence score).
        - Ensures reaction names and abbreviations are unique in VMH.
        - Checks whether metabolites already exist in VMH.
        - Sends metabolite and reaction data to VMH via MATLAB.
        - Logs added metabolites and reactions.

    Parameters:
        request (HttpRequest): 
            The HTTP request containing reaction and metabolite details.

    Returns:
        JsonResponse:
            - Success: Confirmation of reactions and metabolites added to VMH.
            - Error: If validation fails or an error occurs during submission.
    """
    req_body = json.loads(request.body)
    user_id= req_body.get('userID')
    user = User.objects.get(pk=user_id)
    if not user:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid user key.'}, status=404)
    if not user.cred_add_to_vmh:
        return JsonResponse({'status': 'error',
                             'message': 'User does not have permission to add to VMH.'},
                            status=403)
    user_name = user.name
    user_full_name = user.name
    reactions = req_body.get('reactions')
    reaction_ids = []
    not_enough_info, no_comments, not_balanced = [], [], []
    met_added_info = {}

    # Validate reaction fields (missing info and duplicates)
    error_response = validate_reaction_fields(reactions)
    if error_response:
        return error_response

    # Validate reaction existence in VMH (names and abbreviations)
    error_response = validate_reaction_existence(reactions)
    if error_response:
        return error_response

    reactions_new_subs_info = [
        json.loads(
            reaction['substrates_info']) for reaction in reactions]
    reactions_new_prods_info = [json.loads(
        reaction['products_info']) for reaction in reactions]
    reactions_subs_found = [
        json.loads(
            Reaction.objects.get(
                pk=reaction['pk']).subs_found) for reaction in reactions]
    reactions_prods_found = [
        json.loads(
            Reaction.objects.get(
                pk=reaction['pk']).prod_found) for reaction in reactions]

    # Validate metabolite existence in VMH
    error_response = validate_metabolite_existence(reactions_new_subs_info, reactions_new_prods_info,
                                                    reactions_subs_found, reactions_prods_found)
    if error_response:
        return error_response

    subs_abbr = []
    prods_abbr = []
    for reaction in reactions:
        obj = Reaction.objects.get(pk=reaction['pk'])
        obj.short_name = reaction['abbreviation']
        obj.description = reaction['description']
        reaction_ids.append(obj.id)
        # Update substrate names and abbreviations
        subs_info = json.loads(reaction['substrates_info'])
        new_subs_names = [
            capitalize_first_letter(
                sub['name']) for sub in subs_info]
        new_subs_abbrs = [sub['abbreviation'] for sub in subs_info]
        subs_abbr.append(new_subs_abbrs)
        obj.substrates_names = json.dumps(new_subs_names)
        # Update product names and abbreviations
        prods_info = json.loads(reaction['products_info'])
        new_prods_names = [
            capitalize_first_letter(
                prod['name']) for prod in prods_info]
        new_prods_abbrs = [prod['abbreviation'] for prod in prods_info]
        prods_abbr.append(new_prods_abbrs)
        obj.products_names = json.dumps(new_prods_names)
        # Update references, external links, and comments
        references, ext_links, comments = json.loads(
            reaction['references']), json.loads(
            reaction['ext_links']), json.loads(
            reaction['comments'])
        new_references, new_ext_links, new_comments = [], [], []
        for ref in references:
            ref['user_name'] = user_name
            if not ('PMID' in ref['info'] or 'DOI' in ref['info']):
                ref['info'] = f"{ref['ref_type']}:{ref['info']}"
            new_references.append(ref)
        for link in ext_links:
            link['user_name'] = user_name
            new_ext_links.append(link)
        for comment in comments:
            comment['user_name'] = user_name
            new_comments.append(comment)
        if f"Created and Added to VMH via Constructor by: {user_full_name}" not in list(
                map(lambda x: x['info'], new_comments)):
            new_comments.append(
                {'info': f"Created and Added to VMH via Constructor by: {user_full_name}",
                  'user_name': user_name}
                )
        not_enough_info.append(
            len(new_references) < 1 and len(new_ext_links) < 1)
        no_comments.append(len(new_comments) < 2)
        obj.references = new_references if new_references else None
        obj.ext_links = new_ext_links if new_ext_links else None
        obj.comments = new_comments if new_comments else None
        obj.confidence_score = reaction.get(
            'confidence_score', 0)  # Add confidence score
        if not (
            json.loads(
                obj.balanced_charge)[0] and json.loads(
                obj.balanced_count)[0]):
            not_balanced.append(True)
        else:
            not_balanced.append(False)
        # Save the updated reaction object
        obj.save()
    reaction_objs = [
        Reaction.objects.get(
            pk=reaction_id) for reaction_id in reaction_ids]

    error_response = validate_reaction_objects(
        reaction_objs, 
        not_enough_info, 
        no_comments, 
        not_balanced
    )    
    if error_response:
        return error_response

    all_vmh = all(
        all(element == 'VMH' for element in json.loads(reaction.substrates_types)) and
        all(element == 'VMH' for element in json.loads(reaction.products_types))
        for reaction in reaction_objs
    )

    reaction_identifiers, reaction_names = [
        reaction['abbreviation'] for reaction in reactions], [
        reaction.description for reaction in reaction_objs]
    matlab_session = None  if all_vmh  else  MatlabSessionManager() 
    if not all_vmh:
        unique_abbrs, unique_mols, unique_types, unique_names = get_nonfound_metabolites(
            reaction_objs, subs_abbr, prods_abbr, search_func=search_metabolites_vmh)
        if unique_abbrs and unique_mols and unique_types:
            abbrs = unique_abbrs
            smiles, errors = any_to_smiles(
                unique_mols, unique_types, None, side=None)
            smiles = ['' if v is not None else smiles[i]
                      for i, v in enumerate(errors)]
            inchikeys = smiles_to_inchikeys(smiles)
            names = unique_names
            formulas, charges = smiles_to_charged_formula(smiles)
            external_ids_list = get_external_ids(unique_mols, unique_types)
            mol_weights = get_mol_weights(unique_mols, unique_types)
            json_paths = met_prepare_json_paths_and_variables(
                abbrs, names, formulas, charges, inchikeys, smiles, external_ids_list, mol_weights)
            matlab_result = add_metabolites_matlab(json_paths, matlab_session)
            for path in json_paths:
                os.remove(path)
            if matlab_result['status'] == 'success':
                met_ids = matlab_result['met_ids']
                met_added_info = {
                    abbr: [
                        met_id,
                        formula,
                        inchikey] for abbr,
                    met_id,
                    formula,
                    inchikey in zip(
                        abbrs,
                        met_ids,
                        formulas,
                        inchikeys)}
                for abbr, info in met_added_info.items():
                    MetabolitesAddedVMH.objects.create(
                        user=user,
                        user_name=user_name,
                        metabolite_id=info[0],
                        metabolite_formula=info[1],
                        metabolite_abbr=abbr,
                    )
            else:
                return JsonResponse(
                    {'status': 'error', 'message': matlab_result['message']})
    reaction_formulas = [
        construct_vmh_formula(
            reaction_objs[idx],
            subs_abbr[idx],
            prods_abbr[idx]) for idx in range(
            len(reaction_objs))]
    (
        reaction_directions,
        reaction_subsystems,
        reaction_references,
        reaction_external_links,
        reaction_gene_info,
        reaction_comments,
        reaction_confidence_scores
    ) = gather_reaction_details(reaction_objs)

    json_paths = rxn_prepare_json_paths_and_variables(
        reaction_identifiers,
        reaction_names,
        reaction_formulas,
        reaction_directions,
        reaction_subsystems,
        reaction_references,
        reaction_external_links,
        reaction_gene_info,
        reaction_comments,
        reaction_confidence_scores)
    matlab_session = MatlabSessionManager() if not matlab_session else matlab_session
    matlab_result = add_reaction_matlab(json_paths, matlab_session)
    # Cleanup temporary JSON files
    for path in json_paths:
        os.remove(path)

    if matlab_result['status'] == 'success':
        rxn_added_info = {
            abbr: [
                rxn_id,
                reaction_formula] for abbr,
            rxn_id,
            reaction_formula in zip(
                reaction_identifiers,
                matlab_result['rxn_ids'],
                reaction_formulas)}
        for idx, (abbr, info) in enumerate(rxn_added_info.items()):
            ReactionsAddedVMH.objects.create(
                user=user,
                user_name=user_name,
                reaction_id=info[0],         # VMH ID
                reaction_formula=info[1],
                reaction_abbr=abbr,
            )

            workspace = Workspace.objects.get(user=user)
            reaction_obj = reaction_objs[idx]  
            reaction_obj.vmh_found = True
            workspace.reactions.remove(reaction_obj)
        return JsonResponse({'status': 'success',
                             'rxn_added_info': rxn_added_info,
                             'met_added_info': met_added_info})

    return JsonResponse(
        {'status': 'error', 'message': matlab_result['message']})
    
    
@require_POST
def send_to_workspace(request):
    """
    Handles the request to add selected reactions to the user's Workspace.

    Process:
    - Parses JSON data from the POST request body.
    - Retrieves the user ID and a list of selected reaction primary keys (IDs).
    - Fetches the corresponding User object.
    - Retrieves or creates a Workspace object associated with the user.
    - Adds the selected Reaction objects to the Workspace's many-to-many field.
    - Saves the updated Workspace.

    Parameters:
    - request (HttpRequest): The HTTP POST request containing JSON body with:
        - 'userID' (int): ID of the user.
        - 'reactionIds' (list of int): List of reaction primary keys to add.

    Returns:
    - JsonResponse: A JSON response with a success status:
        { "status": "success" }
    """
    data = json.loads(request.body)
    user_id = data.get('userID')
    reaction_ids = data.get('reactionIds', [])
    user = User.objects.get(pk=user_id)
    # Save to a Workspace model (create if not exists)
    workspace, _ = Workspace.objects.get_or_create(user=user)
    workspace.reactions.add(*Reaction.objects.filter(pk__in=reaction_ids))
    workspace.save()
    return JsonResponse({'status': 'success'})         


def vmh_workspace(request):
    """
    View function to render the VMH Workspace page with user's available and added reactions.

    Process:
    - Retrieves the logged-in user's ID from the session.
    - If the user is authenticated:
        - Retrieves the corresponding User object.
        - Gets or creates a Workspace object for the user.
        - Fetches all reactions in the user's workspace.
        - Fetches all reactions the user has already added to VMH, ordered by most recent.
    - If the user is not authenticated:
        - Uses empty QuerySets for both workspace and added reactions.
    - Serializes the data for both available and added reactions to pass to the frontend.

    Parameters:
    - request (HttpRequest): The HTTP request object, which may contain the session data.

    Returns:
    - HttpResponse: Renders 'reactions/VMH_workspace.html' with the following context:
        - 'user_name' (str): Name of the logged-in user (or 'Guest' if unauthenticated).
        - 'userID' (int or str): ID of the user (or empty string if unauthenticated).
        - 'reactions_active' (JSON): List of available (non-added) reactions from the workspace.
        - 'reactions_added' (JSON): List of reactions already added to VMH.
    """
    user_pk = request.session.get('userID')

    if user_pk:
        user_obj = get_object_or_404(User, pk=user_pk)
        user_name = user_obj.name  # or .username depending on your model
        user_id = user_obj.pk
        # Get or create the workspace for this user
        workspace, _ = Workspace.objects.get_or_create(user=user_obj)
        reactions_qs = workspace.reactions.all()

        added = ReactionsAddedVMH.objects.filter(user=user_id).order_by('-created_at')
    else:
        user_name = 'Guest'
        user_id = ''
        reactions_qs = Reaction.objects.none()
        added = ReactionsAddedVMH.objects.none()

    reactions_active = reactions_to_json(reactions_qs)
    added_json = serializers.serialize('json', added)

    ctx = {
        'user_name': user_name,
        'userID': user_id,
        'reactions_active': reactions_active,
        'reactions_added': added_json,
    }
    return render(request, 'reactions/VMH_workspace.html', ctx)