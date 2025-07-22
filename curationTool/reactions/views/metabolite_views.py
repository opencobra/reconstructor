"""
This module provides API endpoints for handling metabolites.

It includes:
- `verify_metabolite`: Checks if a metabolite exists in VMH, SavedMetabolites, or other sources.
- `get_saved_metabolites`: Retrieves all metabolites saved by a user.
- `fetch_rhea_rxn`: Fetches reaction data from the Rhea database.
- `update_metabolite`: Updates a saved metabolite's details.
- `delete_metabolite`: Deletes a metabolite, transferring ownership if necessary.
- `share_metabolites`: Shares a metabolite with another user.
- `generate_abbreviation`: Generates an abbreviation for a given metabolite name.
"""

import json
import os
from collections import defaultdict
import requests

from django.views.decorators.csrf import csrf_exempt
from django.forms.models import model_to_dict
from django.http import JsonResponse
from django.db.models import Q
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula #  pylint: disable=no-name-in-module

from reactions.utils.search_vmh import search_vmh
from reactions.utils.to_mol import any_to_mol
from reactions.utils.utils import capitalize_first_letter, parse_mol_formula
from reactions.utils.get_from_rhea import get_from_rhea
from reactions.models import SavedMetabolite, User

from django.conf import settings

def verify_metabolite(request):
    """
    Verify whether a metabolite exists in VMH or the user's saved database.

    Process:
        - Checks VMH API for the metabolite.
        - Checks if the metabolite is saved in the database for a logged-in user.
        - Extracts molecular information such as formula, charge, and structure.

    Parameters:
        request (HttpRequest): The HTTP request containing metabolite details.

    Returns:
        JsonResponse:
            - Success: Details of the metabolite (name, formula, charge, etc.).
            - Error: If the metabolite is not found or an invalid request is made.
    """
    main_input = request.POST.get('metabolite')
    input_type = request.POST.get('type')

    if not main_input.strip():
        return JsonResponse({'error': True, 'message': 'No input provided'})

    if input_type == 'VMH':
        return _verify_vmh_metabolite(main_input, input_type)
    if input_type == 'Saved':
        return _verify_saved_metabolite(request, main_input, input_type)
    return _verify_other_metabolite(request, main_input, input_type)

def _verify_vmh_metabolite(main_input, input_type):
    """Helper function to verify a metabolite in VMH."""
    base_url = settings.OLD_VMH_BASE_URL
    endpoint = f"{base_url}_api/metabolites/?abbreviation={main_input}"

    try:
        response = requests.get(endpoint, verify=False, timeout=10)
    except Exception:
        response = None

    results = []
    if response and response.status_code == 200:
        data = response.json()
        results = data.get('results', [])

    if results:
        inchi_string = results[0].get('inchiString', '')
        smile = results[0].get('smile', '')

        if not smile and not inchi_string:
            return JsonResponse({
                'error': True,
                'message': f"Metabolite {main_input} does not have SMILES/Inchi in VMH"
            }, status=404)

        mol_obj = Chem.MolFromSmiles(smile) if smile else Chem.MolFromInchi(inchi_string)
        if not mol_obj:
            return JsonResponse({
                'error': True,
                'message': f"Failed to parse structure for {main_input}"
            }, status=500)

        formula = CalcMolFormula(mol_obj)
        atom_counts = parse_mol_formula(formula)
        charge = sum(atom.GetFormalCharge() for atom in mol_obj.GetAtoms())

        return JsonResponse({
            'found': True,
            'abbr': results[0]['abbreviation'],
            'name': results[0]['fullName'],
            'miriam': results[0]['miriam'],
            'input_type': input_type,
            'atom_counts': atom_counts,
            'charge': charge,
            'noStructure': False
        })

    # Fallback: Try to load local .mol file
    mol_path = os.path.join(settings.MOL_FILE_PATH , f"{main_input}.mol")
    if os.path.exists(mol_path):
        mol_obj = Chem.MolFromMolFile(mol_path)
        if mol_obj:
            formula = CalcMolFormula(mol_obj)
            atom_counts = parse_mol_formula(formula)
            charge = sum(atom.GetFormalCharge() for atom in mol_obj.GetAtoms())

            return JsonResponse({
                'found': True,
                'abbr': main_input,
                'name': main_input,
                'miriam': '',
                'input_type': input_type,
                'atom_counts': atom_counts,
                'charge': charge,
                'noStructure': False
            })
        else:
            return JsonResponse({
                'error': True,
                'message': f"Failed to parse local MOL file for `{main_input}`"
            }, status=500)

    return JsonResponse({
        'error': True,
        'message': f"Metabolite `{main_input}` not found in VMH or storage"
    }, status=404)

def _verify_saved_metabolite(request, main_input, input_type):
    """Helper function to verify a metabolite in the saved database."""
    user_id = request.POST.get('userID')

    try:
        user = User.objects.get(pk=user_id) #  pylint: disable=no-member
    except (User.DoesNotExist, ValueError): #  pylint: disable=no-member
        return JsonResponse(
            {'error': True, 'message': "Cannot use saved metabolite without signing in"}
        )

    try:
        saved_met = SavedMetabolite.objects.get( #  pylint: disable=no-member
            Q(owner=user) | Q(shared_with=user),
            pk=main_input
        )
    except SavedMetabolite.DoesNotExist: #  pylint: disable=no-member
        return JsonResponse(
            {
                'error': True,
                'message': f"Saved metabolite with id `{main_input}` not found"
            },
            status=404
        )

    atom_counts = parse_mol_formula(saved_met.mol_formula) if saved_met.mol_formula else {}

    try:
        mol_obj = Chem.MolFromMolBlock(saved_met.mol_file)
        charge = sum(atom.GetFormalCharge() for atom in mol_obj.GetAtoms()) if mol_obj else None
    except Exception:
        charge = None

    return JsonResponse({
        'found': False,
        'abbr': saved_met.vmh_abbr,
        'name': saved_met.name,
        'input_type': input_type,
        'atom_counts': atom_counts,
        'charge': charge
    })


def _verify_other_metabolite(request, main_input, input_type):
    """Helper function to verify metabolites from other sources."""
    mols, errors, names = any_to_mol([main_input], [input_type], request, side=None)
    mol, error, name = mols[0], errors[0], names[0]

    if error:
        return JsonResponse({'error': True, 'message': error})

    found, miriam, abbr, name_vmh = search_vmh(mol, return_abbr=True, return_name=True)
    name = capitalize_first_letter(name_vmh if found else name)
    inchi_key = Chem.MolToInchiKey(mol)
    user_id = request.POST.get('userID')

    user = User.objects.get(pk=user_id) if user_id else None #  pylint: disable=no-member
    saved_exists, name_in_db = False, ""

    if user:
        saved_met_obj = SavedMetabolite.objects.filter( #  pylint: disable=no-member
            Q(owner=user) | Q(shared_with=user),
            inchi_key=inchi_key
        ).distinct()

        if saved_met_obj.exists():
            saved_exists = True
            name_in_db = saved_met_obj.first().name

    formula = CalcMolFormula(mol)
    atom_counts = parse_mol_formula(formula)
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())

    return JsonResponse({
        'found': found,
        'abbr': abbr,
        'name': name,
        'miriam': miriam,
        'input_type': input_type,
        'saved_exists': saved_exists,
        'name_in_db': name_in_db,
        'atom_counts': atom_counts,
        'charge': charge,
        'noStructure': False
    })



def get_saved_metabolites(request):
    """
    Retrieve a list of metabolites saved by a specific user.

    Process:
        - Fetches metabolites owned or shared with the user.
        - Retrieves linked reactions for each metabolite.
        - Formats response with metabolite details.

    Parameters:
        request (HttpRequest): The HTTP request containing the `user_id`.

    Returns:
        JsonResponse:
            - Success: A list of metabolites with relevant details.
            - Error: If the user ID is missing or invalid.
    """
    user_id = request.GET.get('user_id')
    if not user_id:
        return JsonResponse({'metabolites': []})

    try:
        user = User.objects.get(pk=user_id)
    except User.DoesNotExist:
        return JsonResponse({'error': True, 'message': 'User not found'}, status=404)

    saved_metabolites = SavedMetabolite.objects.filter(
        Q(owner=user) | Q(shared_with=user)
    ).distinct()

    metabolite_reactions = _get_metabolite_reactions(user.saved_reactions.all())

    metabolites_list = [
        {
            'id': met.id,
            'name': met.name,
            'vmh_abbr': met.vmh_abbr or '',
            'reactions': metabolite_reactions.get(str(met.id), []),
            'inchi_key': met.inchi_key,
            'inchi': met.inchi,
            'smiles': met.smiles,
            'mol_w': met.mol_w,
            'mol_formula': met.mol_formula,
            'mol_file': met.mol_file,
            'external_links': {
                'keggId': met.keggId,
                'pubChemId': met.pubChemId,
                'cheBlId': met.cheBlId,
                'hmdb': met.hmdb,
                'metanetx': met.metanetx,
            }
        }
        for met in saved_metabolites
    ]

    return JsonResponse({'metabolites': metabolites_list})

def _get_metabolite_reactions(user_reactions):
    """Helper function to process user reactions and link them to metabolites."""
    metabolite_reactions = defaultdict(list)
    for reaction in user_reactions:
        try:
            subs = json.loads(reaction.substrates)
            prods = json.loads(reaction.products)
            subs_types = json.loads(reaction.substrates_types)
            prods_types = json.loads(reaction.products_types)
        except json.JSONDecodeError:
            continue

        reaction_info = {
            'id': reaction.id,
            'name': (
                reaction.short_name
                or f"Reaction {reaction.id}"
            )
        }

        for met_id, met_type in zip(subs + prods, subs_types + prods_types):
            if met_type == 'Saved':
                metabolite_reactions[str(met_id)].append(reaction_info)

    return metabolite_reactions

def fetch_rhea_rxn(request):
    """
    Fetch reaction data from the Rhea database based on a reaction abbreviation.

    Process:
        - Extracts the reaction abbreviation from the request.
        - Calls the Rhea API to retrieve reaction details.

    Parameters:
        request (HttpRequest): The HTTP request containing reaction details.

    Returns:
        JsonResponse:
            - Success: The reaction details from Rhea.
            - Error: If the request is invalid or the reaction abbreviation is missing.
    """
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            reaction_abbreviation = data.get('reactionAbbreviation')
            if reaction_abbreviation:
                result = get_from_rhea(reaction_abbreviation)
                return JsonResponse(result)

            return JsonResponse(
                {'error': 'Missing reaction abbreviation'}, status=400)
        except json.JSONDecodeError:
            return JsonResponse({'error': 'Invalid JSON'}, status=400)
    else:
        return JsonResponse({'error': 'Invalid request'}, status=400)

@csrf_exempt
def update_metabolite(request, metabolite_id):
    """
    Update details of a saved metabolite.

    Process:
        - Retrieves the metabolite from the database.
        - Updates its fields based on the provided JSON data.
        - Ensures the `name` field is not empty.

    Parameters:
        request (HttpRequest): The HTTP request containing updated metabolite details.
        metabolite_id (int): The ID of the metabolite to update.

    Returns:
        JsonResponse:
            - Success: The updated metabolite details.
            - Error: If the metabolite is not found or the request is invalid.
    """
    if request.method != 'PUT':
        return JsonResponse({'status': 'error', 'message': 'Invalid request method'}, status=405)

    try:
        metabolite = SavedMetabolite.objects.get(id=metabolite_id) 
        data = json.loads(request.body)

        for field, value in data.items():
            if field == 'name' and not value.strip():
                return JsonResponse(
                    {'status': 'error', 'message': f"{field} field cannot be empty"},
                    status=400
                )
            setattr(metabolite, field, value)

        metabolite.save()
        updated_metabolite = model_to_dict(metabolite, exclude=['shared_with'])

        return JsonResponse({'status': 'success', 'metabolite': updated_metabolite})

    except SavedMetabolite.DoesNotExist:
        return JsonResponse({'status': 'error', 'message': "Metabolite not found"}, status=404)

    except json.JSONDecodeError:
        return JsonResponse({'status': 'error', 'message': 'Invalid JSON format'}, status=400)

    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)}, status=500)

def _get_affected_reactions(user, metabolite_id):
    """Return reaction objects from user's saved reactions that use the given metabolite."""
    return [
        reaction for reaction in user.saved_reactions.all()
        if str(metabolite_id) in map(str, json.loads(reaction.substrates))
        or str(metabolite_id) in map(str, json.loads(reaction.products))
    ]

@csrf_exempt
def delete_metabolite(request, metabolite_id):
    """
    Delete a metabolite, transferring ownership if it is shared.

    Process:
        - Checks if the user is the owner or a shared user.
        - Transfers ownership if shared, otherwise deletes the metabolite.
        - Removes associated reactions if necessary.

    Parameters:
        request (HttpRequest): The HTTP request containing user details.
        metabolite_id (int): The ID of the metabolite to delete.

    Returns:
        JsonResponse:
            - Success: Confirmation of deletion or ownership transfer.
            - Error: If the metabolite or user is not found.
    """
    if request.method != 'DELETE':
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method'},
            status=405
        )
    try:
        metabolite = SavedMetabolite.objects.get(id=metabolite_id) #  pylint: disable=no-member
        user_id = json.loads(request.body).get('user_id')
        user = User.objects.get(id=user_id) if user_id else None #  pylint: disable=no-member
        # Check if user is not the owner
        if metabolite.owner != user:
            if user in metabolite.shared_with.all():
                metabolite.shared_with.remove(user)
                metabolite.save()
                return JsonResponse(
                    {
                        'status': 'success',
                        'message': "Metabolite removed from your list."
                    }
                )
            # Else, user is not the owner or shared user
            return JsonResponse(
                {
                    'status': 'error',
                    'message': "Permission denied."
                },
                status=403
            )
        # User is the owner
        # Check if metabolite is shared and transfer ownership if possible
        shared_users = list(metabolite.shared_with.all())
        if shared_users:
            new_owner = shared_users[0]
            metabolite.owner = new_owner
            metabolite.shared_with.remove(new_owner)
            if user in metabolite.shared_with.all():
                metabolite.shared_with.remove(user)
            metabolite.save()

            # Delete owner's reactions using this metabolite
            affected_reactions = _get_affected_reactions(user, metabolite_id)
            # Remove reactions from user's saved_reactions
            user.saved_reactions.remove(*affected_reactions)
            return JsonResponse(
                {
                    'status': 'success',
                    'message': "Ownership transferred and reactions removed."
                }
            )
        # Find all reactions containing this metabolite
        affected_reactions = _get_affected_reactions(user, metabolite_id)
        confirm = request.GET.get('confirm') == 'true'
        if not confirm and affected_reactions:
            reactions_list = [
                {'id': r.id, 'short_name': r.short_name} for r in affected_reactions
            ]
            return JsonResponse(
                {
                    'status': 'needs_confirmation',
                    'message': (
                        "IMPORTANT: DELETING THIS METABOLITE WILL DELETE ALL REACTIONS THAT "
                        "CONTAIN IT. ARE YOU SURE YOU WANT TO PROCEED?"
                    ),
                    'reactions': reactions_list
                },
                status=200
            )
        # Delete affected reactions and metabolite
        for reaction in affected_reactions:
            user.saved_reactions.remove(reaction)  # Use the ManyToManyField manager
        metabolite.delete()
        return JsonResponse({'status': 'success'})
    except SavedMetabolite.DoesNotExist: #  pylint: disable=no-member
        return JsonResponse(
            {
                'status': 'error',
                'message': "Metabolite not found"
            },
            status=404
        )
    except User.DoesNotExist: #  pylint: disable=no-member
        return JsonResponse(
            {
                'status': 'error',
                'message': "User not found."
            },
            status=404
        )
    except Exception as e:
        return JsonResponse(
            {
                'status': 'error',
                'message': str(e)
            },
            status=500
        )

def share_metabolites(request):
    """
    Share a metabolite with another user.

    Process:
        - Validates that the user is the owner of the metabolite.
        - Checks that the target user exists and is not the same as the sharing user.
        - Adds the target user to the `shared_with` field of the metabolite.

    Parameters:
        request (HttpRequest): The HTTP request containing metabolite IDs and user details.

    Returns:
        JsonResponse:
            - Success: Confirmation of shared metabolites.
            - Error: If the user is not the owner or the target user is invalid.
    """
    try:
        data = json.loads(request.body)
        metabolite_ids = data.get('metabolite_ids', [])
        target_username = data.get('target_username')
        sharing_user_id = data.get('sharing_user_id')
        sharing_user = User.objects.get(pk=sharing_user_id) #  pylint: disable=no-member
        if not sharing_user:
            return JsonResponse({'status': 'error', 'message': 'Not logged in'}, status=401)
        if not metabolite_ids or not target_username:
            return JsonResponse(
                {
                    'status': 'error',
                    'message': (
                        "Metabolite IDs and target username are required."
                    )
                },
                status=400
            )
        try:
            target_user = User.objects.get(name=target_username) #  pylint: disable=no-member

        except User.DoesNotExist: #  pylint: disable=no-member
            return JsonResponse(
                {
                    'status': 'error',
                    'message': "Target user not found."
                },
                status=404
            )

        if target_user == sharing_user:
            return JsonResponse(
                {
                    'status': 'error',
                    'message': "Cannot share metabolites with yourself."
                },
                status=400
            )
        #Check if the user is shared_with not owner for any of the metabolites
        metabolites_shared_with_user = SavedMetabolite.objects.filter( #  pylint: disable=no-member
            id__in=metabolite_ids,
            shared_with=sharing_user
        )
        if metabolites_shared_with_user.exists():
            names = [met.name for met in metabolites_shared_with_user]
            return JsonResponse(
                {
                    'status': 'error',
                    'message': (
                        "Must be the owner to share metabolites.\n"
                        "You are not the owner of the following metabolites: " 
                        + ", ".join(names)
                    )
                },
                status=400
            )
        # Only share metabolites where the current user is the owner.
        metabolites = SavedMetabolite.objects.filter(id__in=metabolite_ids, owner=sharing_user) #  pylint: disable=no-member
        if not metabolites.exists():
            return JsonResponse(
                {
                    'status': 'error',
                    'message': "No valid metabolites found to share."
                },
                status=400
            )
        shared_count = 0
        for metabolite in metabolites:
            if target_user not in metabolite.shared_with.all():
                metabolite.shared_with.add(target_user)
                shared_count += 1
        already_shared = len(metabolite_ids) - shared_count
        if already_shared > 0:
            msg = (
                f"Shared {shared_count} metabolites with {target_username}. "
                f"{already_shared} metabolites were already shared."
            )
        else:
            msg = f'Shared {shared_count} metabolites with {target_username}.'

        return JsonResponse({'status': 'success', 'message': msg})

    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)}, status=500)

def generate_abbreviation(request):
    """
    Generate an abbreviation for a given metabolite name.

    Process:
        - Extracts the metabolite name from the request.
        - Generates an abbreviation by taking the first three uppercase letters.

    Parameters:
        request (HttpRequest): The HTTP request containing the metabolite name.

    Returns:
        JsonResponse:
            - Success: The generated abbreviation.
            - Error: If no name is provided in the request.
    """
    if request.method != 'POST':
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method'},
            status=405
        )

    # Implement your abbreviation generation logic here
    name = json.loads(request.body).get('name')
    generated_abbr = name[:3].upper()  # Example simple generation
    return JsonResponse({'abbr': generated_abbr})
