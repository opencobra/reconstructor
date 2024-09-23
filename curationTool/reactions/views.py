import json
import re
import pandas as pd 
import sys
import logging
import os
import time
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from django.utils.decorators import method_decorator
from urllib.parse import unquote
from django.shortcuts import render, redirect
from django.urls import reverse
from django.http import JsonResponse, Http404
from django.core import serializers
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth.hashers import make_password
from django.core.exceptions import ObjectDoesNotExist
from .forms import ReactionForm
from .models import User,CreatedReaction,Reaction,Flag
from django.shortcuts import get_object_or_404
from django.views.decorators.http import require_POST
from .utils.get_from_rhea import get_from_rhea
from .models import Reaction, User, MetabolitesAddedVMH, ReactionsAddedVMH, Subsystem,CreatedReaction
from reactions.reaction_info import get_reaction_info, construct_vmh_formula
from reactions.utils.process_strings import construct_reaction_string, construct_reaction_rxnfile, get_mol_names
from reactions.utils.get_mol_info import get_mol_info
from reactions.utils.search_vmh import search_metabolites_vmh, check_reaction_vmh, get_from_vmh, search_vmh, is_name_in_vmh
from reactions.utils.to_smiles import any_to_smiles
from reactions.utils.to_mol import any_to_mol
from reactions.utils.utils import gen_replace_dict, get_fields, seperate_metab_names, capitalize_first_letter, parse_xml,check_edited_keys,clean_dict_keys ,fetch_and_map_gene_expression
from reactions.utils.RDT import RDT
from reactions.utils.gen_vmh_abbrs import gen_reaction_abbr, gen_metabolite_abbr
from reactions.utils.add_to_vmh_utils import check_met_names_abbrs_vmh,check_names_abbrs_vmh, gather_reaction_details, rxn_prepare_json_paths_and_variables, met_prepare_json_paths_and_variables, add_reaction_matlab, add_metabolites_matlab, smiles_to_inchikeys, smiles_to_charged_formula, get_nonfound_metabolites, check_reactions_vmh
from reactions.utils.MatlabSessionManager import MatlabSessionManager
from django.conf import settings
from django.contrib.auth import authenticate
import requests
from reactions.utils.utils import get_subcellular_locations, map_locations_to_wbm
from time import sleep
from urllib3.exceptions import InsecureRequestWarning
import os
from reactions.utils.GPT_functions import get_gene_name, get_vmh_met_from_inchi, metanetx_to_inchi, get_ncbi_gene_id, get_vmh_synonyms, get_gene_reactions, extract_compounds, vmh_to_normal, get_cid_vmh_api, get_cid_api, get_metanetx_id, get_metanetx_id_by_name, parse_metabolic_reactions, parse_metabolic_reactions_gpt, help_format_answer_with_gpt, askGPT4, pubchem_similarity_search, jaccard_similarity, flatten_extend, evaluate_predictions, check_match, print_output
from django.views.decorators.http import require_GET
def about_view(request):
    return render(request, 'reactions/about.html')


def get_user_flags(request, user_id):
    try:
        user = User.objects.get(pk=user_id)
        flags = user.flags.all()
        flags_data = [{'id': flag.id, 'name_flag': flag.name_flag, 'color': flag.color} for flag in flags]
        return JsonResponse({'status': 'success', 'flags': flags_data})
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)}, status=500)

@csrf_exempt
def add_flag(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('user_id')
            flag_name = data.get('name_flag')
            flag_color = data.get('color')
            user = User.objects.get(pk=user_id)
        except User.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Invalid user'})
        
        if flag_name and flag_color:
            flag_name = flag_name.strip()
            flag, created = Flag.objects.get_or_create(name_flag=flag_name, color=flag_color, user=user)
            return JsonResponse({
                'status': 'success',
                'message': 'Flag added successfully',
                'flag': {'id': flag.id, 'name_flag': flag.name_flag, 'color': flag.color}
            })
        else:
            return JsonResponse({'status': 'error', 'message': 'Flag name and color are required'})
    return JsonResponse({'status': 'error', 'message': 'Invalid request method'})

def get_vmh_subsystems():
    BASE_URL = 'https://www.vmh.life/'
    endpoint = f"{BASE_URL}_api/subsystems/"
    subsystems = []

    # Fetch subsystems from VMH
    while True:
        response = requests.get(endpoint, verify=False)
        data = response.json()['results']
        subsystems.extend([subsystem['name'] for subsystem in data])
        endpoint = response.json().get('next')
        if not endpoint:
            break

    return subsystems

def get_subsystems(request):
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
    return JsonResponse({'error': True, 'message': 'Invalid request method'}, status=400)


def get_user(request):
    username = request.POST.get('username', '')
    password = request.POST.get('password', '')
    try:
        user = User.objects.get(name=username)
        if user.check_password(password):
            request.session['userID'] = user.pk  # Store the user ID in the session
            return JsonResponse({'status': 'success', 'userName': user.name, 'userID': user.pk})
        else:
            return JsonResponse({'status': 'error', 'message': 'Invalid password'})
    except User.DoesNotExist:
        return JsonResponse({'status': 'error', 'message': 'User does not exist'})
    

def register_user(request):
    username = request.POST.get('username', '')
    password = request.POST.get('password', '')
    email = request.POST.get('email', '')
    orchid_id = request.POST.get('orchid_id', '')
    full_name = request.POST.get('full_name', '')   # New field
    users = User.objects.all()
    usernames = [user.name for user in users]
    if username in usernames:
        return JsonResponse({'status': 'error', 'message': 'Username already exists'})
    if not username or not password or not email or not full_name:
        return JsonResponse({'status': 'error', 'message': 'Ysername, password, email, and full name are required'})

    try:
        if User.objects.filter(name=username).exists():
            return JsonResponse({'status': 'error', 'message': 'Username already exists'})
        
        user = User(
            name=username,
            email=email,
            password=password,  # Hash the password before saving
            orchid_id=orchid_id,
            full_name=full_name
        )
        user.save()

        request.session['userID'] = user.pk  # Store the user ID in the session
        return JsonResponse({'status': 'success', 'userName': user.name, 'userID': user.pk})
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)})

def set_session_user(request):
    req_body = json.loads(request.body)
    userID = int(req_body.get('userID'))
    request.session['userID'] = userID
    if request.session.get('userID') == userID and validate_user_ID(userID):
        return JsonResponse({'status': 'success', 'message': 'User set in session'})
    else:
        return JsonResponse({'status': 'error', 'message': 'User not set in session'})

def chemdoodle_sketcher(request):
    return render(request, 'reactions/chemdoodle_sketcher.html')


def get_doi_info(request, doi):
    doi = doi.replace('DOI:', '')
    base_url = f'https://api.crossref.org/works/{doi}'
    try:
        response = requests.get(base_url)
        if response.status_code == 200:
            data = response.json()['message']
            
            # Extract needed information
            authors = [f"{author['given']} {author['family']}" for author in data.get('author', [])]
            title = data.get('title', [None])[0]
            abstract = data.get('abstract', None)
            return JsonResponse({
                'status': 'success',
                'author': authors,
                'title': title,
                'abstract': abstract
            })
        else:
            return JsonResponse({'status': 'error', 'message': 'CrossRef API returned an error'}, status=response.status_code)
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': 'Failed to fetch DOI data'}, status=500)


def get_pubmed_info(request, pmid):
    # Base URL for PubMed API
    base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    params = {
        'db': 'pubmed',
        'id': pmid,
        'retmode': 'xml'
    }
    response = requests.get(base_url, params=params)
    
    if response.status_code == 200:
        # Parsing XML response is necessary here to extract needed information
        # This is a placeholder to indicate where parsing should occur
        # You might use libraries like xml.etree.ElementTree or lxml to parse the XML
        
        pubmed_info = parse_xml(response.content) # Implement this function
        if pubmed_info['title'] == None:
            return JsonResponse({'status': 'error', 'message': 'PubMed API returned an empty response'}, status=500)
        return JsonResponse({
            'status': 'success',
            'author': pubmed_info['authors'], 
            'title': pubmed_info['title'],
            'abstract': pubmed_info['abstract']
        })
    else:
        return JsonResponse({'status': 'error', 'message': f"PubMed API returned error {response.status_code}"}, status=500)

import requests
from django.http import JsonResponse

def get_gene_info(request):
    """
    Retrieves gene information based on user input.

    Args:
        request (HttpRequest): The HTTP request object.

    Returns:
        JsonResponse: The JSON response containing the gene information.

    Raises:
        None
    """
    gene_input = request.POST.get('gene')
    type_input = request.POST.get('type')
    
    if gene_input.strip() == '':
        return JsonResponse({'error': True, 'message': 'No input provided'})

    if type_input == 'Entrez ID':
        # Check Entrez ID
        entrez_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        entrez_params = {
            "db": "gene",
            "id": gene_input,
            "retmode": "json"
        }
        
        entrez_response = requests.get(entrez_base_url, params=entrez_params)
        
        if entrez_response.status_code == 200:
            entrez_data = entrez_response.json()
            if "result" in entrez_data and gene_input in entrez_data["result"]:
                gene_data = entrez_data["result"][gene_input]
                if "name" in gene_data:
                    return JsonResponse({'error': False, 'symbol': gene_data["name"]})
        
        # If not found in Entrez, check VMH
        vmh_base_url = 'https://www.vmh.life/'
        vmh_endpoint = f"{vmh_base_url}_api/genes/?gene_number={gene_input}"
        vmh_response = requests.get(vmh_endpoint, verify=False)
        
        if vmh_response.status_code != 200:
            return JsonResponse({'error': True, 'message': f'VMH API returned error {vmh_response.status_code} for gene number `{gene_input}`'}, status=500)
        
        vmh_data = vmh_response.json()
        if vmh_data['count'] == 0:
            return JsonResponse({'error': True, 'message': f'Gene number `{gene_input}` not found in VMH'}, status=404)
        
        gene = vmh_data['results'][0]
        symbol = gene.get('symbol', '')
        return JsonResponse({'error': False, 'symbol': symbol})
    
    elif type_input == 'HGNC Symbol':
        # Check HGNC Symbol
        hgnc_base_url = 'https://rest.genenames.org/search/symbol/'
        hgnc_endpoint = f"{hgnc_base_url}{gene_input}"
        hgnc_headers = {'Accept': 'application/json'}
        hgnc_response = requests.get(hgnc_endpoint, headers=hgnc_headers)
        
        if hgnc_response.status_code != 200:
            return JsonResponse({'error': True, 'message': f'HGNC API returned error {hgnc_response.status_code} for symbol `{gene_input}`'}, status=500)
        
        hgnc_data = hgnc_response.json()
        num_found = hgnc_data['response']['numFound']
        
        if num_found == 0:
            return JsonResponse({'error': True, 'message': f'Gene symbol `{gene_input}` not found in HGNC'}, status=404)
        elif num_found > 1:
            genes = hgnc_data['response']['docs'][:10]
            gene_symbols_and_ids = [{'symbol': gene['symbol'], 'hgnc_id': gene['hgnc_id']} for gene in genes]
            message = f"Multiple genes found for symbol `{gene_input}`. Please specify. Found genes: " + \
                      ", ".join([f"{gene['symbol']} (HGNC ID: {gene['hgnc_id']})" for gene in gene_symbols_and_ids])
            return JsonResponse({'error': True, 'message': message}, status=400)
        
        gene = hgnc_data['response']['docs'][0]
        hgnc_id = gene.get('hgnc_id', '')
        symbol = gene.get('symbol', '')
        return JsonResponse({'error': False, 'hgnc_id': hgnc_id, 'symbol': symbol})

    else:
        return JsonResponse({'error': True, 'message': 'Invalid type input'})


def verify_metabolite(request):
    """ 
    Verifies if a metabolite input from the user is in VMH.
    """
    main_input = request.POST.get('metabolite')
    input_type = request.POST.get('type')
    if main_input.strip() == '' and input_type in ['VMH', 'SwissLipids','ChEBI ID', 'ChEBI Name','PubChem ID']:
        return JsonResponse({'error': True, 'message': 'No input provided'})
    if input_type == 'VMH':
        BASE_URL = 'https://www.vmh.life/'
        endpoint = f"{BASE_URL}_api/metabolites/?abbreviation={main_input}"
        response = requests.get(endpoint, verify=False)
        if response.status_code == 200:
            data = response.json()
            res = data.get('results', [])
            if res:
                inchi_string, smile = res[0].get('inchiString', ''), res[0].get('smile', '')
                if not smile and not inchi_string:
                    return JsonResponse({'error': True, 'message': f'Metabolite {main_input} does not have SMILES or inchi String on VMH'}, status=404)
                else:
                    return JsonResponse({'found': True, 'abbr': data['results'][0]['abbreviation'], 'name': data['results'][0]['fullName'], 'miriam': data['results'][0]['miriam'],'input_type':input_type})
            else:
                return JsonResponse({'error': True, 'message': f'Metabolite `{main_input}` not found in VMH'}, status=404)
        else:
            return JsonResponse({'error': True, 'message': f'VMH API returned error {response.status_code} for metabolite `{main_input}`'}, status=500)
    else:
        mols, errors,names = any_to_mol([main_input], [input_type], request, side=None)
        mol, error,name = mols[0], errors[0],names[0]
        if error:
            return JsonResponse({'error': True, 'message': error})
        found, miriam,abbr, nameVMH = search_vmh(mol, return_abbr=True, return_name=True)
        name = nameVMH if found else name
        name = capitalize_first_letter(name)
        return JsonResponse({'found': found, 'abbr': abbr, 'name': name, 'miriam': miriam,'input_type':input_type})   

def get_user_created_reactions(user_id):

    user = User.objects.get(pk=user_id)
    reactions = CreatedReaction.objects.filter(user=user)
    reaction_ids = [reaction.reaction.id for reaction in reactions]
    return reaction_ids

def input_reaction(request):
    """
    Handles the POST request for a reaction input form.
    Processes the reaction data, RDT, Checks for the reaction in VMH,
    and returns the processed data.

    Input:
    - request: The HTTP request object from Django.

    Output:
    - HttpResponse: Renders a template with the reaction form or returns a JsonResponse with reaction data.
    """
    if request.method == 'POST':
        action = request.POST.get('action')  # Action to perform (either 'create' or 'edit')
        form = ReactionForm(request.POST, request.FILES)
        user_id = request.POST.get('userID')
        if action == 'edit':
            reaction_id = request.POST.get('reaction_id')
            try:
                reaction = Reaction.objects.get(id=reaction_id)
            except Reaction.DoesNotExist:
                return JsonResponse({'message': 'Reaction not found', 'status': 'error'})
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

        subs_sch = request.POST.getlist('subs_sch')  # Stoichiometry for substrates
        prod_sch = request.POST.getlist('prod_sch')  # Stoichiometry for products
        subs_comp = request.POST.getlist('subs_comps')  # Compartments for substrates
        prod_comp = request.POST.getlist('prod_comps')
        substrates_types = request.POST.getlist('substrates_type')
        products_types = request.POST.getlist('products_type')
        direction = request.POST.get('direction')
        subs_sch = [int(s) for s in subs_sch]
        prod_sch = [int(s) for s in prod_sch]

        subs_mols, subs_errors, _ = any_to_mol(substrates_list, substrates_types, request, side='substrates')
        prod_mols, prod_errors, _ = any_to_mol(products_list, products_types, request, side='products')
        subsystem = request.POST.get('subsystem')

        all_errors = subs_errors + prod_errors
        if any(elem is not None for elem in all_errors):
            print(all_errors)

            error_message = "\n".join([error for error in all_errors if error is not None])
            if request.headers.get('X-Requested-With') == 'XMLHttpRequest':
                return JsonResponse({'status': 'error', 'message': error_message})
            else:
                # Return error message in context for non-AJAX requests
                context = {'form': form, 'error_message': error_message}
                return render(request, 'reactions/Home_page.html', context)

        metabolite_formulas, metabolite_charges, metabolite_mol_file_strings = get_mol_info(subs_mols + prod_mols)
        metabolite_names = substrates_names + products_names
        reaction_rxn_file = construct_reaction_rxnfile(subs_mols, subs_sch, prod_mols, prod_sch, substrates_names, products_names)
        reaction.save()

        # Skip atom mapping if any product fields are empty or only one substrate is provided
        skip_atom_mapping = request.POST.get('skipAtomMapping') == 'true' or (len(substrates_list) == 1 and len(products_list) == 0)
        if skip_atom_mapping:
            response_data = {'visualizations': ['/images/atom_mapping_skip.png']}
            balanced_count, (subs_atoms, prods_atoms), balanced_charge, (subs_charge, prods_charge), molc_formula, symb_to_name = get_reaction_info(reaction_rxn_file, direction)
        else:
            response_data = RDT(reaction_rxn_file, destination_path_png=f'media/images/visual{reaction.id}.png', destination_path_rxn=f'media/rxn_files/rxn{reaction.id}.rxn')
            balanced_count, (subs_atoms, prods_atoms), balanced_charge, (subs_charge, prods_charge), molc_formula, symb_to_name = get_reaction_info(f'media/rxn_files/rxn{reaction.id}.rxn', direction)

        vmh_found = check_reaction_vmh(substrates_list, products_list, subs_sch, prod_sch, substrates_types, products_types, subs_mols, prod_mols, direction, subsystem, subs_comp, prod_comp)

        if 'error' in response_data:
            context = {'form': form, 'error_message': response_data['error']}
            return render(request, 'reactions/Home_page.html', context)

        subs_found, subs_miriams = search_metabolites_vmh(substrates_list, substrates_types, request, side='substrates')
        prod_found, prod_miriams = search_metabolites_vmh(products_list, products_types, request, side='products')
        substrates_list = get_fields(request, substrates_list, substrates_types, settings.MEDIA_ROOT, settings.MEDIA_URL, side='substrates')
        products_list = get_fields(request, products_list, products_types, settings.MEDIA_ROOT, settings.MEDIA_URL, side='products')
        
        # Assign the values directly to the reaction instance
        reaction.Organs=json.dumps(organs)
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
        reaction.vmh_url = json.dumps(vmh_found['url']) if vmh_found['found'] else None
        reaction.vmh_formula = json.dumps(vmh_found['formula']) if vmh_found['found'] else None
        reaction.metabolite_names = json.dumps(metabolite_names)
        reaction.metabolite_formulas = json.dumps(metabolite_formulas)
        reaction.metabolite_charges = json.dumps(metabolite_charges)
        reaction.metabolite_mol_file_strings = json.dumps(metabolite_mol_file_strings)
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
            }
            if vmh_found['found']:
                data['vmh_found'] = vmh_found['found']
                data['vmh_found_similar'] = vmh_found['similar']
                data['vmh_url'] = vmh_found['url']
                data['vmh_formula'] = vmh_found['formula']
            data['status'] = 'success'
            return JsonResponse(data)
    else:
        form = ReactionForm()
    return render(request, 'reactions/Home_page.html', {'form': form})



def safe_json_loads(data):
    return json.loads(data) if data is not None else None


@require_POST
def clone_reaction_view(request):
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
        return JsonResponse({'status': 'success', 'message': 'Reaction cloned successfully'})
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)})





def get_reaction(request, reaction_id):
    """
    Fetches and returns details of a reaction by its ID.
    
    :param request: The HTTP request object.
    :param reaction_id: The ID of the reaction to retrieve.
    :return: JsonResponse containing the reaction details or an error message.
    """

    try:
        
        reaction = Reaction.objects.get(pk=reaction_id)
        reaction_data = {
            'Organs': reaction.Organs,
            'reaction_id': reaction.id,
            'short_name': reaction.short_name,
            'substrates': safe_json_loads(reaction.substrates),
            'products': safe_json_loads(reaction.products),
            'substrates_names': safe_json_loads(reaction.substrates_names),
            'products_names': safe_json_loads(reaction.products_names),
            'direction': reaction.direction,
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
        }

        return JsonResponse(reaction_data)

    except Reaction.DoesNotExist:
        return JsonResponse({'error': 'Reaction not found'}, status=404)
    
    except Exception as e:
        return JsonResponse({'error': 'An unexpected error occurred'}, status=500)


@csrf_exempt
def add_info_to_reaction(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('userID')
            info_type = data.get('infoType')
            info_text = data.get('infoText')
            ext_link_type = data.get('extLinkType', '')
            ref_type = data.get('refType', '')
            reaction_id = data.get('reactionId')


            if not all([user_id, info_type, info_text]) or (info_type != 'Gene Info' and not reaction_id):
                return JsonResponse({'status': 'error', 'message': 'All fields are required.', 'info_type': info_type}, status=400)
            
            user = User.objects.get(pk=user_id)
            username = user.name

            if info_type == 'Gene Info' and not reaction_id:
                # Store gene info in session if no reaction_id is provided
                if 'gene_info' not in request.session:
                    request.session['gene_info'] = []
                
                info_data = {'info': info_text, 'user_name': username}
                request.session['gene_info'].append(info_data)
                request.session.modified = True

                return JsonResponse({'status': 'success', 'message': 'Gene information added to session.', 'info_type': info_type})

            # Fetch the reaction if reaction_id is provided
            reaction = Reaction.objects.get(id=reaction_id)

            if info_type == 'Reference':
                if reaction.references is None:
                    reaction.references = []
                info_data_template = {'user_name': username, 'ref_type': ref_type}

                # Determine the delimiter based on the reference type
                if 'PMID' in ref_type:
                    ref_list = [ref.strip() for ref in info_text.split(';')]
                    existing_refs = {ref['info'] for ref in reaction.references if 'PMID' in ref['info']}
                elif 'DOI' in ref_type:
                    ref_list = [ref.strip() for ref in info_text.split(',')]
                    existing_refs = {ref['info'] for ref in reaction.references if 'DOI' in ref['info']}
                else:
                    ref_list = [info_text.strip()]
                    existing_refs = {ref['info'] for ref in reaction.references}

                for ref in ref_list:
                    if ref not in existing_refs:
                        info_data = info_data_template.copy()
                        info_data['info'] = ref
                        reaction.references.append(info_data)
                        existing_refs.add(ref)

            elif info_type == 'External Link':
                if reaction.ext_links is None:
                    reaction.ext_links = []
                info_data = {'info': info_text, 'user_name': username, 'ext_link_type': ext_link_type}
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
            return JsonResponse({'status': 'success', 'message': 'Information added successfully.', 'reaction_id': reaction_id, 'info_type': info_type})

        except User.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Invalid user key.', 'info_type': info_type}, status=404)
        except Reaction.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Invalid reaction ID.', 'info_type': info_type}, status=404)
        except json.JSONDecodeError:
            return JsonResponse({'status': 'error', 'message': 'Invalid JSON.', 'info_type': info_type}, status=400)
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e), 'info_type': info_type}, status=500)
    else:
        return JsonResponse({'status': 'error', 'message': 'Invalid request method.', 'info_type': info_type}, status=405)


@csrf_exempt
def update_gene_info(request):
    try:
        data = json.loads(request.body)
        user_id = data.get('userID')
        reaction_id = data.get('reactionID')
        gene = data.get('gene')
        field_type = data.get('fieldType')
        updated_value = data.get('updatedValue')

        # Debugging statements to trace values

        # Check if any required field is missing or updated_value is empty
        if not all([user_id, reaction_id, gene, field_type]) or updated_value is None or updated_value.strip() == "":
            return JsonResponse({'status': 'error', 'message': 'All fields are required and updated value cannot be empty.'}, status=400)

        user = User.objects.get(pk=user_id)
        reaction = Reaction.objects.get(id=reaction_id)

        if reaction.gene_info is None:
            return JsonResponse({'status': 'error', 'message': 'No gene information found.'}, status=404)

        gene_info_updated = False
        for gene_info in reaction.gene_info:
            if gene in gene_info['info']:
                if field_type == 'Organs':
                    organ_match = re.search(r'ORGAN\(([^)]+)\)', gene_info['info'])
                    if organ_match:
                        old_organs = organ_match.group(1)
                        new_info = gene_info['info'].replace(f'ORGAN({old_organs})', f'ORGAN({updated_value})')
                        gene_info['info'] = new_info
                        gene_info_updated = True
                elif field_type == 'SubcellularLocations':
                    subcellular_match = re.search(r'SUBCELLULAR\(([^)]+)\)', gene_info['info'])
                    if subcellular_match:
                        old_subcellular = subcellular_match.group(1)
                        new_info = gene_info['info'].replace(f'SUBCELLULAR({old_subcellular})', f'SUBCELLULAR({updated_value})')
                        gene_info['info'] = new_info
                        gene_info_updated = True

        if gene_info_updated:
            reaction.save()
            return JsonResponse({'status': 'success', 'message': 'Gene information updated successfully.'})
        else:
            return JsonResponse({'status': 'error', 'message': 'Gene information not found or not updated.'}, status=404)

    except User.DoesNotExist:
        return JsonResponse({'status': 'error', 'message': 'Invalid user key.'}, status=404)
    except Reaction.DoesNotExist:
        return JsonResponse({'status': 'error', 'message': 'Invalid reaction ID.'}, status=404)
    except json.JSONDecodeError:
        return JsonResponse({'status': 'error', 'message': 'Invalid JSON.'}, status=400)
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)}, status=500)    
    
@csrf_exempt
def delete_reaction_info(request):
    if request.method == 'POST':
        try:
            req_body = json.loads(request.body)
            reaction_id = req_body.get('reaction_id')
            tab_id = req_body.get('tab_id')
            item_to_delete = req_body.get('item_to_delete')
            # Fetch the reaction object
            react_obj = Reaction.objects.get(pk=reaction_id)

            # Handle deletion based on the tab_id
            if tab_id == 'refs-content':
                for ref in react_obj.references:
                    if ref['info'] == item_to_delete['info'] and ref['ref_type'] == item_to_delete['ref_type']:
                        react_obj.references.remove(ref)
                        react_obj.save()
                        return JsonResponse({'status': 'success', 'message': 'Reference deleted successfully'})

            elif tab_id == 'ext-links-content':
                for ext_link in react_obj.ext_links:
                    if ext_link['info'] == item_to_delete['info'] and ext_link['ext_link_type'] == item_to_delete['ext_link_type']:
                        react_obj.ext_links.remove(ext_link)
                        react_obj.save()
                        return JsonResponse({'status': 'success', 'message': 'External link deleted successfully'})

            elif tab_id == 'comments-content':
                for comment in react_obj.comments:
                    if comment['info'] == item_to_delete['info']:
                        react_obj.comments.remove(comment)
                        react_obj.save()
                        return JsonResponse({'status': 'success', 'message': 'Comment deleted successfully'})

            elif tab_id == 'gene-info-content':
                for gene_info in react_obj.gene_info:
                    if gene_info['info'] == item_to_delete['info']:
                        react_obj.gene_info.remove(gene_info)
                        react_obj.save()
                        return JsonResponse({'status': 'success', 'message': 'Gene info deleted successfully'})

            # If no valid tab_id matched
            return JsonResponse({'error': True, 'message': 'Invalid tab_id or item not found'}, status=400)
        
        except Reaction.DoesNotExist:
            return JsonResponse({'error': True, 'message': 'Reaction not found'}, status=404)
        except Exception as e:
            return JsonResponse({'error': True, 'message': str(e)}, status=500)
    else:
        return JsonResponse({'error': True, 'message': 'Invalid request method'}, status=400)


def check_session_data(request):
    # Check if 'gene_info' is in session
    gene_info = request.session.get('gene_info', None)
    
    if gene_info:
        return JsonResponse({'status': 'success', 'gene_info': gene_info})
    else:
        return JsonResponse({'status': 'error', 'message': 'No gene info in session.'})
    
@csrf_exempt  # Use csrf_exempt if you don't want to deal with CSRF tokens in development. However, be cautious with this in production.
def delete_gene_info_from_session(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            info_to_delete = data.get('info_to_delete')
            
            # Get the gene info list from the session
            gene_info = request.session.get('gene_info', [])
            
            # Filter out the item to delete
            updated_gene_info = [info for info in gene_info if info['info'] != info_to_delete]
            
            # Update the session with the filtered list
            request.session['gene_info'] = updated_gene_info
            request.session.modified = True
            
            return JsonResponse({'status': 'success', 'message': 'Gene info deleted from session.'})
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e)}, status=400)
    else:
        return JsonResponse({'status': 'error', 'message': 'Invalid request method.'}, status=405)


@csrf_exempt  # Use csrf_exempt if you don't want to deal with CSRF tokens in development. However, be cautious with this in production.
def clear_session(request):
    if request.method == 'POST':
        try:
            # Clear the entire session
            request.session.flush()  # This will clear all session data

            return JsonResponse({'status': 'success', 'message': 'Session cleared successfully.'})
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e)}, status=400)
    else:
        return JsonResponse({'status': 'error', 'message': 'Invalid request method.'}, status=405)

def get_reaction_details(request):
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
    

def validate_user_ID(user_id):
    """Validate the user ID."""
    try:
        user = User.objects.get(pk=user_id)
        return user
    except User.DoesNotExist:
        return None
    
def save_user_reaction(request):
    if request.method == 'POST':
        reaction_id = request.POST.get('reaction_id')
        userID = request.POST.get('userID')
        short_name = request.POST.get('short_name')
        flag_name = request.POST.get('flag_name')
        flag_color = request.POST.get('flag_color')


        try:
            reaction = Reaction.objects.get(pk=reaction_id)
        except Reaction.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Invalid reaction'})

        user = validate_user_ID(userID)

        if user and reaction:
            reaction.short_name = short_name
            if flag_name.strip() not in ['None', 'Choose a flag'] and flag_color != 'null':
                # Get or create the flag
                flag = Flag.objects.get(name_flag=flag_name, color=flag_color)

                # Associate the new flag with the reaction
                reaction.flags.add(flag)

            reaction.save()
            user.saved_reactions.add(reaction)
            return JsonResponse({'status': 'success'})

        return JsonResponse({'status': 'error', 'message': 'Invalid user or reaction'})

    return JsonResponse({'status': 'error', 'message': 'Invalid request'})


@csrf_exempt 
def save_flags_in_saved_reactions(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('userID')
            reaction_ids = data.get('reaction_ids', [])
            flag_name = data.get('flag_name')
            flag_color = data.get('flag_color')

            user = User.objects.get(pk=user_id)
            flag = Flag.objects.get(name_flag=flag_name, color=flag_color, user=user)
            for reaction_id in reaction_ids:
                reaction = get_object_or_404(Reaction, pk=reaction_id)
                reaction.flags.add(flag)

            return JsonResponse({'status': 'success'})
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e)})
    return JsonResponse({'status': 'error', 'message': 'Invalid request method'})

def saved_reactions(request, modal=False):
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
            
            # Construct details strings
            subs_details = " + ".join(f"{float(sch)} {json.loads(reaction.substrates_names)[idx]} [{comp}]" for idx, (sch, comp) in enumerate(zip(subs_sch, subs_comps)))
            prods_details = " + ".join(f"{float(sch)} {json.loads(reaction.products_names)[idx]} [{comp}]" for idx, (sch, comp) in enumerate(zip(prods_sch, prods_comps)))

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
            flag_details = [{"name": flag.name_flag, "color": flag.color} for flag in flags]

            combined_reactions_details.append({
                'reaction': reaction,
                'details': {
                    'short_name': short_name,
                    'subs_details': subs_details,
                    'prods_details': prods_details,
                    'molc_formula': reaction.molc_formula,
                    'balanced_count': json.loads(reaction.balanced_count)[0] if reaction.balanced_count else None,
                    'balanced_charge': json.loads(reaction.balanced_charge)[0] if reaction.balanced_charge else None,
                    'subsystem': reaction.subsystem,
                    'direction': reaction.direction,
                    'gene_info': gene_info_list,
                    'flags': flag_details,  # Include flag details with name and color
                }
            })
        
        context = {
            'reactions': reactions, 
            'reactions_json': reactions_json, 
            'userID': userID, 
            'user_name': user.name, 
            'combined_reactions_details': combined_reactions_details
        }

        if modal:
            return render(request, 'reactions/saved_reactions_modal.html', context)
        else:
            return render(request, 'reactions/saved_reactions.html', context)
    else:
        return render(request, 'reactions/error.html', {'error_message': 'Invalid key'})


def get_ai_response(request):
    if request.method == 'POST':
        userID = request.session.get('userID')
        user = validate_user_ID(userID)
        if not user.cred_add_to_vmh:
            return JsonResponse({'status': 'error', 'error_message': 'User does not have permission to access this feature'}, status=403)
        if user:
            try:
                data = json.loads(request.body)  # Parse the JSON data from the request body
                input_text = data.get('key')  # This should match the key you send in the JSON
                temperature = data.get('temperature', 0.5)  # Get temperature value, default to 0.5 if not provided
                genes = input_text.split(' ')
                genes = [gene for gene in genes if gene not in ['and', 'or', 'AND', 'OR']]
                predictions = {}
                for gene in genes:
                    llm_response_html = get_gpt_predictions(gene, temperature)
                    predictions[gene] = llm_response_html
                return JsonResponse({'status': 'success', 'predictions': predictions})
            except json.JSONDecodeError:
                return JsonResponse({'status': 'error', 'error_message': 'Invalid JSON'}, status=400)
        else:
            return render(request, 'reactions/error.html', {'error_message': 'Invalid key'})
    else:
        return JsonResponse({'status': 'error', 'error_message': 'Unsupported method'}, status=405)

    
def get_gpt_predictions(gene, temperature):
    gpt_predictions_names = []  # Parsed GPT output into separate metabolite and gene reactions
    reactions = askGPT4(gene, temperature)  # Function that asks ChatGPT to predict the reactions associated with the gene
    parsed_reactions = parse_metabolic_reactions_gpt(reactions)  # Function to parse ChatGPT's answer
    predicted_reaction_names = [reaction for reaction in parsed_reactions]  # Collect parsed reactions
    gpt_predictions_names.append(predicted_reaction_names)
    vmh_abbreviations = [
        [get_vmh_met_from_inchi(metanetx_to_inchi(met), met) for met in reaction]
        for reaction in gpt_predictions_names[0]
    ]
    return json.dumps(vmh_abbreviations)
   # return JsonResponse({'status': 'success', 'data': vmh_abbreviations})
    
    
def delete_reaction(request):
    if request.method == 'POST':
        userID = request.POST.get('userID')
        reaction_id = request.POST.get('reaction_id')
        user = validate_user_ID(userID)
        if user:
            reaction = Reaction.objects.get(pk=reaction_id)
            user.saved_reactions.remove(reaction)
            saved_reactions_url = reverse('saved_reactions')
            return redirect(saved_reactions_url)
    return render(request, '404.html')

logger = logging.getLogger(__name__)

@csrf_exempt
def prepare_add_to_vmh(request):
    if request.method != 'POST':
        return JsonResponse({'status': 'error', 'message': 'Invalid request method. Use POST instead.'}, status=400)

    try:
        request_data = json.loads(request.body)
        reactionIds = request_data['reactionIds']
    except json.JSONDecodeError:
        logger.error("Invalid JSON data.")
        return JsonResponse({'status': 'error', 'message': 'Invalid JSON data.'}, status=400)
    except KeyError:
        logger.error("Missing reactionIds in the request data.")
        return JsonResponse({'status': 'error', 'message': 'Missing reactionIds in the request data.'}, status=400)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return JsonResponse({'status': 'error', 'message': 'An unexpected error occurred.'}, status=500)

    try:
        reaction_objs = [Reaction.objects.get(pk=int(reactionId)) for reactionId in reactionIds]

        in_vmh = [reaction.vmh_found and not reaction.vmh_found_similar for reaction in reaction_objs]
        if True in in_vmh:
            reaction_objs_in_vmh = [reaction for reaction, found in zip(reaction_objs, in_vmh) if found]
            names = [reaction.short_name for reaction in reaction_objs_in_vmh]
            return JsonResponse({'status': 'error', 'message': f'The following reactions are already in VMH: {", ".join(names)}'})

        subs_in_vmh = [json.loads(reaction.subs_found) if reaction.subs_found else [] for reaction in reaction_objs]
        prods_in_vmh = [json.loads(reaction.prod_found) if reaction.prod_found else [] for reaction in reaction_objs]


        matlab_session = MatlabSessionManager()
        subs_abbr = [[gen_metabolite_abbr(sub, sub_type, sub_name, search_metabolites_vmh, matlab_session) for sub, sub_type, sub_name in zip(json.loads(reaction.substrates), json.loads(reaction.substrates_types), json.loads(reaction.substrates_names))] for reaction in reaction_objs]
        prods_abbr = [[gen_metabolite_abbr(prod, prod_type, prod_name, search_metabolites_vmh, matlab_session) for prod, prod_type, prod_name in zip(json.loads(reaction.products), json.loads(reaction.products_types), json.loads(reaction.products_names))] for reaction in reaction_objs]

        subs_need_new_names = [[] for _ in reactionIds]
        prods_need_new_names = [[] for _ in reactionIds]

        for idx, in_vmh_list in enumerate(subs_in_vmh):
            for j, sub_in_vmh in enumerate(in_vmh_list):
                if not sub_in_vmh:
                    need_new_name = is_name_in_vmh(json.loads(reaction_objs[idx].substrates_names)[j])
                    subs_need_new_names[idx].append(need_new_name)
                else:
                    subs_need_new_names[idx].append(False)
        for idx, in_vmh_list in enumerate(prods_in_vmh):
            for j, prod_in_vmh in enumerate(in_vmh_list):
                if prod_in_vmh:
                    prods_need_new_names[idx].append(False)
                else:
                    need_new_name = is_name_in_vmh(json.loads(reaction_objs[idx].products_names)[j])
                    prods_need_new_names[idx].append(need_new_name)

        reaction_abbrs = [gen_reaction_abbr(sub_abbr, prod_abbr, reaction) for sub_abbr, prod_abbr, reaction in zip(subs_abbr, prods_abbr, reaction_objs)]

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
        logger.error(f"Error processing reactions: {e}")
        return JsonResponse({'status': 'error', 'message': 'An error occurred while processing reactions.'}, status=500)
    
@csrf_exempt  # Remove this if you want to deal with CSRF tokens later
def create_formula_abbr(request):
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
        if not all([metabolite, mtype, metabolite_name]):
            return JsonResponse({'error': 'Missing data'}, status=400)

        # Initialize Matlab session manager
        matlab_session = MatlabSessionManager()

        # Use the bypass function to force the else clause
        abbr = gen_metabolite_abbr(metabolite, mtype, metabolite_name, bypass_search_func, matlab_session)

        # Return the abbreviation as JSON
        return JsonResponse({'abbr': abbr})

    return JsonResponse({'error': 'Invalid request method'}, status=400)


def bypass_search_func(metabolites, types, *args, **kwargs):
    # Always return False for found and None for abbreviation
    return [False], [None]

def add_to_vmh(request):
    """
    Main function to handle the request for adding reactions to VMH.
    """
    req_body = json.loads(request.body)

    userID = req_body.get('userID')
    user = User.objects.get(pk=userID)
    if not user:
        return JsonResponse({'status': 'error', 'message': 'Invalid user key.'}, status=404)
    if not user.cred_add_to_vmh:
        return JsonResponse({'status': 'error', 'message': 'User does not have permission to add to VMH.'}, status=403)
    user_name = user.name
    user_full_name = user.name
    reactions = req_body.get('reactions')
    reaction_ids = []
    not_enough_info, no_comments, not_balanced = [], [], []
    met_added_info = {}
    names_list = [reaction['short_name'] for reaction in reactions]
    abbr_list = [reaction['abbreviation'] for reaction in reactions]

    # Check for abbreviations, names, and confidence scores
    missing_names = [reaction['short_name'] == '' for reaction in reactions]
    missing_abbrs = [reaction['abbreviation'] == '' for reaction in reactions]
    missing_conf_scores = [reaction['confidence_score'] == '" "' for reaction in reactions]

    if True in missing_names:
        return JsonResponse({'status': 'error', 'message': 'Please enter a description for reaction'})

    if True in missing_abbrs:
        return JsonResponse({'status': 'error', 'message': 'Please enter an abbreviation'})

    if True in missing_conf_scores:
        return JsonResponse({'status': 'error', 'message': 'Please enter a confidence score for all reactions'})

    for name in names_list:
        if names_list.count(name) > 1:
            return JsonResponse({'status': 'error', 'message': f'Reaction with name `{name}` is repeated in the list.'})
    for abbr in abbr_list:
        if abbr_list.count(abbr) > 1:
            return JsonResponse({'status': 'error', 'message': f'Reaction with abbreviation `{abbr}` is repeated in the list.'})

    name_in_vmh, abbr_in_vmh = check_names_abbrs_vmh([(reaction['short_name'], reaction['abbreviation']) for reaction in reactions])
    if True in list(name_in_vmh.values()):
        name_in_vmh_reactions = [reaction for reaction, in_vmh in zip(reactions, name_in_vmh.values()) if in_vmh]
        return JsonResponse({'status': 'error', 'message': f'The following reaction descriptions are already in VMH: {", ".join([reaction["short_name"] for reaction in name_in_vmh_reactions])}'})
    if True in list(abbr_in_vmh.values()):
        abbr_in_vmh_reactions = [reaction for reaction, in_vmh in zip(reactions, abbr_in_vmh.values()) if in_vmh]
        return JsonResponse({'status': 'error', 'message': f'The following reaction abbreviations are already in VMH: {", ".join([reaction["abbreviation"] for reaction in abbr_in_vmh_reactions])}'})

    reactions_new_subsInfo = [json.loads(reaction['substrates_info']) for reaction in reactions]
    reactions_new_prodsInfo = [json.loads(reaction['products_info']) for reaction in reactions]
    reactions_subs_found = [json.loads(Reaction.objects.get(pk=reaction['pk']).subs_found) for reaction in reactions]
    reactions_prods_found = [json.loads(Reaction.objects.get(pk=reaction['pk']).prod_found) for reaction in reactions]
    subs_names_vmh, subs_abbr_vmh, prods_names_vmh, prods_abbr_vmh = check_met_names_abbrs_vmh(reactions_new_subsInfo, reactions_new_prodsInfo, reactions_subs_found, reactions_prods_found)
    if True in list(subs_names_vmh.values()):
        subs_names_in_vmh = [sub for sub in subs_names_vmh.keys() if subs_names_vmh[sub]]
        return JsonResponse({'status': 'error', 'message': f'The following substrates have metabolite names that are already in VMH: {", ".join(subs_names_in_vmh)}'})
    if True in list(subs_abbr_vmh.values()):
        subs_abbrs_in_vmh = [sub for sub in subs_abbr_vmh.keys() if subs_abbr_vmh[sub]]
        return JsonResponse({'status': 'error', 'message': f'The following substrates have metabolite abbreviations that are already in VMH: {", ".join(subs_abbrs_in_vmh)}'})
    if True in list(prods_names_vmh.values()):
        prods_names_in_vmh = [prod for prod in prods_names_vmh.keys() if prods_names_vmh[prod]]
        return JsonResponse({'status': 'error', 'message': f'The following products have metabolite names that are already in VMH: {", ".join(prods_names_in_vmh)}'})
    if True in list(prods_abbr_vmh.values()):
        prods_abbrs_in_vmh = [prod for prod in prods_abbr_vmh.keys() if prods_abbr_vmh[prod]]
        return JsonResponse({'status': 'error', 'message': f'The following products have metabolite abbreviations that are already in VMH: {", ".join(prods_abbrs_in_vmh)}'})
    subs_abbr = []
    prods_abbr = []
    for reaction in reactions:
        
        obj = Reaction.objects.get(pk=reaction['pk'])
        obj.short_name = reaction['short_name']
        reaction_ids.append(obj.id)
        # Update substrate names and abbreviations
        subs_info = json.loads(reaction['substrates_info'])
        new_subs_names = [capitalize_first_letter(sub['name']) for sub in subs_info]
        new_subs_abbrs = [sub['abbreviation'] for sub in subs_info]
        subs_abbr.append(new_subs_abbrs)
        obj.substrates_names = json.dumps(new_subs_names)
        # Update product names and abbreviations
        prods_info = json.loads(reaction['products_info'])
        new_prods_names = [capitalize_first_letter(prod['name']) for prod in prods_info]
        new_prods_abbrs = [prod['abbreviation'] for prod in prods_info]
        prods_abbr.append(new_prods_abbrs)
        obj.products_names = json.dumps(new_prods_names)
        # Update references, external links, and comments
        references, ext_links, comments = json.loads(reaction['references']), json.loads(reaction['ext_links']), json.loads(reaction['comments'])
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
        if f"Created and Added to VMH by: {user_full_name}" not in list(map(lambda x: x['info'], new_comments)):
            new_comments.append({'info': f"Created and Added to VMH by: {user_full_name}", 'user_name': user_name})
        not_enough_info.append(len(new_references) < 1 and len(new_ext_links) < 1)
        no_comments.append(len(new_comments) < 2)
        obj.references = new_references if new_references else None
        obj.ext_links = new_ext_links if new_ext_links else None
        obj.comments = new_comments if new_comments else None
        obj.confidence_score = reaction.get('confidence_score', 0)  # Add confidence score
        if not (json.loads(obj.balanced_charge)[0] and json.loads(obj.balanced_count)[0]):
            not_balanced.append(True)
        else:
            not_balanced.append(False)
        # Save the updated reaction object
        obj.save()
    reaction_objs = [Reaction.objects.get(pk=reaction_id) for reaction_id in reaction_ids]

    if True in not_enough_info:
        not_enough_info_reactions = [reaction for reaction, not_enough in zip(reaction_objs, not_enough_info) if not_enough]
        return JsonResponse({'status': 'error', 'message': f'The following reactions do not have at least one reference or external link: {", ".join([reaction.short_name for reaction in not_enough_info_reactions])}'})
    if True in no_comments:
        no_comments_reactions = [reaction for reaction, no_comment in zip(reaction_objs, no_comments) if no_comment]
        return JsonResponse({'status': 'error', 'message': f'The following reactions do not have at least one comment: {", ".join([reaction.short_name for reaction in no_comments_reactions])}'})
    if True in not_balanced:
        not_balanced_reactions = [reaction for reaction, not_balanced in zip(reaction_objs, not_balanced) if not_balanced]
        return JsonResponse({'status': 'error', 'message': f'The following reactions are not balanced: {", ".join([reaction.short_name for reaction in not_balanced_reactions])}'})
    all_vmh = False if any([any(element != 'VMH' for element in json.loads(reaction.substrates_types)) or any(element != 'VMH' for element in json.loads(reaction.products_types)) for reaction in reaction_objs]) else True
    reaction_identifiers, reaction_names = [reaction['abbreviation'] for reaction in reactions], [reaction.short_name for reaction in reaction_objs]
    matlab_session = MatlabSessionManager()
    if not all_vmh:
        unique_abbrs, unique_mols, unique_types, unique_names = get_nonfound_metabolites(reaction_objs, subs_abbr, prods_abbr, search_func=search_metabolites_vmh)
        if unique_abbrs and unique_mols and unique_types:
            abbrs = unique_abbrs
            smiles, errors = any_to_smiles(unique_mols, unique_types, None, side=None)
            smiles = ['' if v != None else smiles[i] for i, v in enumerate(errors)]
            inchikeys = smiles_to_inchikeys(smiles)
            names = unique_names
            formulas, charges = smiles_to_charged_formula(smiles)
            json_paths = met_prepare_json_paths_and_variables(abbrs, names, formulas, charges, inchikeys, smiles)
            matlab_result = add_metabolites_matlab(json_paths, matlab_session)
            for path in json_paths:
                os.remove(path)
            if matlab_result['status'] == 'success':
                met_ids = matlab_result['met_ids']
                met_added_info = {abbr: [met_id, formula, inchikey] for abbr, met_id, formula, inchikey in zip(abbrs, met_ids, formulas, inchikeys)}
                for abbr, info in met_added_info.items():
                    MetabolitesAddedVMH.objects.create(
                        user=user,
                        user_name=user_name,
                        metabolite_id=info[0],
                        metabolite_formula=info[1],
                        metabolite_abbr=abbr,
                    )
            else:
                matlab_session.quit()
                return JsonResponse({'status': 'error', 'message': matlab_result['message']})
    reaction_formulas = [construct_vmh_formula(reaction_objs[idx], subs_abbr[idx], prods_abbr[idx]) for idx in range(len(reaction_objs))]
    reaction_directions, reaction_subsystems, reaction_references, reaction_external_links, reaction_gene_info, reaction_comments, reaction_confidence_scores = gather_reaction_details(reaction_objs)
    json_paths = rxn_prepare_json_paths_and_variables(reaction_identifiers, reaction_names, reaction_formulas, reaction_directions, reaction_subsystems, reaction_references, reaction_external_links, reaction_gene_info, reaction_comments, reaction_confidence_scores)
    matlab_result = add_reaction_matlab(json_paths, matlab_session)
    # Cleanup temporary JSON files
    for path in json_paths:
        os.remove(path)

    if matlab_result['status'] == 'success':
        rxn_added_info = {abbr: [rxn_id, reaction_formula] for abbr, rxn_id, reaction_formula in zip(reaction_identifiers, matlab_result['rxn_ids'], reaction_formulas)}
        for abbr, info in rxn_added_info.items():
            ReactionsAddedVMH.objects.create(
                user=user,
                user_name=user_name,
                reaction_id=info[0],
                reaction_formula=info[1],
                reaction_abbr=abbr,
            )
        matlab_session.quit()
        return JsonResponse({'status': 'success', 'rxn_added_info': rxn_added_info, 'met_added_info': met_added_info})
    else:
        matlab_session.quit()
        return JsonResponse({'status': 'error', 'message': matlab_result['message']})

    
@csrf_exempt    
def gene_parsing(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        statement = data.get('geneinfo', '')

        # Patterns for various parts of the logical statement
        alphanumeric_pattern = r'[A-Za-z0-9]+'
        operator_pattern = r'(AND|OR)'
        substatement_pattern = fr'\(\s*{alphanumeric_pattern}\s*{operator_pattern}\s*{alphanumeric_pattern}\s*\)'

        # Full pattern combining the subpatterns
        full_pattern = fr'^{alphanumeric_pattern}\s*({operator_pattern}\s*{alphanumeric_pattern}\s*)*$'

        # Compile the regular expression
        pattern = re.compile(full_pattern)

        # Match the statement against the pattern
        match = pattern.match(statement)

        response_data = {
            'processed_string': statement if match else None,
            'error': None
        }

        # If there's no match, determine why
        if not match:
            if not re.match(alphanumeric_pattern, statement):
                response_data['error'] = 'The statement must start with a Gene.'
            elif not re.search(fr'\s*{operator_pattern}\s*', statement):
                response_data['error'] = 'The statement must contain at least one AND/OR operator after an alphanumeric string.'
            elif re.search(r'[^A-Za-z0-9\s\(\)ANDOR]', statement):
                response_data['error'] = 'The statement contains invalid characters. Only alphanumeric characters, spaces, parentheses, and the words AND/OR are allowed.'
            else:
                response_data['error'] = 'The statement does not match the required logical pattern. Ensure it follows the structure: alphanumeric (AND/OR alphanumeric).'

        return JsonResponse(response_data)
    return JsonResponse({'error': 'Invalid request method'}, status=405)


def parse_genes(data_string):
    # Split by spaces and logical operators, retain only alphanumeric strings
    genes = re.split(r'\s+(?:AND|OR)\s+|\s+', data_string)
    genes = [gene for gene in genes if gene.isalnum()]
    return genes


ORGAN_MAPPING = {
    "adipose tissue": "Adipocytes",
    "adrenal gland": "Agland",
    "amygdala": "Brain",
    "appendix": "/",
    "basal ganglia": "Brain",
    "bone marrow": "/",
    "breast": "Breast",
    "cerebellum": "Brain",
    "cerebral cortex": "Brain",
    "cervix": "Cervix",
    "choroid plexus": "Brain",
    "colon": "Colon",
    "duodenum": "sIEC",
    "endometrium 1": "Uterus",
    "epididymis": "Testis",
    "esophagus": "Lung",
    "fallopian tube": "Uterus",
    "gallbladder": "Gall",
    "heart muscle": "Heart",
    "hippocampal formation": "Brain",
    "hypothalamus": "Brain",
    "kidney": "Kidney",
    "liver": "Liver",
    "lung": "Lung",
    "lymph node": "/",
    "midbrain": "Brain",
    "ovary": "Ovary",
    "pancreas": "Pancreas",
    "parathyroid gland": "Pthyroidgland",
    "pituitary gland": "/",
    "placenta": "/",
    "prostate": "Prostate",
    "rectum": "Colon",
    "retina": "Retina",
    "salivary gland": "/",
    "seminal vesicle": "Testis",
    "skeletal muscle": "Muscle",
    "skin 1": "Skin",
    "small intestine": "sIEC",
    "smooth muscle": "/",
    "spinal cord": "Scord",
    "spleen": "Spleen",
    "stomach 1": "Stomach",
    "testis": "Testis",
    "thymus": "/",
    "thyroid gland": "Thyroidgland",
    "tongue": "/",
    "tonsil": "/",
    "urinary bladder": "Urinarybladder",
    "vagina": "Cervix"
}
# Determine the base directory three levels up from the current file's directory
# Construct the base directory path
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
base_dir = os.path.join(base_dir, 'reconstructor')
# Construct the full path to the config.json file
config_path = os.path.join(base_dir, 'config.json')
# Load the config file
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

# Extract the file_path from the config
file_path = config.get('file_path')

# Construct the full path to the file using base_dir and file_path from the config
full_file_path = os.path.join(base_dir, file_path)

# Load the CSV file using the full file path
df = pd.read_csv(full_file_path, sep='\t')


@csrf_exempt
def gene_details_view(request):
    if request.method != 'POST':
        return JsonResponse({"error": "Invalid request method"}, status=405)

    try:
        # Log the raw request body

        data = json.loads(request.body)
        # Extract the nested infoText from the data object
        data_info = data.get("infoText", {})
        data_string = data_info.get("infoText", "")
        if not isinstance(data_string, str):
            raise TypeError("infoText must be a string")

    except json.JSONDecodeError:
        return JsonResponse({"error": "Invalid JSON"}, status=400)
    except TypeError as e:
        return JsonResponse({"error": str(e)}, status=400)

    genes = parse_genes(data_string)
    
    gene_details = []
    info_text_parts = [f"GPR: {data_string}"]

    all_organs = set()
    all_subcellular_locations = set()

    for gene in genes:
        if not gene.isalnum():  # Skip if gene is not alphanumeric
            continue

        gene_info = {}
        # Fetch and map gene expression data
        unique_organs, error = fetch_and_map_gene_expression(gene, df, ORGAN_MAPPING)
        if error:
            gene_info["ORGAN"] = "Error fetching organs"
        else:
            gene_info["ORGAN"] = unique_organs
            all_organs.add(unique_organs)

        # Get subcellular locations
        subcellular_locations = get_subcellular_locations(gene)
        if subcellular_locations:
            mapped_locations = map_locations_to_wbm(subcellular_locations)
            gene_info["SUBCELLULAR LOCATION"] = mapped_locations
            mapped_locations = extract_unique_elements(mapped_locations)
            all_subcellular_locations.update(mapped_locations)
        else:
            gene_info["SUBCELLULAR LOCATION"] = "Subcellular locations not found"

        gene_details.append(gene_info)
        
    all_subcellular_locations = list(set(all_subcellular_locations))

    # Combine and deduplicate organ and subcellular location information
    combined_organs_text = ", ".join(all_organs)
    combined_subcellular_text = ", ".join(all_subcellular_locations)

    # Format the final infoText part
    info_text_parts.append(f"ORGAN({combined_organs_text}), SUBCELLULAR({combined_subcellular_text})")

    organized_result = {
        "userID": data_info.get("userID", ""),
        "infoType": data_info.get("infoType", ""),
        "extLinkType": data_info.get("extLinkType", ""),
        "refType": data_info.get("refType", ""),
        "reactionId": data_info.get("reactionId", ""),
        "infoText": "; ".join(info_text_parts)
    }

    return JsonResponse(organized_result)


def extract_unique_elements(input_set):
    unique_elements = set()
    for item in input_set:
        elements = item.split(',')
        unique_elements.update(elements)
    return unique_elements

def parse_gene_info(request):
    info = request.GET.get('info', '')

    if not info:
        return JsonResponse({'error': 'No info provided'}, status=400)

    response_data = {}

    try:
        # Split the input by `; GENE:` to separate each gene's information
        gene_sections = re.split(r';\s*GENE:', info)

        # The first part contains GPR, we'll ignore it as we focus on genes
        for section in gene_sections:
            gene_match = re.match(r'([^;]+); ORGAN\(([^)]+)\), SUBCELLULAR\(([^)]+)\)', section)

            if not gene_match:
                continue

            gene = gene_match.group(1).strip()
            organs = [organ.strip() for organ in gene_match.group(2).split(',')]
            subcellular_locations = [loc.strip() for loc in gene_match.group(3).split(',')]

            if gene not in response_data:
                response_data[gene] = {'Organs': [], 'SubcellularLocations': []}

            response_data[gene]['Organs'].extend(organs)
            response_data[gene]['SubcellularLocations'].extend(subcellular_locations)
    
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)
    return JsonResponse(response_data)

@csrf_exempt
def temp_gene_details(request):
    if request.method != 'POST':
        return JsonResponse({"error": "Invalid request method"}, status=405)

    try:
        # Log the raw request body

        data = json.loads(request.body)
        genes_string = data.get("genes", "")

        # Parse genes using the parse_genes function
        genes = parse_genes(genes_string)

        gene_details = []
        info_text_parts = [f"GPR: {genes_string}"]

        all_organs = set()
        all_subcellular_locations = set()

        for gene in genes:
            if not gene.isalnum():
                continue

            gene_info = {}
            unique_organs, error = fetch_and_map_gene_expression(gene, df, ORGAN_MAPPING)
            if error:
                gene_info["ORGAN"] = "Error fetching organs"
            else:
                gene_info["ORGAN"] = unique_organs
                all_organs.update(unique_organs)

            subcellular_locations = get_subcellular_locations(gene)
            if subcellular_locations:
                mapped_locations = map_locations_to_wbm(subcellular_locations)
                gene_info["SUBCELLULAR LOCATION"] = mapped_locations
                mapped_locations = extract_unique_elements(mapped_locations)
                all_subcellular_locations.update(mapped_locations)
            else:
                gene_info["SUBCELLULAR LOCATION"] = "Subcellular locations not found"

            gene_details.append(gene_info)

        response_data = {
            "gene_details": gene_details,
            "all_subcellular_locations": list(all_subcellular_locations),
            "all_organs": list(all_organs)
        }
        return JsonResponse(response_data)

    except json.JSONDecodeError:
        return JsonResponse({"error": "Invalid JSON"}, status=400)
    except Exception as e:
        return JsonResponse({"error": str(e)}, status=500)
        
def get_all_reaction_ids(request):
    """
    Fetches and returns all reaction IDs.
    """
    reaction_ids = Reaction.objects.values_list('id', flat=True)
    return JsonResponse(list(reaction_ids), safe=False)



def get_user_saved_reaction_ids(request):

    reaction_ids = User.objects.values_list('saved_reactions__id', flat=True).distinct()
    return JsonResponse(list(reaction_ids), safe=False)



@csrf_exempt
def is_reaction_in_user_saved_reactions(request):
    """Check if the given reaction_id is in the user's saved reactions."""
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('user_id')
            reaction_id = int(data.get('reaction_id'))  # Convert reaction_id to integer
            
            user = User.objects.get(pk=user_id)
            saved_reaction_ids = list(user.saved_reactions.all().values_list('id', flat=True))

            is_reaction_saved = reaction_id in saved_reaction_ids
            
            return JsonResponse({'is_reaction_saved': is_reaction_saved})
        except User.DoesNotExist:
            return JsonResponse({'error': 'User does not exist'}, status=404)
        except (KeyError, ValueError):
            return JsonResponse({'error': 'Invalid data'}, status=400)
    return JsonResponse({'error': 'Invalid request method'}, status=405)


def get_available_reactions(request):
    """Return the IDs of available reactions for the given user."""
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('user_id')

            user = User.objects.get(pk=user_id)
            saved_reaction_ids = list(user.saved_reactions.all().values_list('id', flat=True))
            available_reactions = Reaction.objects.exclude(id__in=saved_reaction_ids)
            available_reaction_ids = list(available_reactions.values_list('id', flat=True))
            last_index = saved_reaction_ids[-1] if saved_reaction_ids else None

             
            return JsonResponse({'available_reaction_ids': available_reaction_ids, 'last_index': last_index})
        except User.DoesNotExist:
            return JsonResponse({'error': 'User does not exist'}, status=404)
        except (KeyError, ValueError):
            return JsonResponse({'error': 'Invalid data'}, status=400)
    return JsonResponse({'error': 'Invalid request method'}, status=405)



def reaction_view(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            reaction_abbreviation = data.get('reactionAbbreviation')
            if reaction_abbreviation:
                result = get_from_rhea(reaction_abbreviation)
                return JsonResponse(result)
            else:
                return JsonResponse({'error': 'Missing reaction abbreviation'}, status=400)
        except json.JSONDecodeError:
            return JsonResponse({'error': 'Invalid JSON'}, status=400)
    else:
        return JsonResponse({'error': 'Invalid request'}, status=400)

@csrf_exempt
def get_rxn_template(request):
    """
    Django view to return substrates and products templates based on the reaction type.

    Input:
    - request: The Django HTTP request object.

    Output:
    - JsonResponse: A JsonResponse object containing the substrates and products lists.
    """
    if request.method == 'POST':
        data = json.loads(request.body)  # Parse JSON data from request body
        reaction_type = data.get('reaction_type', '')

        def create_entry(components):
            substrates = [component if component != 'empty' else 'empty' for component in components['substrates']]
            subs_sch = ['1' if component != 'empty' else '' for component in components['substrates']]
            subs_comps = ['-' if component != 'empty' else '' for component in components['substrates']]
            subs_types = components['subs_types']
            products = [component if component != 'empty' else 'empty' for component in components['products']]
            prod_sch = ['1' if component != 'empty' else '' for component in components['products']]
            prods_comps = ['-' if component != 'empty' else '' for component in components['products']]
            prods_types = components['prods_types']

            return {
                'substrates': substrates,
                'subs_sch': subs_sch,
                'subs_comps': subs_comps,
                'subs_types': subs_types,
                'products': products,
                'prod_sch': prod_sch,
                'prods_comps': prods_comps,
                'prods_types': prods_types
            }

        reaction_templates = {
            'hydrolysis': {
                'substrates': ['empty', 'h2o', 'h'],
                'subs_types': ['', '', ''],
                'products': ['empty', 'empty', 'h'],
                'prods_types': ['', '', '']
            },
            'O2NADPHOX': {
                'substrates': ['empty', 'o2', 'nadph', 'h'],
                'subs_types': ['', '', '', ''],
                'products': ['empty', 'h2o', 'nadp'],
                'prods_types': ['', '', '']
            },
            'SULT': {
                'substrates': ['empty', 'paps'],
                'subs_types': ['', ''],
                'products': ['empty', 'pap', 'h'],
                'prods_types': ['', '', '']
            },
            'UGT': {
                'substrates': ['empty', 'udpglcur'],
                'subs_types': ['', ''],
                'products': ['empty', 'udp', 'h'],
                'prods_types': ['', '', '']
            },
            'UGT glucose': {
                'substrates': ['empty', 'udpg'],
                'subs_types': ['', ''],
                'products': ['empty', 'udp', 'h'],
                'prods_types': ['', '', '']
            },
            'UGT carb glucur': {
                'substrates': ['empty', 'udpglcur', 'co2'],
                'subs_types': ['', '', ''],
                'products': ['empty', 'udp', 'h'],
                'prods_types': ['', '', '']
            },
            'CoA': {
                'substrates': ['empty', 'coa', 'atp'],
                'subs_types': ['', '', ''],
                'products': ['empty', 'amp', 'ppi'],
                'prods_types': ['', '', '']
            },
            'FAOXhd': {
                'substrates': ['empty', 'h2o'],
                'subs_types': ['', ''],
                'products': ['empty'],
                'prods_types': ['']
            },
            'FAOXnad': {
                'substrates': ['empty', 'nad'],
                'subs_types': ['', ''],
                'products': ['empty', 'nadh', 'h'],
                'prods_types': ['', '', '']
            },
            'FAOXcoa': {
                'substrates': ['empty', 'coa'],
                'subs_types': ['', ''],
                'products': ['empty', 'accoa'],
                'prods_types': ['', '']
            },
            'Nad ox': {
                'substrates': ['empty', 'nad'],
                'subs_types': ['', ''],
                'products': ['empty', 'nadh', 'h'],
                'prods_types': ['', '', '']
            },
            'NADH red': {
                'substrates': ['empty', 'nadh', 'h'],
                'subs_types': ['', '', ''],
                'products': ['empty', 'nad'],
                'prods_types': ['', '']
            },
            'cycl': {
                'substrates': ['empty', 'h'],
                'subs_types': ['', ''],
                'products': ['empty', 'h2o'],
                'prods_types': ['', '']
            },
            'NADPH red': {
               'substrates': ['empty', 'nadph', 'h'],
                'subs_types': ['', '', ''],
                'products': ['empty', 'nadp'],
                'prods_types': ['', '']
            },
            'AT':{
                'substrates': ['taur', 'empty'],
                'subs_types': ['', ''],
                'products': ['empty', 'coa','h'],
                'prods_types': ['', '', '']
            }         
        }

        components = reaction_templates.get(reaction_type, {
            'substrates': [],
            'subs_types': [],
            'products': [],
            'prods_types': []
        })

        result = create_entry(components)
        return JsonResponse(result)
    else:
        return JsonResponse({'error': 'Invalid request'}, status=400)

        
def search_reactions(request):
    user = request.user
    query = request.GET.get('q', '')
    if query:
        reactions = user.saved_reactions.filter(
            substrates__icontains=query) | user.saved_reactions.filter(
            products__icontains=query) | user.saved_reactions.filter(
            short_name__icontains=query) | user.saved_reactions.filter(
            direction__icontains=query)
    else:
        reactions = user.saved_reactions.all()
    reactions_data = [
        {
            'substrates': reaction.substrates,
            'products': reaction.products,
            'short_name': reaction.short_name,
            'direction': reaction.direction,
            # Add more fields if needed
        }
        for reaction in reactions
    ]
    return JsonResponse({'reactions': reactions_data})

def leader_board(request):
    return render(request, 'reactions/leader_board.html')

@require_GET
def get_user_reactions_and_vmh(request):
    """Return all existing users' full names, the number of reactions saved, reactions added to VMH, and reactions created."""
    try:
        users_data = []
        users = User.objects.all()

        for user in users:
            saved_reactions_count = user.saved_reactions.count()
            reactions_added_vmh_count = ReactionsAddedVMH.objects.filter(user=user).count()
            created_reactions_count = CreatedReaction.objects.filter(user=user).count()

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
    

def create_reaction(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        user_id = data.get('user_id')
        reaction_id = data.get('reaction_id')
        
        # Validate the user ID using the provided function
        user = validate_user_ID(user_id)
        
        if not user:
            return JsonResponse({'success': False, 'error': 'User does not exist'}, status=400)
        
        # Fetch the reaction based on reaction_id
        reaction = get_object_or_404(Reaction, id=reaction_id)
        
        created_reaction = CreatedReaction.objects.create(user=user, reaction=reaction)
        
        
        return JsonResponse({'success': True, 'created_reaction_id': created_reaction.id})
    else:
        return JsonResponse({'success': False, 'error': 'Invalid request method'}, status=400)
    
@csrf_exempt

def parse_formula_with_compartments(request):
    if request.method == 'POST':
        try:
            body = json.loads(request.body)
        except json.JSONDecodeError as e:
            return JsonResponse({'error': f'JSON decode error: {str(e)}'}, status=400)

        formulas = body.get('formulas', [])
        subs_comps = body.get('subs_comps', [])
        prods_comps = body.get('prods_comps', [])

        # Log received data for debugging

        detailed_formulas = []

        def extract_components(comp_str):
            """Extracts quantity and name from a component string."""
            match = re.match(r'(\d*)\s*(\w+)', comp_str.strip())
            if match:
                quantity = match.group(1) if match.group(1) else '1'
                name = match.group(2)
                return quantity, name
            return None, None

        for formula in formulas:
            # Normalize formula by ensuring spaces around the arrows
            formula = formula.replace('->', ' -> ').replace('<=>', ' <=> ')

            if ' -> ' in formula:
                direction = 'forward'
                substrates, products = formula.split(' -> ')
            elif ' <=> ' in formula:
                direction = 'reversible'
                substrates, products = formula.split(' <=> ')
            else:
                return JsonResponse({'error': f'Invalid formula format: {formula}'}, status=400)

            subs_list = [comp.strip() for comp in substrates.split(' + ')] if substrates else []
            prods_list = [comp.strip() for comp in products.split(' + ')] if products else []

            # Remove empty components
            subs_list = [comp for comp in subs_list if comp]
            prods_list = [comp for comp in prods_list if comp]



            # Check if the lengths of subs_comps and prods_comps match the lengths of subs_list and prods_list
            if len(subs_list) != len(subs_comps):
                return JsonResponse({'error': 'Number of substrate compartments does not match number of substrates'}, status=400)

            if len(prods_list) != len(prods_comps):
                return JsonResponse({'error': 'Number of product compartments does not match number of products'}, status=400)

            subs_detailed = []
            for j, comp in enumerate(subs_list):
                quantity, name = extract_components(comp)
                subs_detailed.append(f"{quantity} {name}[{subs_comps[j]}]")

            prods_detailed = []
            for j, comp in enumerate(prods_list):
                quantity, name = extract_components(comp)
                prods_detailed.append(f"{quantity} {name}[{prods_comps[j]}]")

            if direction == 'forward':
                detailed_formula = ' + '.join(subs_detailed) + ' -> ' + ' + '.join(prods_detailed)
            else:  # reversible
                detailed_formula = ' + '.join(subs_detailed) + ' <=> ' + ' + '.join(prods_detailed)

            detailed_formulas.append(detailed_formula)

        return JsonResponse({'detailed_formulas': detailed_formulas})

    return JsonResponse({'error': 'Invalid request method'}, status=400)

@csrf_exempt
def convert_to_smiles(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        substrates = [item for item in data if item['side'] == 'substrate']
        products = [item for item in data if item['side'] == 'product']
        
        substrates_mols = [substrate['mol'] for substrate in substrates]
        substrates_types = [substrate['type'] for substrate in substrates]
        
        products_mols = [product['mol'] for product in products]
        products_types = [product['type'] for product in products]
        
        substrates_smiles, substrates_errors = any_to_smiles(substrates_mols, substrates_types, request, side='substrates')
        products_smiles, products_errors = any_to_smiles(products_mols, products_types, request, side='products')
        
        response_data = {
            'substrates_smiles': substrates_smiles,
            'substrates_errors': substrates_errors,
            'products_smiles': products_smiles,
            'products_errors': products_errors
        }
        
        return JsonResponse(response_data)
    else:
        return JsonResponse({'error': 'Invalid request method'}, status=400)


@csrf_exempt  # Use csrf_exempt if CSRF token isn't being managed
def save_formula(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            formula = data.get('formula', '')
            reaction_id = data.get('reaction_id', '')

            # Ensure the formula is saved with quotes
            quoted_formula = f'"{formula}"'

            # Find the reaction by the provided reaction ID
            try:
                reaction = Reaction.objects.get(id=reaction_id)
            except Reaction.DoesNotExist:
                return JsonResponse({'status': 'failed', 'error': 'Reaction not found'}, status=404)

            # Update the formulas field with the provided formula string with quotes
            reaction.rxn_formula = quoted_formula
            reaction.save()

            return JsonResponse({'status': 'success'})
        except Exception as e:
            return JsonResponse({'status': 'failed', 'error': str(e)}, status=500)
    return JsonResponse({'status': 'failed', 'error': 'Invalid request method'}, status=400)
