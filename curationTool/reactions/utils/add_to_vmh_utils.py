import json
import requests
import random
import os

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from django.http import JsonResponse

from reactions.utils.to_mol import any_to_mol
from reactions.utils.search_vmh import check_reaction_vmh
# Function to gather additional reaction details

from django.conf import settings


def gather_reaction_details(reaction_objs):
    """
    Gathers additional details for reactions, such as direction, references, external links, gene info, comments, and confidence scores.
    """
    reaction_directions = [reaction.direction for reaction in reaction_objs]
    reaction_subsystems = [reaction.subsystem for reaction in reaction_objs]
    reaction_references = [
        reaction.references if reaction.references else [] for reaction in reaction_objs]
    reaction_external_links = [
        reaction.ext_links if reaction.ext_links else [] for reaction in reaction_objs]
    reaction_gene_info = [
        reaction.gene_info if reaction.gene_info else [] for reaction in reaction_objs]
    reaction_comments = [reaction.comments if reaction.comments else []
                         for reaction in reaction_objs]
    reaction_confidence_scores = [
        reaction.confidence_score for reaction in reaction_objs]
    return reaction_directions, reaction_subsystems, reaction_references, reaction_external_links, reaction_gene_info, reaction_comments, reaction_confidence_scores


def save_json(data, filepath):
    """
    Saves the provided data into a JSON file at the given filepath.
    """
    with open(filepath, 'w') as f:
        json.dump(data, f)


def rxn_prepare_json_paths_and_variables(
        reaction_identifiers,
        reaction_names,
        reaction_formulas,
        reaction_directions,
        reaction_subsystems,
        reaction_references,
        reaction_external_links,
        reaction_gene_info,
        reaction_comments,
        reaction_confidence_scores):
    """
    Prepares JSON paths and variables for MATLAB execution, saving them to temporary files.
    """
    rand_float = random.uniform(0, 10000000)
    json_paths = [
        f'reactionIds.json{rand_float}',
        f'reactionNames.json{rand_float}',
        f'reactionFormulas.json{rand_float}',
        f'reactionDirections.json{rand_float}',
        f'reactionSubsystems.json{rand_float}',
        f'reactionReferences.json{rand_float}',
        f'reactionExternalLinks.json{rand_float}',
        f'reactionGeneInfo.json{rand_float}',
        f'reactionComments.json{rand_float}',
        f'reactionConfidenceScores.json{rand_float}']
    variables = [
        reaction_identifiers,
        reaction_names,
        reaction_formulas,
        reaction_directions,
        reaction_subsystems,
        reaction_references,
        reaction_external_links,
        reaction_gene_info,
        reaction_comments,
        reaction_confidence_scores]
    for idx, (path, variable) in enumerate(zip(json_paths, variables)):
        path = os.path.join(os.getcwd(), path)
        json_paths[idx] = path
        save_json(variable, path)
    return json_paths


def met_prepare_json_paths_and_variables(
        met_abbrs,
        met_names,
        met_formulas,
        met_charges,
        met_inchikeys,
        met_smiles,
        met_external_links,
        met_weights):
    """
    Prepares JSON paths and variables for MATLAB execution, saving them to temporary files.
    """
    rand_float = random.uniform(0, 10000000)
    json_paths = [
        f'metAbbrs.json{rand_float}',
        f'metNames.json{rand_float}',
        f'metFormulas.json{rand_float}',
        f'metInchikeys.json{rand_float}',
        f'metSmiles.json{rand_float}',
        f'metCharges.json{rand_float}',
        f'metExternalLinks.json{rand_float}',
        f'metWeights.json{rand_float}']
    variables = [
        met_abbrs,
        met_names,
        met_formulas,
        met_inchikeys,
        met_smiles,
        met_charges,
        met_external_links,
        met_weights]
    for idx, (path, variable) in enumerate(zip(json_paths, variables)):
        path = os.path.join(os.getcwd(), path)
        json_paths[idx] = path
        save_json(variable, path)
    return json_paths


def add_reaction_matlab(json_paths, matlab_session):
    """
    Executes MATLAB operations to add reactions to VMH, using the provided JSON paths.
    """
    result = matlab_session.execute('add_rxn_python', *json_paths)
    result['rxn_ids'] = result['result'] if result['status'] == 'success' else []
    return result


def add_metabolites_matlab(json_paths, matlab_session):
    """
    Executes MATLAB operations to add metabolites to VMH, using the provided JSON paths.
    """
    result = matlab_session.execute('add_metab_python', *json_paths)
    result['met_ids'] = result['result'] if result['status'] == 'success' else []
    return result


def smiles_to_inchikeys(smiles_list):
    inchi_list = []
    for smiles in smiles_list:
        # Convert the SMILES string to a molecule object
        mol = Chem.MolFromSmiles(smiles)
        if mol:  # Ensure the molecule was created successfully
            # Convert the molecule to an InChI string
            inchi = Chem.MolToInchi(mol)
            inchi_key = Chem.inchi.InchiToInchiKey(inchi)
            inchi_list.append(inchi_key)
        else:
            inchi_list.append(None)  # Append None if the conversion fails
        inchi_list = ['' if x is None else x for x in inchi_list]
    return inchi_list


def smiles_to_charged_formula(smiles_list):
    charged_formulas = []
    charges = []

    for smiles in smiles_list:
        # Convert the SMILES string to a molecule object
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Adding hydrogen to the molecule to ensure the formula includes
            # hydrogens
            mol = Chem.AddHs(mol)
            # Calculate the molecular formula
            formula = CalcMolFormula(mol)
            # Calculate the net charge
            charge = Chem.GetFormalCharge(mol)
            charged_formulas.append(formula)
            charges.append(charge)
        else:
            # Handle the case where SMILES conversion fails
            charged_formulas.append(None)
            charges.append(None)
    charged_formulas = ['' if x is None else x for x in charged_formulas]
    charges = ['' if x is None else x for x in charges]
    return charged_formulas, charges


def get_nonfound_metabolites(
        reaction_objs,
        subs_abbr,
        prods_abbr,
        search_func):
    react_metab_notfound = {}
    for react_idx, reaction in enumerate(reaction_objs):
        subs, subs_types, subs_names = json.loads(
            reaction.substrates), json.loads(
            reaction.substrates_types), json.loads(
            reaction.substrates_names)

        prods, prods_types, prods_names = json.loads(
            reaction.products), json.loads(
            reaction.products_types), json.loads(
            reaction.products_names)

        subs_founds, subs_miriams = search_func(
            subs, subs_types, None, side='substrates', nofile=True)
        prods_founds, prods_miriams = search_func(
            prods, prods_types, None, side='products', nofile=True)
        subs_not_found = [(sub,
                           subs_types[idx],
                           subs_abbr[react_idx][idx],
                           subs_names[idx]) for idx,
                          (sub,
                           found) in enumerate(zip(subs,
                                                   subs_founds)) if not found]
        prods_not_found = [(prod,
                            prods_types[idx],
                            prods_abbr[react_idx][idx],
                            prods_names[idx]) for idx,
                           (prod,
                            found) in enumerate(zip(prods,
                                                    prods_founds)) if not found]
        react_metab_notfound[reaction.id] = {
            'subs': subs_not_found, 'prods': prods_not_found}

    # Add the ones not found to the list of metabolites to add to the VMH

    abbrs_all_subs_not_found, abbrs_all_prods_not_found, mols_all_subs_not_found, mols_all_prods_not_found, types_all_subs_not_found, types_all_prods_not_found, names_all_subs_not_found, names_all_prods_not_found = [], [], [], [], [], [], [], []
    # Iterate through each reaction in react_metab_notfound
    for _, metabolites_info in react_metab_notfound.items():
        # Iterate through the subs_not_found for the current reaction
        for sub_info in metabolites_info['subs']:
            mols_all_subs_not_found.append(sub_info[0])
            types_all_subs_not_found.append(sub_info[1])
            abbrs_all_subs_not_found.append(sub_info[2])
            names_all_subs_not_found.append(sub_info[3])
        # Iterate through the prods_not_found for the current reaction
        for prod_info in metabolites_info['prods']:
            mols_all_prods_not_found.append(prod_info[0])
            types_all_prods_not_found.append(prod_info[1])
            abbrs_all_prods_not_found.append(prod_info[2])
            names_all_prods_not_found.append(prod_info[3])

    combined_abbrs = abbrs_all_subs_not_found + abbrs_all_prods_not_found
    combined_mols = mols_all_subs_not_found + mols_all_prods_not_found
    combined_types = types_all_subs_not_found + types_all_prods_not_found
    combined_names = names_all_subs_not_found + names_all_prods_not_found
    # Step 2: Identify unique abbreviations while maintaining order
    unique_abbrs = []
    indices_unique_abbrs = []

    for i, abbr in enumerate(combined_abbrs):
        if abbr not in unique_abbrs:
            unique_abbrs.append(abbr)
            indices_unique_abbrs.append(i)

    # Step 3: Map unique abbreviations to corresponding molecules and types
    unique_mols = [combined_mols[i] for i in indices_unique_abbrs]
    unique_types = [combined_types[i] for i in indices_unique_abbrs]
    unique_names = [combined_names[i] for i in indices_unique_abbrs]
    return unique_abbrs, unique_mols, unique_types, unique_names
# dO SAME FOR VMH DB


def check_reactions_vmh(reaction_objs):
    in_vmh = []
    for reaction in reaction_objs:
        substrates = json.loads(reaction.substrates)
        products = json.loads(reaction.products)
        substrates_types = json.loads(reaction.substrates_types)
        products_types = json.loads(reaction.products_types)
        subs_sch = json.loads(reaction.subs_sch)
        prod_sch = json.loads(reaction.prods_sch)
        direction = reaction.direction
        subsystem = reaction.subsystem
        subs_comps = json.loads(reaction.subs_comps)
        prods_comps = json.loads(reaction.prods_comps)
        subs_mols, subs_errors = any_to_mol(
            substrates, substrates_types, request=None, side='substrates')
        prod_mols, prod_errors = any_to_mol(
            products, products_types, request=None, side='products')
        all_errors = subs_errors + prod_errors
        if any(elem is not None for elem in all_errors):
            in_vmh.append(False)
            continue
        # Check if the reaction is found in VMH
        vmh_found = check_reaction_vmh(
            substrates,
            products,
            subs_sch,
            prod_sch,
            substrates_types,
            products_types,
            subs_mols,
            prod_mols,
            direction,
            subsystem,
            subs_comps,
            prods_comps)
        if vmh_found['found'] and not vmh_found['similar']:
            in_vmh.append(True)
        else:
            in_vmh.append(False)
    return in_vmh


def check_names_abbrs_vmh(names_abbr_list):
    names_vmh = {}
    abbr_vmh = {}
    BASE_URL = settings.OLD_VMH_BASE_URL
    for name, abbr in names_abbr_list:
        endpoint = f"{BASE_URL}_api/reactions/?abbreviation={abbr}"
        response = requests.get(endpoint, verify=False)
        found_abbr = False
        if response.json().get('count', 0) > 0:
            for result in response.json().get('results', []):
                if result['abbreviation'].lower() == abbr.lower():
                    found_abbr = True
                    break
        abbr_vmh[abbr] = found_abbr
        endpoint = f"{BASE_URL}_api/reactions/?description={name}"
        response = requests.get(endpoint, verify=False)
        found_name = False
        if response.json().get('count', 0) > 0:
            for result in response.json().get('results', []):
                if result['description'].lower() == name.lower():
                    found_name = True
                    break
        names_vmh[name] = found_name
    return names_vmh, abbr_vmh


def make_request_names_abbrs(name, abbr):
    BASE_URL = settings.OLD_VMH_BASE_URL
    endpoint = f"{BASE_URL}_api/metabolites/?abbreviation={abbr}"
    response = requests.get(endpoint, verify=False)
    abbr_found = False
    name_found = False
    if response.json().get('count', 0) > 0:
        for result in response.json().get('results', []):
            if result['abbreviation'].lower() == abbr.lower():
                abbr_found = True
                break

    endpoint = f"{BASE_URL}_api/metabolites/?fullName={name}"
    response = requests.get(endpoint, verify=False)
    if response.json().get('count', 0) > 0:
        for result in response.json().get('results', []):
            if result['fullName'].lower() == name.lower():
                name_found = True
                break

    return name_found, abbr_found


def check_met_names_abbrs_vmh(
        subs_info,
        prods_info,
        subs_founds,
        prods_founds):
    def update_vmh_info(items_info, names_vmh, abbr_vmh, found_info):
        for idx, item in enumerate(items_info):
            if found_info[idx]:
                continue
            name_found, abbr_found = make_request_names_abbrs(
                item['name'], item['abbreviation'])
            names_vmh[item['name']] = name_found
            abbr_vmh[item['abbreviation']] = abbr_found
    n_reactions = len(subs_info)
    subs_names_vmh = {}
    subs_abbr_vmh = {}
    prods_names_vmh = {}
    prods_abbr_vmh = {}
    for i in range(n_reactions):
        subs_found = subs_founds[i]
        prods_found = prods_founds[i]
        update_vmh_info(
            subs_info[i],
            subs_names_vmh,
            subs_abbr_vmh,
            subs_found)
        update_vmh_info(
            prods_info[i],
            prods_names_vmh,
            prods_abbr_vmh,
            prods_found)
    return subs_names_vmh, subs_abbr_vmh, prods_names_vmh, prods_abbr_vmh

def validate_reaction_fields(reactions):
    """
    Validate that each reaction has non-empty description, abbreviation, and confidence_score.
    Also check for duplicates in the list.
    """
    missing_names = [reaction['description'] == '' for reaction in reactions]
    if True in missing_names:
        return JsonResponse({'status': 'error',
                             'message': 'Please enter a description for reaction'})
    
    missing_abbrs = [reaction['abbreviation'] == '' for reaction in reactions]
    if True in missing_abbrs:
        return JsonResponse({'status': 'error',
                             'message': 'Please enter an abbreviation'})
    
    missing_conf_scores = [reaction['confidence_score'] == '" "' for reaction in reactions]
    if True in missing_conf_scores:
        return JsonResponse({'status': 'error',
                             'message': 'Please enter a confidence score for reaction'})
    
    names_list = [reaction['description'] for reaction in reactions]
    for name in names_list:
        if names_list.count(name) > 1:
            return JsonResponse({'status': 'error',
                                 'message': f'Reaction with name `{name}` is repeated in the list.'})
    
    abbr_list = [reaction['abbreviation'] for reaction in reactions]
    for abbr in abbr_list:
        if abbr_list.count(abbr) > 1:
            return JsonResponse({'status': 'error',
                                 'message': f'Reaction with abbreviation `{abbr}` is repeated in the list.'})
    return None


def validate_reaction_existence(reactions):
    """
    Validate that the reaction names and abbreviations do not already exist in VMH.
    """
    name_in_vmh, abbr_in_vmh = check_names_abbrs_vmh(
        [(reaction['description'], reaction['abbreviation']) for reaction in reactions]
    )
    if True in list(name_in_vmh.values()):
        name_in_vmh_reactions = [
            reaction for reaction, in_vmh in zip(reactions, name_in_vmh.values()) if in_vmh
        ]
        reaction_names = ", ".join([reaction["description"] for reaction in name_in_vmh_reactions])
        return JsonResponse({'status': 'error',
                             'message': f'The following reaction descriptions are already in VMH: {reaction_names}'})
    
    if True in list(abbr_in_vmh.values()):
        abbr_in_vmh_reactions = [
            reaction for reaction, in_vmh in zip(reactions, abbr_in_vmh.values()) if in_vmh
        ]
        reaction_abbrs = ", ".join([reaction["abbreviation"] for reaction in abbr_in_vmh_reactions])
        return JsonResponse({'status': 'error',
                             'message': f'The following reaction abbreviations are already in VMH: {reaction_abbrs}'})
    return None


def validate_metabolite_existence(reactions_new_subsInfo, reactions_new_prodsInfo, reactions_subs_found, reactions_prods_found):
    """
    Validate that substrate and product names and abbreviations do not already exist in VMH.
    """
    subs_names_vmh, subs_abbr_vmh, prods_names_vmh, prods_abbr_vmh = check_met_names_abbrs_vmh(
        reactions_new_subsInfo, reactions_new_prodsInfo, reactions_subs_found, reactions_prods_found
    )
    if True in list(subs_names_vmh.values()):
        subs_names_in_vmh = [sub for sub in subs_names_vmh.keys() if subs_names_vmh[sub]]
        return JsonResponse({'status': 'error',
                             'message': f'The following substrates have metabolite names that are already in VMH: {", ".join(subs_names_in_vmh)}'})
    
    if True in list(subs_abbr_vmh.values()):
        subs_abbrs_in_vmh = [sub for sub in subs_abbr_vmh.keys() if subs_abbr_vmh[sub]]
        return JsonResponse({'status': 'error',
                             'message': f'The following substrates have metabolite abbreviations that are already in VMH: {", ".join(subs_abbrs_in_vmh)}'})
    
    if True in list(prods_names_vmh.values()):
        prods_names_in_vmh = [prod for prod in prods_names_vmh.keys() if prods_names_vmh[prod]]
        return JsonResponse({'status': 'error',
                             'message': f'The following products have metabolite names that are already in VMH: {", ".join(prods_names_in_vmh)}'})
    
    if True in list(prods_abbr_vmh.values()):
        prods_abbrs_in_vmh = [prod for prod in prods_abbr_vmh.keys() if prods_abbr_vmh[prod]]
        return JsonResponse({'status': 'error',
                             'message': f'The following products have metabolite abbreviations that are already in VMH: {", ".join(prods_abbrs_in_vmh)}'})
    return None


def validate_reaction_objects(reaction_objs, not_enough_info, no_comments, not_balanced):
    """
    Validate that each updated reaction object has at least one reference/external link,
    at least one comment, and is balanced.
    """
    if True in not_enough_info:
        not_enough_info_reactions = [
            reaction for reaction, not_enough in zip(reaction_objs, not_enough_info) if not_enough
        ]
        return JsonResponse({
            'status': 'error',
            'message': f'The following reactions do not have at least one reference or external link: {", ".join([reaction.short_name for reaction in not_enough_info_reactions])}'
        })
    
    if True in no_comments:
        no_comments_reactions = [
            reaction for reaction, no_comment in zip(reaction_objs, no_comments) if no_comment
        ]
        return JsonResponse({
            'status': 'error',
            'message': f'The following reactions do not have at least one comment: {", ".join([reaction.short_name for reaction in no_comments_reactions])}'
        })
    
    if True in not_balanced:
        not_balanced_reactions = [
            reaction for reaction, not_bal in zip(reaction_objs, not_balanced) if not_bal
        ]
        return JsonResponse({
            'status': 'error',
            'message': f'The following reactions are not balanced: {", ".join([reaction.short_name for reaction in not_balanced_reactions])}'
        })
    return None
