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

def update_vmh_from_constructor(json_dir, matlab_session, update_existing=False, dry_run=False):
    """
    Execute the unified MATLAB function updateVMHFromConstructor to add/update 
    reactions and metabolites in VMH.
    
    This replaces the need to call both add_rxn_python and add_metab_python separately.
    The MATLAB function:
    - Checks which metabolites in the formulas are new (not in VMH)
    - Adds new metabolites to metabolites table + compartment associations
    - For each reaction:
        - If new: Adds to reactions, recon, and reconws_SMatrix tables
        - If exists and update_existing=True: Updates existing entries
        - If exists and update_existing=False: Skips (default behavior)
    
    Args:
        json_dir: Directory containing the 10 required JSON files:
            - reactionIds.json
            - reactionNames.json
            - reactionFormulas.json
            - reactionDirections.json
            - reactionSubsystems.json
            - reactionReferences.json
            - reactionExternalLinks.json
            - reactionGeneInfo.json
            - reactionComments.json
            - reactionConfidenceScores.json
        matlab_session: MatlabHTTPClient instance
        update_existing: Whether to update reactions that already exist in VMH (default False)
        dry_run: Preview changes without committing to DB (default False)
        
    Returns:
        Dict with:
            - 'status': 'success' or 'error'
            - 'addedMets': List of abbreviations of newly added metabolites
            - 'addedRxns': List of abbreviations of newly added reactions
            - 'updatedRxns': List of abbreviations of updated reactions
            - 'message': Error message if status is 'error'
    """
    # Build options struct for MATLAB
    options = {
        'updateExisting': update_existing,
        'dryRun': dry_run,
        'verbose': True,
        'model_id': 8564  # RECON4IMD
    }
    
    result = matlab_session.execute('updateVMHFromConstructor', json_dir, options)
    
    if result.get('status') == 'success':
        matlab_result = result.get('result', {})
        return {
            'status': 'success',
            'addedMets': matlab_result.get('addedMets', []),
            'addedRxns': matlab_result.get('addedRxns', []),
            'updatedRxns': matlab_result.get('updatedRxns', [])
        }
    else:
        return {
            'status': 'error',
            'message': result.get('message', 'Unknown error from MATLAB'),
            'addedMets': [],
            'addedRxns': [],
            'updatedRxns': []
        }

def parse_gene_info(info):
    """
    Parse a gene info string to extract just the clean GPR rule.
    
    The input may contain additional metadata like:
    "GPR: ABL1 AND AOC1; ORGAN(Adipocytes_), SUBCELLULAR([e], [c])"
    
    We want to extract just: "ABL1 and AOC1"
    """
    if not info:
        return ""
    
    # First, take only the part before the semicolon (removes ORGAN, SUBCELLULAR metadata)
    first_part = info.split(';')[0].strip()
    
    # Remove "GPR: " prefix if present
    if first_part.startswith("GPR: "):
        first_part = first_part[5:]  # Remove "GPR: " prefix
    
    # Normalize AND/OR to lowercase for MATLAB compatibility
    if " AND " in first_part:
        first_part = first_part.replace(" AND ", " and ")
    if " OR " in first_part:
        first_part = first_part.replace(" OR ", " or ")
    
    return first_part

def merge_gene_infos(gene_infos):
    """
    if more than 1 gene info for a reaction, merge them into a single GPR string (add "or" between them)
    """
    infos = [g['info'] for g in gene_infos]
    gprs = [parse_gene_info(info) for info in infos]
    if len(gprs) == 1:
        return gprs[0]
    if len(gprs) == 0:
        return ""
    merged_gpr = " or ".join(gprs)
    return merged_gpr

def prepare_vmh_update_json_files(
        reaction_identifiers,
        reaction_names,
        reaction_formulas,
        reaction_directions,
        reaction_subsystems,
        reaction_references,
        reaction_external_links,
        reaction_gene_info,
        reaction_comments,
        reaction_confidence_scores,
        output_dir=None):
    """
    Prepares all 10 JSON files required by updateVMHFromConstructor in a dedicated directory.
    
    Args:
        reaction_identifiers: List of reaction abbreviations
        reaction_names: List of reaction names/descriptions
        reaction_formulas: List of reaction formulas (e.g., "glc_D[e] => glc_D[c]")
        reaction_directions: List of directions ("forward", "bidirectional", "reverse")
        reaction_subsystems: List of subsystem assignments
        reaction_references: List of reference structs [{info, ref_type}, ...]
        reaction_external_links: List of external link structs [{ext_link_type, info}, ...]
        reaction_gene_info: List of GPR rule structs [{info}, ...]
        reaction_comments: List of comment structs [{info}, ...]
        reaction_confidence_scores: List of confidence scores
        output_dir: Optional directory path. If None, creates a temp directory.
        
    Returns:
        str: Path to the directory containing all JSON files
    """
    import tempfile
    import uuid
    
    if output_dir is None:
        # Create a unique temp directory
        output_dir = os.path.join(tempfile.gettempdir(), f'vmh_update_{uuid.uuid4().hex}')
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Process gene info to merge GPRs
    processed_gene_info = []
    for gene_info in reaction_gene_info:
        if gene_info:
            merged_gpr = merge_gene_infos(gene_info) 
            processed_gene_info.append({'info': merged_gpr} if merged_gpr else {})
        else:
            processed_gene_info.append({})

    # Map of filename to data
    files_data = {
        'reactionIds.json': reaction_identifiers,
        'reactionNames.json': reaction_names,
        'reactionFormulas.json': reaction_formulas,
        'reactionDirections.json': reaction_directions,
        'reactionSubsystems.json': reaction_subsystems,
        'reactionReferences.json': reaction_references,
        'reactionExternalLinks.json': reaction_external_links,
        'reactionGeneInfo.json': processed_gene_info,
        'reactionComments.json': reaction_comments,
        'reactionConfidenceScores.json': reaction_confidence_scores
    }
    
    for filename, data in files_data.items():
        filepath = os.path.join(output_dir, filename)
        save_json(data, filepath)
    
    return output_dir


def cleanup_vmh_update_json_files(json_dir):
    """
    Remove the temporary directory and all JSON files created for VMH update.
    
    Args:
        json_dir: Path to the directory containing the JSON files
    """
    import shutil
    if os.path.exists(json_dir):
        shutil.rmtree(json_dir)


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
        # Use abbreviation or pk to identify reactions with missing names
        missing_reactions = [
            reaction.get('abbreviation') or f"Reaction #{reaction.get('pk', 'unknown')}"
            for reaction, missing in zip(reactions, missing_names) if missing
        ]
        return JsonResponse({'status': 'error',
                             'message': f'Missing description for: {", ".join(missing_reactions)}'})
    
    missing_abbrs = [reaction['abbreviation'] == '' for reaction in reactions]
    if True in missing_abbrs:
        missing_reactions = [
            reaction['description'] or f"Reaction #{reaction.get('pk', 'unknown')}"
            for reaction, missing in zip(reactions, missing_abbrs) if missing
        ]
        return JsonResponse({'status': 'error',
                             'message': f'Missing abbreviation for: {", ".join(missing_reactions)}'})
    
    missing_conf_scores = [reaction['confidence_score'] == '" "' for reaction in reactions]
    if True in missing_conf_scores:
        missing_reactions = [
            reaction.get('abbreviation') or reaction.get('description') or f"Reaction #{reaction.get('pk', 'unknown')}"
            for reaction, missing in zip(reactions, missing_conf_scores) if missing
        ]
        return JsonResponse({'status': 'error',
                             'message': f'Missing confidence score for: {", ".join(missing_reactions)}'})
    
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
