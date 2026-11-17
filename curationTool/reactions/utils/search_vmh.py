# This includes a collection of functions for handling chemical compounds and reactions.
# It provides functionality to convert various molecular identifiers to SMILES strings, fetch information from VMH,
# and process chemical formulas. The functions are designed to integrate
# with a Django web application, handling both HTTP requests and file
# operations.
#
# All VMH access goes through reactions.utils.vmh_api (the VMH public API at
# https://www.vmh.life/api/public).  Functions here keep their historical
# "never raise" behaviour unless called with strict=True, in which case
# vmh_api.VMHUnavailableError propagates when VMH cannot be asked.

import json
import logging
import os
import re
from urllib.parse import unquote

import requests
from django.conf import settings
from django.core.files.temp import NamedTemporaryFile
from django.http import JsonResponse
from rdkit import Chem

from reactions.utils.to_mol import any_to_mol
from reactions.utils.to_smiles import any_to_smiles, smiles_with_explicit_hydrogens
from reactions.utils.vmh_api import (
    VMHUnavailableError,
    find_metabolite_by_full_name,
    lookup_metabolite,
    lookup_reaction,
    match_metabolite_structure,
    search_reactions_by_participants,
    vmh_metabolite_url,
    vmh_reaction_url,
)

logger = logging.getLogger(__name__)


def is_name_in_vmh(name):
    """
    Checks if a molecule name is found in the VMH database.

    Input:
    - name (str): The molecule name.

    Output:
    - (bool): True if found, False otherwise.
    """
    return bool(find_metabolite_by_full_name(name))


def any_to_vmh(mols, types, smiles):
    """
    Converts a list of molecules to their VMH database abbreviations.

    Input:
    - mols (list): List of molecule identifiers.
    - types (list): List of types corresponding to each molecule identifier.
    - smiles (list): List of SMILES strings for the molecules.

    Output:
    - (list): List of VMH database abbreviations or 'error' for unsuccessful conversions.

    Raises VMHUnavailableError when VMH cannot be asked.
    """
    mols_list = []
    for idx, mol in enumerate(mols):
        if types[idx] == 'VMH':
            mols_list.append(mol)
        else:
            try:
                smiles_explicit = smiles_with_explicit_hydrogens(str(smiles[idx]))
                m = Chem.MolFromSmiles(smiles_explicit, sanitize=False)
                if m is None:
                    mols_list.append('error')
                    continue
                inchi = Chem.MolToInchi(m)
                inchi_key = Chem.MolToInchiKey(m)
            except Exception:
                mols_list.append('error')
                continue

            match = match_metabolite_structure(
                inchi_key=inchi_key,
                inchi_string=inchi,
                smiles=smiles_explicit,
            )

            abbr = (match or {}).get('abbreviation', '')
            mols_list.append(abbr if abbr else 'error')
    return mols_list


def _no_vmh_match(error=None):
    result = {"found": False, "similar": False, "url": None, "formula": None}
    if error:
        result["error"] = error
    return result


def check_reaction_vmh(
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
        prods_comps):
    """
    Checks if a chemical reaction is found in the VMH database based on substrates and products.

    Input:
    - substrates, products (list): Lists of substrate and product molecules.
    - subs_sch, prod_sch (list): Lists of stoichiometries for substrates and products.
    - substrates_types, products_types (list): Lists of types for substrates and products.
    - subs_mols, prod_mols (list): RDKit molecules for substrates and products.

    Output:
    - (dict): {'found', 'similar', 'url', 'formula'}.  'found' means VMH has a
      reaction with the same metabolites and stoichiometry on each side;
      'similar' means it differs in compartments, direction or subsystem.  When
      VMH could not be asked the result is "not found" plus an 'error' message.
    """
    try:
        return _check_reaction_vmh(
            substrates, products, subs_sch, prod_sch, substrates_types,
            products_types, subs_mols, prod_mols, direction, subsystem,
            subs_comps, prods_comps)
    except VMHUnavailableError as exc:
        logger.warning("VMH reaction check skipped: %s", exc.reason)
        return _no_vmh_match(exc.user_message("check the reaction against VMH"))


def _check_reaction_vmh(substrates, products, subs_sch, prod_sch, substrates_types,
                        products_types, subs_mols, prod_mols, direction, subsystem,
                        subs_comps, prods_comps):
    subs_smiles = any_to_smiles(
        subs_mols, ['MDL Mol file' for _ in range(len(subs_mols))],
        request=None, side='substrates')
    prod_smiles = any_to_smiles(
        prod_mols, ['MDL Mol file' for _ in range(len(prod_mols))],
        request=None, side='products')
    if None in subs_smiles[0] or None in prod_smiles[0]:
        return _no_vmh_match()
    substrates = any_to_vmh(substrates, substrates_types, subs_smiles[0])
    products = any_to_vmh(products, products_types, prod_smiles[0])
    if not substrates and not products:
        return _no_vmh_match()
    if 'error' in substrates or 'error' in products:
        # A metabolite that is not in VMH cannot be part of a VMH reaction.
        return _no_vmh_match()

    subs_sch = [float(s) for s in subs_sch]
    prod_sch = [float(p) for p in prod_sch]
    substrate_to_stoich = dict(zip(substrates, subs_sch))
    product_to_stoich = dict(zip(products, prod_sch))
    substrate_to_comps = dict(zip(substrates, subs_comps))
    product_to_comps = dict(zip(products, prods_comps))

    similar_match = None
    # One query for reactions with exactly these metabolites on each side; the
    # comparison below checks stoichiometry, compartments, direction, subsystem.
    for candidate in search_reactions_by_participants(substrates, products):
        formula = candidate.get("formula", "")
        try:
            (this_substrates, this_products, this_subs_stoich, this_prod_stoich,
             this_subs_comps, this_prods_comps, this_direction) = decode_formula(formula)
        except ValueError:
            continue

        left_mole_to_stoich = dict(zip(this_substrates, this_subs_stoich))
        right_mole_to_stoich = dict(zip(this_products, this_prod_stoich))
        left_mole_to_comps = dict(zip(this_substrates, this_subs_comps))
        right_mole_to_comps = dict(zip(this_products, this_prods_comps))

        substrates_match = (
            all(m in left_mole_to_stoich and substrate_to_stoich[m] == left_mole_to_stoich[m]
                for m in substrates)
            and all(m in substrate_to_stoich and left_mole_to_stoich[m] == substrate_to_stoich[m]
                    for m in this_substrates))
        products_match = (
            all(m in right_mole_to_stoich and product_to_stoich[m] == right_mole_to_stoich[m]
                for m in products)
            and all(m in product_to_stoich and right_mole_to_stoich[m] == product_to_stoich[m]
                    for m in this_products))
        if not (substrates_match and products_match):
            continue

        substrates_comps_match = (
            all(m in left_mole_to_comps and substrate_to_comps[m] == left_mole_to_comps[m]
                for m in substrates)
            and all(m in substrate_to_comps and left_mole_to_comps[m] == substrate_to_comps[m]
                    for m in this_substrates))
        products_comps_match = (
            all(m in right_mole_to_comps and product_to_comps[m] == right_mole_to_comps[m]
                for m in products)
            and all(m in product_to_comps and right_mole_to_comps[m] == product_to_comps[m]
                    for m in this_products))
        comps_match = substrates_comps_match and products_comps_match

        abbr = candidate.get("abbreviation") or ""
        match = {
            "found": True,
            "similar": not (candidate.get("subsystem", "") == subsystem
                            and this_direction == direction and comps_match),
            "url": vmh_reaction_url(abbr) if abbr else "",
            "formula": formula,
        }
        if not match["similar"]:
            return match
        if similar_match is None:
            similar_match = match
    return similar_match or _no_vmh_match()


def match_mol_in_vmh(mol):
    """
    Finds the VMH metabolite with the same structure as an RDKit molecule.

    Tries the explicit-hydrogen form first, then the molecule as given.
    Returns the VMH record or None; raises VMHUnavailableError when VMH
    cannot be asked.
    """
    smiles = Chem.MolToSmiles(mol)
    smiles = smiles_with_explicit_hydrogens(smiles)
    m = Chem.MolFromSmiles(str(smiles), sanitize=False)
    try:
        Chem.SanitizeMol(m)
    except BaseException:
        pass
    try:
        Chem.AssignStereochemistry(m)
    except BaseException:
        pass
    inchi = Chem.MolToInchi(m)
    if inchi.strip() == '':
        return None
    inchi_key = Chem.MolToInchiKey(m) if m else ''
    match = match_metabolite_structure(
        inchi_key=inchi_key,
        inchi_string=inchi,
        smiles=smiles,
    )
    if match:
        return match

    try:
        raw_inchi = Chem.MolToInchi(mol)
        raw_inchi_key = Chem.MolToInchiKey(mol)
    except Exception:
        raw_inchi = ''
        raw_inchi_key = ''
    if raw_inchi and raw_inchi != inchi:
        try:
            raw_smiles = smiles_with_explicit_hydrogens(Chem.MolToSmiles(mol))
        except Exception:
            raw_smiles = ""
        return match_metabolite_structure(
            inchi_key=raw_inchi_key,
            inchi_string=raw_inchi,
            smiles=raw_smiles,
        )
    return None


def _search_result(match, return_abbr=False, return_name=False):
    abbr = (match or {}).get('abbreviation', '')
    name = (match or {}).get('fullName', '')
    found = bool(abbr)
    miriam = vmh_metabolite_url(abbr) if abbr else ''
    if return_abbr and return_name:
        return found, miriam, abbr, name
    if return_abbr:
        return found, abbr
    if return_name:
        return found, name
    return found, miriam


def search_vmh(mol, return_abbr=False, return_name=False, strict=False):
    """
    Searches the VMH database for a molecule.

    Input:
    - mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object.
    - strict (bool): raise VMHUnavailableError instead of reporting "not found"
      when VMH cannot be asked.

    Output:
    - (tuple): (found, miriam), or (found, abbr) / (found, name) /
      (found, miriam, abbr, name) depending on return_abbr / return_name.
    """
    try:
        match = match_mol_in_vmh(mol)
    except VMHUnavailableError as exc:
        if strict:
            raise
        logger.warning("VMH structure search skipped: %s", exc.reason)
        match = None
    return _search_result(match, return_abbr=return_abbr, return_name=return_name)


def get_vmh_miriam(abbr, strict=False):
    """
    Retrieves the VMH link for a molecule abbreviation.

    Input:
    - abbr (str): The molecule abbreviation.
    - strict (bool): raise VMHUnavailableError when VMH cannot be asked.

    Output:
    - (str): The VMH link if the metabolite is found, else empty string.
    """
    try:
        match = lookup_metabolite(abbr)
    except VMHUnavailableError as exc:
        if strict:
            raise
        logger.warning("VMH lookup of %r skipped: %s", abbr, exc.reason)
        match = None
    if match and not match.get('is_former'):
        resolved_abbr = match.get("abbreviation") or abbr
        return vmh_metabolite_url(resolved_abbr)

    # Fallback: local mol file
    mol_path = os.path.join(settings.MOL_FILE_PATH, f"{abbr}.mol")
    if os.path.exists(mol_path):
        try:
            mol = Chem.MolFromMolFile(mol_path, sanitize=False, removeHs=False)
            found, miriam = search_vmh(mol, strict=strict)
            return miriam if found else ''
        except VMHUnavailableError:
            raise
        except Exception:
            return ''
    return ''


def search_metabolites_vmh(
        mols,
        types,
        request,
        side='substrates',
        nofile=False,
        return_abbr=False,
        strict=False):
    """
    Searches for metabolites in the VMH database based on a list of molecule identifiers.

    Input:
    - mols (list): List of molecule identifiers.
    - type (list): Corresponding types for each molecule identifier.
    - request: The Django request object (for file handling).
    - side (str): Indicating 'substrates' or 'products'.
    - strict (bool): raise VMHUnavailableError instead of reporting "not found"
      when VMH cannot be asked.

    Output:
    - (tuple): Tuple of two lists - one indicating if each molecule was found, and another with MIRIAM IDs.
    """
    def vmh_search(mol_obj):
        return search_vmh(mol_obj, return_abbr=return_abbr, strict=strict)

    founds = []
    miriams = []
    file_idx = 0  # Initialize file index
    for idx, this_type in enumerate(types):
        if this_type == 'VMH':
            miriam = get_vmh_miriam(mols[idx], strict=strict)
            found = bool(miriam)
        elif this_type == 'SwissLipids':
            base_url = 'https://www.swisslipids.org/api/index.php/entity/'
            endpoint = f"{base_url}{mols[idx]}"
            response = requests.get(endpoint, verify=False)
            if response.status_code != 200:
                found = False
                miriam = ''
            else:
                data = response.json()
                structures = data.get('structures', [])
                swiss_smiles = structures.get('smiles', '')
                inchi = structures.get('inchi', '')
                if swiss_smiles:
                    m = Chem.MolFromSmiles(swiss_smiles, sanitize=False)
                    found, miriam = vmh_search(m)
                elif inchi != 'InChI=none' and inchi:
                    m = Chem.MolFromInchi(
                        inchi, sanitize=False, removeHs=False)
                    found, miriam = vmh_search(m)
                else:
                    found, miriam = False, ''
        elif this_type == 'PubChem ID':
            base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
            endpoint = f"{base_url}{mols[idx]}/property/CanonicalSMILES,InChI/JSON"
            response = requests.get(endpoint, verify=False)
            if response.status_code != 200:
                found = False
                miriam = ''
            else:
                data = response.json()
                properties = data.get(
                    'PropertyTable', {}).get(
                    'Properties', [
                        {}])[0]
                pubchem_smiles = properties.get('CanonicalSMILES', '')
                inchi = properties.get('InChI', '')
                if pubchem_smiles:
                    m = Chem.MolFromSmiles(pubchem_smiles, sanitize=False)
                    found, miriam = vmh_search(m)
                elif inchi != 'InChI=none' and inchi:
                    m = Chem.MolFromInchi(
                        inchi, sanitize=False, removeHs=False)
                    found, miriam = vmh_search(m)
                else:
                    found, miriam = False, ''
        elif this_type == 'MDL Mol file':
            if nofile:
                temp_file_path = mols[idx]
            else:
                file_input = request.FILES.getlist(side)[file_idx]
                file_idx += 1  # Increment file index
                with NamedTemporaryFile(delete=False, suffix='.mol') as temp_file:
                    # Write the contents of the uploaded file to the temporary
                    # file
                    for chunk in file_input.chunks():
                        temp_file.write(chunk)

                    temp_file_path = temp_file.name
            try:
                mol_obj = Chem.MolFromMolFile(
                    temp_file_path, sanitize=False, removeHs=False)
                found, miriam = vmh_search(mol_obj)
            except VMHUnavailableError:
                raise
            except BaseException:
                found = False
                miriam = ''
        elif this_type == 'Draw':
            mol_file = unquote(mols[idx])
            with NamedTemporaryFile(delete=False, suffix='.mol', mode='w+') as temp_file:
                temp_file.write(mol_file)
                temp_file.flush()  # Ensure data is written to disk
            temp_file_path = temp_file.name
            try:
                mol_obj = Chem.MolFromMolFile(
                    temp_file_path, sanitize=False, removeHs=False)
                found, miriam = vmh_search(mol_obj)
            except VMHUnavailableError:
                raise
            except BaseException:
                found = False
                miriam = ''
        elif this_type == 'ChEBI ID' or this_type == 'ChEBI Name':
            mol_objs, errors, _ = any_to_mol([mols[idx]], [this_type], None, None)
            mol_obj, error = mol_objs[0], errors[0]
            if error:
                found = False
                miriam = ''
            else:
                found, miriam = vmh_search(mol_obj)
        elif this_type == 'SMILES':
            mol_obj = Chem.MolFromSmiles(mols[idx])
            if mol_obj is None:
                found, miriam = False, ''
            else:
                found, miriam = vmh_search(mol_obj)
        elif this_type == 'Saved':
            found = False
            miriam = ''
        else:
            raise Exception('Invalid type')
        founds.append(found)
        miriams.append(miriam)
    return founds, miriams


_FORMULA_ARROW = re.compile(r"\s*(<=>|->)\s*")


def decode_formula(formula):
    """
    Decodes a chemical formula into its constituent substrates and products with stoichiometry, along with their compartments.

    Input:
    - formula (str): The chemical reaction formula, e.g.
      ``2.0 h[c] + 2.0 r_o2[c] -> h2o2[c] + o2[c]`` or ``h2o[e] <=>``.

    Output:
    - (tuple): A tuple containing lists of substrates, products, their respective stoichiometries, and compartments.

    Raises ValueError when the formula has no reaction arrow.
    """
    # Function to extract the compartment identifier from a component
    def extract_compartment(component):
        match = re.search(r"\[(.*?)\]", component)
        return match.group(1) if match else ""

    # Function to parse each part into components, stoichiometry, and
    # compartments
    def parse_part(part):
        part = (part or '').strip()
        if not part:
            return [], [], []
        components = [component.strip()
                      for component in re.split(r"\s+\+\s+", part) if component.strip()]
        names = []
        stoichiometries = []
        compartments = []
        for component in components:
            # Extracting compartment
            compartment = extract_compartment(component)

            # Removing brackets and contents
            component = re.sub(r"\[.*?\]", "", component).strip()

            # Regex pattern to find stoichiometry (if separated by space) and
            # compound name
            match = re.match(r'(\d+\.?\d*)\s+(\S+)', component)
            if match:
                stoichiometry, name = match.groups()
                stoichiometry = float(stoichiometry)
            else:
                name = component
                stoichiometry = 1.0  # Default to 1.0 if no stoichiometry is provided

            names.append(name)
            stoichiometries.append(stoichiometry)
            compartments.append(compartment)
        return names, stoichiometries, compartments

    # Split the formula into substrates and products parts
    arrow = _FORMULA_ARROW.search(formula or '')
    if not arrow:
        raise ValueError(f"Not a reaction formula: {formula!r}")

    substrates, subs_stoich, subs_comps = parse_part(formula[:arrow.start()])
    products, prod_stoich, prods_comps = parse_part(formula[arrow.end():])
    direction = 'forward' if arrow.group(1) == '->' else 'bidirectional'
    return substrates, products, subs_stoich, prod_stoich, subs_comps, prods_comps, direction


def get_from_vmh(request):
    """
    Handles a POST request to fetch reaction data from the VMH database based on reaction abbreviation.

    Input:
    - request: The Django HTTP request object.

    Output:
    - (JsonResponse): A JsonResponse object containing reaction data or an error message.
    """
    if request.method != 'POST':
        return JsonResponse({'error': 'Invalid request'}, status=400)
    try:
        data = json.loads(request.body or b'{}')
    except ValueError:
        return JsonResponse({'error': 'Invalid request'}, status=400)
    reaction_abbreviation = str(data.get('reactionAbbreviation') or '').strip()
    if not reaction_abbreviation:
        return JsonResponse({"error": "Enter a VMH reaction abbreviation"})
    try:
        reaction_data = lookup_reaction(reaction_abbreviation)
    except VMHUnavailableError as exc:
        return JsonResponse(
            {"error": exc.user_message(f"fetch `{reaction_abbreviation}` from VMH")})
    if not reaction_data:
        return JsonResponse(
            {"error": "No reaction found with this abbreviation"})
    formula = reaction_data.get('formula', '')
    if not formula:
        return JsonResponse(
            {"error": f"Reaction `{reaction_abbreviation}` has no formula in VMH"})
    try:
        substrates, products, subs_sch, prod_sch, subs_comps, prod_comps, direction = decode_formula(
            formula)
    except ValueError:
        return JsonResponse(
            {"error": f"Could not read the VMH formula of `{reaction_abbreviation}`: {formula}"})
    # Construct the response data
    subs_types, prods_types = [
        'VMH' for _ in substrates], [
        'VMH' for _ in products]
    subsystem = reaction_data.get('subsystem', '')
    data = {
        'substrates': substrates,
        'products': products,
        'subs_sch': subs_sch,
        'prod_sch': prod_sch,
        'subs_comps': subs_comps,
        'prods_comps': prod_comps,
        'direction': direction,
        'subs_types': subs_types,
        'prods_types': prods_types,
        'subsystem': subsystem,
        'abbreviation': reaction_data.get('abbreviation') or reaction_abbreviation,
    }
    return JsonResponse(data)
