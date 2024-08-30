# This includes a collection of functions for handling chemical compounds and reactions.
# It provides functionality to convert various molecular identifiers to SMILES strings, fetch information from VMH,
# and process chemical formulas. The functions are designed to integrate with a Django web application, handling both HTTP requests and file operations.

import requests
from rdkit import Chem
from reactions.utils.to_smiles import smiles_with_explicit_hydrogens
from django.core.files.temp import NamedTemporaryFile
from urllib.parse import quote, unquote
import json
import re
from django.http import JsonResponse
from reactions.utils.to_smiles import any_to_smiles
from reactions.utils.to_mol import any_to_mol

def is_name_in_vmh(name):
    """
    Checks if a molecule name is found in the VMH database.

    Input:
    - name (str): The molecule name.

    Output:
    - (bool): True if found, False otherwise.
    """
    BASE_URL = 'https://www.vmh.life/'
    encoded_name = quote(name)
    endpoint = f"{BASE_URL}_api/metabolites/?fullName={encoded_name}"
    response = requests.get(endpoint, verify=False)
    if response.status_code != 200:
        return False
    elif response.json().get('count', 0) == 0:
        return False
    else: 
        match = False
        for result in response.json().get('results', []):
            if result.get('fullName', '').lower() == name.lower():
                match = True
                break
        return match

def any_to_vmh(mols,types,smiles):
    """
    Converts a list of molecules to their VMH database abbreviations.

    Input:
    - mols (list): List of molecule identifiers.
    - types (list): List of types corresponding to each molecule identifier.
    - smiles (list): List of SMILES strings for the molecules.

    Output:
    - (list): List of VMH database abbreviations or 'error' for unsuccessful conversions.
    """
    BASE_URL = 'https://www.vmh.life/'
    mols_list = []
    for idx,mol in enumerate(mols):
        if types[idx] == 'VMH':
            mols_list.append(mol)
        else:
            m = Chem.MolFromSmiles(str(smiles[idx]), sanitize=False)
            inchi = Chem.MolToInchi(m)
            encoded_inchi = quote(inchi)
            endpoint = f"{BASE_URL}_api/metabolites/?inchiString={encoded_inchi}"
            # Make the GET request
            response = requests.get(endpoint, verify=False)
            if response.status_code != 200:
                mols_list.append('error')
            if response.json().get('count', 0) == 0:
                mols_list.append('error')
            else:
                abbr = response.json().get('results',[[]])[0].get('abbreviation','')
                mols_list.append(abbr)
    return mols_list


def check_reaction_vmh(substrates, products, subs_sch, prod_sch,substrates_types,products_types,subs_mols,prod_mols,direction, subsystem, subs_comps, prods_comps):
    """
    Checks if a chemical reaction is found in the VMH database based on substrates and products.

    Input:
    - substrates, products (list): Lists of substrate and product molecules.
    - subs_sch, prod_sch (list): Lists of stoichiometries for substrates and products.
    - substrates_types, products_types (list): Lists of types for substrates and products.
    - subs_smiles, prod_smiles (list): Lists of SMILES strings for substrates and products.

    Output:
    - (dict): Dictionary containing reaction information if found, otherwise an 'error' or 'not found' message.
    """
    # Checking lo
    BASE_URL = 'https://www.vmh.life/'
    subs_smiles = any_to_smiles(subs_mols,['MDL Mol file' for _ in range(len(subs_mols))],request=None,side='substrates')
    prod_smiles = any_to_smiles(prod_mols,['MDL Mol file' for _ in range(len(prod_mols))],request=None,side='products')
    if None in subs_smiles[0] or None in prod_smiles[0]:
         return {"found": False}
    subs_smiles = subs_smiles[0]
    prod_smiles = prod_smiles[0]
    # Use the first substrate and product for the API call
    substrates = any_to_vmh(substrates,substrates_types,subs_smiles)
    products = any_to_vmh(products,products_types,prod_smiles)
    subs_sch = [float(s) for s in subs_sch]
    prod_sch = [float(p) for p in prod_sch]
    found = False
    miriams = []
    formulas = []  
    for ir in range(len(substrates)):
        if found and not similar:
            break
        for ip in range(len(products)):
            if found and not similar:
                break
            abbrReactant = substrates[ir]
            abbrProduct = products[ip]
        
            # Construct the API endpoint URL
            endpoint = f"{BASE_URL}_api/reactionfromreactandproduct/{abbrReactant}/{abbrProduct}/"
            
            # Make the GET request
            response = requests.get(endpoint, verify=False)

            if response.status_code != 200:
                return json.dumps({"error": f"Failed to fetch data from API. Status code: {response.status_code}"})

            # Parse the API response
            data = response.json()
            for result in data.get("results", []):
                if found and not similar:
                    break
                formula = result.get("rxn", {}).get("formula", "")
                this_subsystem = result.get("rxn", {}).get("subsystem", "")
                # Split the formula into substrates (left) and products (right) sides
                this_substrates, this_products, this_subs_stoich, this_prod_stoich, this_subs_comps, this_prods_comps, this_direction = decode_formula(formula)
                
                left_mole_to_stoich = dict(zip(this_substrates, this_subs_stoich))
                right_mole_to_stoich = dict(zip(this_products, this_prod_stoich))
                
                substrate_to_stoich = dict(zip(substrates, subs_sch))
                product_to_stoich = dict(zip(products, prod_sch))
                # Mapping molecules to compartments for both sides of the equation
                left_mole_to_comps = dict(zip(this_substrates, this_subs_comps))
                right_mole_to_comps = dict(zip(this_products, this_prods_comps))

                # Mapping molecules to compartments for the comparison sets
                substrate_to_comps = dict(zip(substrates, subs_comps))
                product_to_comps = dict(zip(products, prods_comps))

                # Check if substrates' compartments match
                substrates_comps_match = all(molecule in left_mole_to_comps and
                                            substrate_to_comps[molecule] == left_mole_to_comps[molecule]
                                            for molecule in substrates) and all(molecule in substrate_to_comps and
                                            left_mole_to_comps[molecule] == substrate_to_comps[molecule]
                                            for molecule in this_substrates)

                # Check if products' compartments match
                products_comps_match = all(molecule in right_mole_to_comps and
                                            product_to_comps[molecule] == right_mole_to_comps[molecule]
                                            for molecule in products) and all(molecule in product_to_comps and
                                            right_mole_to_comps[molecule] == product_to_comps[molecule]
                                            for molecule in this_products)

                # Update comps_match based on substrates and products compartments match
                comps_match = substrates_comps_match and products_comps_match

                substrates_match = all(molecule in left_mole_to_stoich and
                       substrate_to_stoich[molecule] == left_mole_to_stoich[molecule]
                       for molecule in substrates) and all(molecule in substrate_to_stoich and
                       left_mole_to_stoich[molecule] == substrate_to_stoich[molecule]
                       for molecule in this_substrates)

                products_match = all(molecule in right_mole_to_stoich and
                                    product_to_stoich[molecule] == right_mole_to_stoich[molecule]
                                    for molecule in products) and all(molecule in product_to_stoich and
                       right_mole_to_stoich[molecule] == product_to_stoich[molecule]
                       for molecule in this_products)
                if substrates_match and products_match:
                    found = True
                    img_code = result.get("rxn", {}).get("miriam")
                    miriams.append(img_code)
                    formulas.append(formula)
                    similar = False if this_subsystem == subsystem and this_direction == direction and comps_match else True

    if found:
        return {
            "found": found,
            "similar": similar,
            "url": miriams[0] if miriams else None,
            "formula": formulas[0] if formulas else None
        }
    else:
        return {"found": found,
                "similar": False}
    
def search_vmh(mol,return_abbr=False,return_name=False):
    """
    Searches the VMH database for a molecule.

    Input:
    - mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object.
                                                                                
    Output:
    - (tuple): A tuple containing a boolean indicating if found, and the miriam ID if found.
    """
    BASE_URL = 'https://www.vmh.life/'
    found = False
    smiles = Chem.MolToSmiles(mol)
    smiles = smiles_with_explicit_hydrogens(smiles)
    m = Chem.MolFromSmiles(str(smiles), sanitize=False)
    try:
        Chem.SanitizeMol(m)
    except:
        pass
    try:
        Chem.AssignStereochemistry(m)
    except:
        pass 
    inchi = Chem.MolToInchi(m)
    encoded_inchi = quote(inchi)
    if inchi.strip() == '':
        found = False
        miriam = None
        abbr = ''
        name = ''
    else:
        endpoint = f"{BASE_URL}_api/metabolites/?inchiString__icontains={encoded_inchi}"
        # Make the GET request
        response = requests.get(endpoint, verify=False)
        if response.status_code != 200 or response.json().get('count', 0) == 0: 
            encoded_inchi = quote(Chem.MolToInchi(mol))
            endpoint = f"{BASE_URL}_api/metabolites/?inchiString={encoded_inchi}"
            response = requests.get(endpoint, verify=False)
            if response.status_code != 200 or response.json().get('count', 0) == 0:
                if return_abbr and return_name:
                    return found,'','',''
                return found,''
        miriam = response.json().get('results',[[]])[0].get('miriam','')
        abbr = response.json().get('results',[[]])[0].get('abbreviation','')
        name = response.json().get('results',[[]])[0].get('fullName','')
    if miriam:
        found = True
    if return_abbr and return_name:
        return found,miriam, abbr, name
    elif return_abbr:
        return found,abbr
    elif return_name:
        return found,name
    return found,miriam

def get_vmh_miriam(abbr):
    """
    Retrieves the MIRIAM ID for a molecule abbreviation from the VMH database.

    Input:
    - abbr (str): The molecule abbreviation.

    Output:
    - (str): The MIRIAM ID for the molecule.
    """
    BASE_URL = 'https://www.vmh.life/'
    encoded_abbr = quote(abbr)
    endpoint = f"{BASE_URL}_api/metabolites/?abbreviation={encoded_abbr}"
    response = requests.get(endpoint, verify=False)
    miriam = response.json().get('results',[[]])[0].get('miriam','')
    return miriam

def search_metabolites_vmh(mols, type, request, side='substrates',nofile=False,return_abbr=False):
    """
    Searches for metabolites in the VMH database based on a list of molecule identifiers.

    Input:
    - mols (list): List of molecule identifiers.
    - type (list): Corresponding types for each molecule identifier.
    - request: The Django request object (for file handling).
    - side (str): Indicating 'substrates' or 'products'.

    Output:
    - (tuple): Tuple of two lists - one indicating if each molecule was found, and another with MIRIAM IDs.
    """
    founds = []
    miriams = []
    file_idx = 0  # Initialize file index
    for idx,type in enumerate(type):
        if type == 'VMH':
            found = True
            miriam = get_vmh_miriam(mols[idx])
        elif type == 'SwissLipids':
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
                    found,miriam = search_vmh(m,return_abbr=return_abbr)
                elif inchi != 'InChI=none' and inchi:
                    m = Chem.MolFromInchi(inchi,sanitize=False,removeHs=False)
                    found,miriam = search_vmh(m,return_abbr=return_abbr)
                else:
                    found,miriam = False,''
        elif type == 'PubChem ID':
            base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
            endpoint = f"{base_url}{mols[idx]}/property/CanonicalSMILES,InChI/JSON"
            response = requests.get(endpoint, verify=False)
            if response.status_code != 200:
                found = False
                miriam = ''
            else:
                data = response.json()
                properties = data.get('PropertyTable', {}).get('Properties', [{}])[0]
                pubchem_smiles = properties.get('CanonicalSMILES', '')
                inchi = properties.get('InChI', '')
                if pubchem_smiles:
                    m = Chem.MolFromSmiles(pubchem_smiles, sanitize=False)
                    found, miriam = search_vmh(m, return_abbr=return_abbr)
                elif inchi != 'InChI=none' and inchi:
                    m = Chem.MolFromInchi(inchi, sanitize=False, removeHs=False)
                    found, miriam = search_vmh(m, return_abbr=return_abbr)
                else:
                    found, miriam = False, ''
        elif type == 'MDL Mol file':
            if nofile:
                temp_file_path = mols[idx]
            else:
                file_input = request.FILES.getlist(side)[file_idx]
                file_idx += 1  # Increment file index
                with NamedTemporaryFile(delete=False, suffix='.mol') as temp_file:
                    # Write the contents of the uploaded file to the temporary file
                    for chunk in file_input.chunks():
                        temp_file.write(chunk)

                    temp_file_path = temp_file.name
            try:
                mol_obj = Chem.MolFromMolFile(temp_file_path,sanitize=False,removeHs=False)
                found,miriam = search_vmh(mol_obj,return_abbr=return_abbr)
            except:
                found = False
                miriam = ''
        elif type == 'Draw':
            mol_file = unquote(mols[idx])
            with NamedTemporaryFile(delete=False, suffix='.mol',mode='w+') as temp_file:
                temp_file.write(mol_file)
                temp_file.flush()  # Ensure data is written to disk
            temp_file_path = temp_file.name
            try:
                mol_obj = Chem.MolFromMolFile(temp_file_path,sanitize=False,removeHs=False)
                found,miriam = search_vmh(mol_obj,return_abbr=return_abbr)
            except:
                found = False
                miriam = ''
        elif type == 'ChEBI ID' or type == 'ChEBI Name':
            mol_objs, errors,_ = any_to_mol([mols[idx]],[type],None,None)
            mol_obj, error = mol_objs[0], errors[0]
            if error:
                found = False
                miriam = ''
            else:
                found,miriam = search_vmh(mol_obj,return_abbr=return_abbr)
        else:
            raise Exception('Invalid type')
        founds.append(found)
        miriams.append(miriam)
    return founds, miriams

def decode_formula(formula):
    """
    Decodes a chemical formula into its constituent substrates and products with stoichiometry, along with their compartments.

    Input:
    - formula (str): The chemical reaction formula.

    Output:
    - (tuple): A tuple containing lists of substrates, products, their respective stoichiometries, and compartments.
    """
    # Function to extract the compartment identifier from a component
    def extract_compartment(component):
        match = re.search(r"\[(.*?)\]", component)
        return match.group(1) if match else ""

    # Function to parse each part into components, stoichiometry, and compartments
    def parse_part(part):
        components = part.split(' + ')
        names = []
        stoichiometries = []
        compartments = []
        for component in components:
            # Extracting compartment
            compartment = extract_compartment(component)

            # Removing brackets and contents
            component = re.sub(r"\[.*?\]", "", component).strip()

            # Regex pattern to find stoichiometry (if separated by space) and compound name
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
    parts = formula.split(' -> ') if ' -> ' in formula else formula.split(' <=> ')

    substrates, subs_stoich, subs_comps = parse_part(parts[0])
    products, prod_stoich, prods_comps = parse_part(parts[1])
    direction = 'forward' if '->' in formula else 'bidirectional'
    return substrates, products, subs_stoich, prod_stoich, subs_comps, prods_comps, direction

def get_from_vmh(request):
    """
    Handles a POST request to fetch reaction data from the VMH database based on reaction abbreviation.

    Input:
    - request: The Django HTTP request object.

    Output:
    - (JsonResponse): A JsonResponse object containing reaction data or an error message.
    """
    BASE_URL = 'https://www.vmh.life/'
    if request.method == 'POST':
        data = json.loads(request.body)  # Parse JSON data from request body
        reaction_abbreviation = data.get('reactionAbbreviation', '')
        endpoint = f"{BASE_URL}_api/reactions/?abbreviation={reaction_abbreviation}"
        response = requests.get(endpoint, verify=False)
        data = response.json()
        if response.status_code != 200:
            return JsonResponse({"error": f"Failed to fetch data from API. Status code: {response.status_code}"})
        if data.get('count', 0) == 0:
            return JsonResponse({"error": "No reaction found with this abbreviation"})
        formula = data.get('results', [])[0].get('formula', '')
        substrates, products, subs_sch, prod_sch, subs_comps,prod_comps, direction= decode_formula(formula)
        # Construct the response data
        subs_types, prods_types = ['VMH' for _ in substrates], ['VMH' for _ in products]
        subsystem = data.get('results', [])[0].get('subsystem', '')
        data = {
            'substrates': substrates,
            'products': products,
            'subs_sch': subs_sch,
            'prod_sch': prod_sch,
            'subs_comps': subs_comps,
            'prods_comps': prod_comps,
            'direction': direction,
            'subs_types':subs_types,
            'prods_types':prods_types,
            'subsystem': subsystem
        }
        print(data)
        return JsonResponse(data)
    else:
        # Handle non-POST requests if necessary
        return JsonResponse({'error': 'Invalid request'}, status=400)

