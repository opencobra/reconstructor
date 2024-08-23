# This provides a set of functions to convert various types of molecular identifiers 
# (like abbreviations from different databases and MDL Mol files) into SMILES strings. 
# It supports conversions from VMH, SwissLipids, CHeBI, and MDL Mol files, 
# ensuring that all hydrogen atoms are explicitly represented in the resulting SMILES strings.

from rdkit import Chem
import os
from django.core.files.temp import NamedTemporaryFile
from urllib.parse import quote, unquote
import requests
from reactions_project.settings import MEDIA_ROOT,MEDIA_URL
from zeep import Client

def smiles_with_explicit_hydrogens(smiles):
    """
    Converts a SMILES string to a version with all hydrogen atoms explicitly represented.
    
    Input:
    - smiles (str): The SMILES string to be processed.

    Output:
    - (str): The SMILES string with all hydrogen atoms explicitly represented.
    """
    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    # Add explicit hydrogens
    try:
        mol_with_h = Chem.AddHs(mol)
    except:
        mol_with_h = mol
    
    # Convert back to SMILES with explicit hydrogens
    try:
        smiles_with_h = Chem.MolToSmiles(mol_with_h, allHsExplicit=True)
    except:
        return smiles
    return smiles_with_h

def vmh_to_smiles(abbreviation):
    """
    Fetches the SMILES representation of a metabolite from the VMH database using its abbreviation.

    Input:
    - abbreviation (str): The abbreviation of the metabolite to query in the VMH database.

    Output:
    - (tuple): A tuple containing the SMILES string and an error message (if any).
    """
    BASE_URL = 'https://www.vmh.life/'
    encoded_abbr = quote(abbreviation)
    endpoint = f"{BASE_URL}_api/metabolites/?abbreviation={encoded_abbr}"
    response = requests.get(endpoint, verify=False)
    
    if response.status_code != 200:
        return None, f"VMH API returned error {response.status_code} for metabolite {abbreviation}"
    
    data = response.json()
    res = data.get('results', [])
    if len(res) == 0:
        return None, f"Metabolite {abbreviation} does not exist in VMH"
    
    smile = res[0].get('smile', '')
    if not smile:
        return None, f"Metabolite {abbreviation} does not have SMILES on VMH"
    smile = smiles_with_explicit_hydrogens(smile)
    return smile, None

def chebi_to_smiles(chebi_id,name):
    """
    Fetches the MOL rdkit object of a metabolite from ChEBI using its ID.
    input:
    - chebi_id (str): The ID of the metabolite in ChEBI.
    - name (bool): If True, the chebi ID is expected to be a name, otherwise it is expected to be an ID.
    Output:
    - (tuple): A tuple containing - the MOL object,- an error message (if any) and - the name of the metabolite.
    """
    wsdl = 'https://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl'
    client = Client(wsdl)

    if name:
        # Setup parameters to search for "group"
        params = {'search': chebi_id, 'searchCategory': 'CHEBI NAME', 'maximumResults': 10, 'stars': 'ALL'}

        # Call getLiteEntity method
        response = client.service.getLiteEntity(**params)
        response = response[0]
        chebi_id = response['chebiId']
        
    complete_entity = client.service.getCompleteEntity(chebiId=chebi_id)
    try:
        smiles = complete_entity['smiles']
    except:
        return None, f"Metabolite {chebi_id} does not have InChI on ChEBI"
    try:
        smiles = smiles_with_explicit_hydrogens(smiles)
    except:
        pass
    return smiles, None
    
def mdl_to_smiles(temp_file_path,file_input_name=None):
    """
    Converts a molecule from an MDL Mol file to a SMILES string with all hydrogens explicitly represented.

    Input:
    - file_input: The MDL Mol file input.

    Output:
    - (tuple): A tuple containing the SMILES string and an error message (if any).
    """
    try:
        m = Chem.MolFromMolFile(temp_file_path,sanitize=False,removeHs=False)
        smiles = Chem.MolToSmiles(m,allHsExplicit=True,allBondsExplicit=True)
        smiles = smiles_with_explicit_hydrogens(smiles)
    except Exception as e:
        print(e)
        return None, f"{file_input_name} is not a valid MDL Mol file"
    # Clean up: Remove the temporary file
    if os.path.exists(temp_file_path):
        os.remove(temp_file_path)
    return smiles, None
def swisslipids_to_smiles(id):
    """
    Fetches the SMILES representation of a metabolite from SwissLipids using its ID.

    Input:
    - id (str): The ID of the metabolite in SwissLipids.

    Output:
    - (tuple): A tuple containing the SMILES string and an error message (if any).
    """
    base_url = 'https://www.swisslipids.org/api/index.php/entity/'
    endpoint = f"{base_url}{id}"
    response = requests.get(endpoint)
    if response.status_code != 200:
        return None, f"SwissLipids API returned error {response.status_code} for metabolite {id}"
    else:
        data = response.json()
        structures = data.get('structures', [])
        smiles = structures.get('smiles', '')
        inchi = structures.get('inchi', '')
        if smiles:
            found_smiles = smiles
            found_smiles = smiles_with_explicit_hydrogens(found_smiles)
        elif inchi != 'InChI=none' and inchi:
            m = Chem.MolFromInchi(inchi,sanitize=False,removeHs=False)
            found_smiles = Chem.MolToSmiles(m, allHsExplicit=True)
            found_smiles = smiles_with_explicit_hydrogens(found_smiles)
        else:
            return None, f"Metabolite {id} does not have SMILES or Inchi on SwissLipids"
    return found_smiles, None

def pubchem_id_to_smiles(pubchem_id):
    """
    Fetches the SMILES representation of a compound from PubChem using its ID.

    Input:
    - pubchem_id (str): The ID of the compound in PubChem.

    Output:
    - (tuple): A tuple containing the SMILES string and an error message (if any).
    """
    base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
    endpoint = f"{base_url}{pubchem_id}/property/CanonicalSMILES,InChI/JSON"
    response = requests.get(endpoint)
    
    if response.status_code != 200:
        return None, f"PubChem API returned error {response.status_code} for compound {pubchem_id}"
    
    data = response.json()
    try:
        properties = data['PropertyTable']['Properties'][0]
        smiles = properties.get('CanonicalSMILES', '')
        inchi = properties.get('InChI', '')
        
        if smiles:
            found_smiles = smiles
            found_smiles = smiles_with_explicit_hydrogens(found_smiles)
        elif inchi != 'InChI=none' and inchi:
            m = Chem.MolFromInchi(inchi, sanitize=False, removeHs=False)
            found_smiles = Chem.MolToSmiles(m, allHsExplicit=True)
            found_smiles = smiles_with_explicit_hydrogens(found_smiles)
        else:
            return None, f"Compound {pubchem_id} does not have SMILES or InChI in PubChem"
        
        return found_smiles, None
    
    except (KeyError, IndexError) as e:
        return None, f"Error parsing PubChem response for compound {pubchem_id}: {str(e)}"


def draw_to_smiles(inp):
    mol_file = unquote(inp)
    name = mol_file.split('\n')[0]+ ' drawn'
    with NamedTemporaryFile(delete=False, suffix='.mol',mode='w+') as temp_file:
        temp_file.write(mol_file)
        temp_file.flush()  # Ensure data is written to disk
    
    smile, error = mdl_to_smiles(temp_file.name,file_input_name=name)
    return smile, error

def any_to_smiles(mols, types,request,side='substrates'):
    """
    Converts a list of molecular identifiers from various sources to SMILES strings.

    Input:
    - mols (list): A list of molecular identifiers.
    - type (list): A list of types corresponding to each molecular identifier.
    - request: The Django request object (used for file handling).
    - side (str): A string indicating whether the molecules are substrates or products.

    Output:
    - (tuple): A tuple containing a list of SMILES strings and a list of error messages.
    """
    smiles = []
    errors = []
    file_idx = 0  # Initialize file index
    for idx,type in enumerate(types):
        if type == 'VMH':
            smile, error = vmh_to_smiles(mols[idx])
        elif type == 'ChEBI ID':
            smile, error = chebi_to_smiles(mols[idx],name=False)
        elif type == 'ChEBI Name':
            smile, error = chebi_to_smiles(mols[idx],name=True)
        elif type == 'PubChem ID':
            smile, error = pubchem_id_to_smiles(mols[idx])
        elif type == 'SwissLipids':
            smile, error = swisslipids_to_smiles(mols[idx])
        elif type == 'MDL Mol file':
            if request:
                file_input = request.FILES.getlist(side)[file_idx]
                with NamedTemporaryFile(delete=False, suffix='.mol') as temp_file:
                    # Write the contents of the uploaded file to the temporary file
                    for chunk in file_input.chunks():
                        temp_file.write(chunk)
                    temp_file.flush()
                    temp_file_path = temp_file.name
                file_idx += 1  # Increment file index
                smile, error = mdl_to_smiles(temp_file_path,file_input_name=file_input)
            else:
                try:
                    mol_block = Chem.MolToMolBlock(mols[idx])
                except:
                    save_path = mols[idx].split(MEDIA_URL)[1]
                    save_path = os.path.join(MEDIA_ROOT, save_path)
                    mol = Chem.MolFromMolFile(save_path, sanitize=False, removeHs=False)
                    mol_block = Chem.MolToMolBlock(mol)
                with NamedTemporaryFile(delete=False, suffix='.mol',mode='w+') as temp_file:
                    temp_file.write(mol_block)
                    temp_file.flush()
                smile, error = mdl_to_smiles(temp_file.name,file_input_name='Temp Mol file')
        elif type == 'Draw':
            smile, error = draw_to_smiles(mols[idx])
        else:
            raise Exception('Invalid type')
        smiles.append(smile)
        errors.append(error)
    return smiles, errors
