# This provides a set of functions to convert various types of molecular identifiers
# (like abbreviations from different databases and MDL Mol files) into MOL files.
# It supports conversions from VMH, SwissLipids, CHeBI (not yet), and MDL Mol files,
# ensuring that all hydrogen atoms are explicitly represented in the
# resulting SMILES strings.

from rdkit import Chem
import os
from django.core.files.temp import NamedTemporaryFile
from urllib.parse import unquote
import requests
import uuid
from reactions_project.settings import MEDIA_ROOT, MEDIA_URL
from zeep import Client
from reactions.models import SavedMetabolite
from reactions.utils.vmh_api import find_metabolite_by_abbreviation

from django.conf import settings

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
    mol.UpdatePropertyCache(strict=False)

    # Add explicit hydrogens
    mol_with_h = Chem.AddHs(mol)

    # Convert back to SMILES with explicit hydrogens
    smiles_with_h = Chem.MolToSmiles(mol_with_h, allHsExplicit=True)
    return smiles_with_h


def vmh_to_mol(abbreviation):
    """
    Creates a mol file of a metabolite from the VMH database using its abbreviation.
    Falls back to /home/sandeep/ctf/mets/molFiles if not found in API.

    Input:
    - abbreviation (str): The abbreviation of the metabolite to query.

    Output:
    - (tuple): (RDKit Mol object or None, error message or None, metabolite name or '')
    """
    row = find_metabolite_by_abbreviation(abbreviation)
    if not row:
        # Fallback to local mol file
        MOL_FILE_PATH = settings.MOL_FILE_PATH
        local_path = os.path.join(MOL_FILE_PATH, f"{abbreviation}.mol")
        if os.path.exists(local_path):
            mol = Chem.MolFromMolFile(local_path, sanitize=False, removeHs=False)
            if mol:
                mol.UpdatePropertyCache(strict=False)
                return mol, None, abbreviation
            else:
                return None, f"Failed to parse .mol file for {abbreviation}", abbreviation
        return None, f"Metabolite {abbreviation} does not exist in VMH", abbreviation

    # If VMH returns results
    inchi_string = row.get('inchiString', '')
    smile = row.get('smile', '')
    name = row.get('fullName', '') or abbreviation

    if not smile and not inchi_string:
        return None, f"Metabolite {abbreviation} does not have SMILES or InChI string on VMH", name

    if inchi_string:
        mol = Chem.MolFromInchi(inchi_string, sanitize=False, removeHs=False)
        if mol:
            mol.UpdatePropertyCache(strict=False)
        else:
            return None, f"InChI parsing failed for {abbreviation}", name
    elif smile:
        if '[*]' in smile:
            return None, f"Metabolite {abbreviation} has a placeholder ([*]) in the SMILES string", name
        smile = smiles_with_explicit_hydrogens(smile)
        mol = Chem.MolFromSmiles(smile, sanitize=False, removeHs=False)
        if mol:
            mol.UpdatePropertyCache(strict=False)
        else:
            return None, f"SMILES parsing failed for {abbreviation}", name

    return mol, None, name


def pubchem_id_to_mol(pubchem_id):
    """
    Fetches the MOL rdkit object of a compound from PubChem using its ID.

    Input:
    - pubchem_id (str): The ID of the compound in PubChem.

    Output:
    - (tuple): A tuple containing the MOL object, an error message (if any), and the compound name.
    """
    base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
    endpoint = f"{base_url}{pubchem_id}/property/CanonicalSMILES,InChI,Title/JSON"
    response = requests.get(endpoint)

    if response.status_code != 200:
        return None, f"PubChem API returned error {response.status_code} for compound {pubchem_id}", ''

    data = response.json()
    try:
        properties = data['PropertyTable']['Properties'][0]
        smiles = properties.get('CanonicalSMILES', '')
        inchi = properties.get('InChI', '')
        name = properties.get('Title', '')

        inchi = None if inchi == 'InChI=none' else inchi

        if not smiles and not inchi:
            return None, f"Compound {pubchem_id} does not have SMILES or InChI in PubChem", ''

        if inchi:
            mol = Chem.MolFromInchi(inchi, sanitize=False, removeHs=False)
            mol.UpdatePropertyCache(strict=False)
            if not mol:
                inchi = None

        if smiles and not inchi:
            try:
                smiles = smiles_with_explicit_hydrogens(smiles)
            except BaseException:
                pass
            try:
                mol = Chem.MolFromSmiles(
                    smiles, sanitize=False, removeHs=False)
            except BaseException:
                mol = Chem.MolFromSmiles(smiles)
            mol.UpdatePropertyCache(strict=False)

        return mol, None, name

    except (KeyError, IndexError) as e:
        return None, f"Error parsing PubChem response for compound {pubchem_id}: {str(e)}", ''


def swisslipids_to_mol(id):
    """
    Fetches the MOL rdkit object of a metabolite from SwissLipids using its ID.

    Input:
    - id (str): The ID of the metabolite in SwissLipids.

    Output:
    - (tuple): A tuple containing the MOL object and an error message (if any).
    """
    base_url = 'https://www.swisslipids.org/api/index.php/entity/'
    endpoint = f"{base_url}{id}"
    response = requests.get(endpoint)
    if response.status_code != 200:
        return None, f"SwissLipids API returned error {response.status_code} for metabolite {id}", ''
    else:
        data = response.json()
        structures = data.get('structures', {})
        smiles = structures.get('smiles', '')
        inchi = structures.get('inchi', '')
        name = data.get('entity_name', '')
        mol_str = structures.get('structure', '')
        inchi = None if inchi == 'InChI=none' else inchi
        if not smiles and not inchi and not mol_str:
            return None, f"Metabolite {id} does not have SMILES or Inchi on SwissLipids", ''
        elif inchi:
            mol = Chem.MolFromInchi(inchi, sanitize=False, removeHs=False)
            mol.UpdatePropertyCache(strict=False)
            if not mol:
                inchi = None
        elif smiles and not inchi:
            # if '[*]' in smiles:
            #     return None, f"Metabolite {id} has a placeholder ([*]) in the SMILES string",''
            # else:
            try:
                smiles = smiles_with_explicit_hydrogens(smiles)
            except BaseException:
                pass
            try:
                mol = Chem.MolFromSmiles(
                    smiles, sanitize=False, removeHs=False)
            except BaseException:
                mol = Chem.MolFromSmiles(smiles)
            mol.UpdatePropertyCache(strict=False)
        elif mol_str:
            mol = Chem.MolFromMolBlock(mol_str)
            mol.UpdatePropertyCache(strict=False)
            smiles = Chem.MolToSmiles(mol)
            smiles = smiles_with_explicit_hydrogens(smiles)
            mol = Chem.MolFromSmiles(smiles)
    return mol, None, name


def draw_to_mol(inp):
    mol_file = unquote(inp)
    name = mol_file.split('\n')[0] + ' drawn'
    with NamedTemporaryFile(delete=False, suffix='.mol', mode='w+') as temp_file:
        temp_file.write(mol_file)
        temp_file.flush()  # Ensure data is written to disk
    mol = Chem.MolFromMolFile(temp_file.name, sanitize=False, removeHs=False)
    mol.UpdatePropertyCache(strict=False)
    if not mol:
        return None, f"{name} is not a valid MDL Mol file"

    return mol, None


def chebi_to_mol(chebi_id, name):
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
        params = {
            'search': chebi_id,
            'searchCategory': 'CHEBI NAME',
            'maximumResults': 10,
            'stars': 'ALL'}

        # Call getLiteEntity method
        response = client.service.getLiteEntity(**params)
        if not response or not response[0]:
            return None, f"Metabolite `{chebi_id}` does not exist in ChEBI", ''
        response = response[0]
        chebi_id = response['chebiId']
        name = response['chebiAsciiName']
    try:
        complete_entity = client.service.getCompleteEntity(chebiId=chebi_id)
    except BaseException:
        return None, f"Metabolite {chebi_id} does not exist in ChEBI", ''
    try:
        inchi = complete_entity['inchi']
    except BaseException:
        return None, f"Metabolite {chebi_id} does not have InChI on ChEBI", ''
    mol = Chem.MolFromInchi(inchi, sanitize=False, removeHs=False)
    mol.UpdatePropertyCache(strict=False)
    name = complete_entity['chebiAsciiName'] if not name else name
    try:
        mol = Chem.AddHs(mol)
    except BaseException:
        pass
    if mol:
        return mol, None, name
    else:
        return None, f"Invalid InChI {inchi}", ''

def saved_met_to_mol(id):
    """
    Fetches the MOL rdkit object of a metabolite from the saved metabolites database using its ID.

    Input:
    - id (int): The ID of the metabolite in the saved metabolites database.

    Output:
    - (tuple): A tuple containing the MOL object, an error message (if any), and the metabolite name.
    """
    try:
        saved_met = SavedMetabolite.objects.get(pk=id)
    except SavedMetabolite.DoesNotExist:
        return None, f"Saved metabolite with ID {id} not found", ''
    
    mol = Chem.MolFromInchi(saved_met.inchi, sanitize=False, removeHs=False)
    mol.UpdatePropertyCache(strict=False)
    if not mol:
        return None, f"Saved metabolite with ID {id} has an invalid InChI key", ''
    return mol, None, saved_met.name

def any_to_mol(mols, types, request, side='substrates'):
    """
    Converts a list of molecular identifiers from various sources to RDKit molecule objects.

    Input:
    - mols (list): A list of molecular identifiers.
    - type (list): A list of types corresponding to each molecular identifier.
    - request: The Django request object (used for file handling).
    - side (str): A string indicating whether the molecules are substrates or products.

    Output:
    - (tuple): A tuple containing a list of RDKIT objects and a list of error messages (None if no error).
    """
    mol_objs = []
    errors = []
    names = []
    file_idx = 0  # Initialize file index
    for idx, type in enumerate(types):
        if type == 'VMH':
            mol, error, name = vmh_to_mol(mols[idx])
        elif type == 'ChEBI ID':
            mol, error, name = chebi_to_mol(mols[idx], name=False)
        elif type == 'ChEBI Name':
            mol, error, name = chebi_to_mol(mols[idx], name=True)
        elif type == 'PubChem ID':
            mol, error, name = pubchem_id_to_mol(mols[idx])
        elif type == 'SwissLipids':
            mol, error, name = swisslipids_to_mol(mols[idx])
        elif type == 'Saved':
            mol, error, name = saved_met_to_mol(mols[idx])
        elif type == 'MDL Mol file':
            if request:
                try:
                    file_input = request.FILES.getlist(side)[file_idx]
                except BaseException:
                    file_input = request.FILES.getlist('file')[file_idx]
                if not file_input:
                    mol, error, name = None, f"No file was uploaded", ''
                else:
                    original_filename = file_input.name
                    unique_filename = f"{uuid.uuid4()}"
                    # Define the permanent file save path
                    save_path = os.path.join(
                        MEDIA_ROOT, 'mol_files', unique_filename)
                    with open(save_path, 'wb+') as permanent_file:
                        # Write the contents of the uploaded file to the
                        # permanent file
                        for chunk in file_input.chunks():
                            permanent_file.write(chunk)
            else:
                save_path = mols[idx].split(MEDIA_URL)[1]
                save_path = os.path.join(MEDIA_ROOT, save_path)
            mol = Chem.MolFromMolFile(
                save_path, sanitize=False, removeHs=False)
            if not mol:
                mol, error = None, f"Invalid MDL Mol file {original_filename}"
            else:
                mol.UpdatePropertyCache(strict=False)
            if request:
                os.remove(save_path)
            error = None if mol else f"Invalid MDL Mol file {original_filename}"
            file_idx += 1
            name = ''
        elif type == 'Draw':
            mol, error = draw_to_mol(mols[idx])
            name = ''
        else:
            raise Exception('Invalid type')
        try:
            mol.UpdatePropertyCache(strict=False)
            mol = Chem.AddHs(mol)
        except BaseException:
            pass
        mol_objs.append(mol)
        errors.append(error)
        names.append(name)
    return mol_objs, errors, names
