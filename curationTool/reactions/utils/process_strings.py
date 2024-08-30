from rdkit import Chem
from urllib.parse import unquote
import time

def get_mol_names(mols_list, mols_types):
    mols_names = [unquote(mol).split('\n')[0].split(' ')[1] if mol_type == 'Draw' else mol for mol, mol_type in zip(mols_list, mols_types)]
    mols_names = [mol.split('.mol')[0] if mol_type == 'MDL Mol file' and 'mol' in mol else mol for mol, mol_type in zip(mols_names, mols_types)]
    mols_names = [mol.split('.sdf')[0] if mol_type == 'MDL Mol file' and 'sdf' in mol else mol for mol, mol_type in zip(mols_names, mols_types)]
    return mols_names
def construct_reaction_string(substrate_smiles, substrate_stoichiometry, product_smiles, product_stoichiometry):
    """
    Constructs a reaction string in SMILES format from the given substrates and products with their stoichiometries.

    Input:
    - substrate_smiles (list): List of SMILES strings representing the substrates in the reaction.
    - substrate_stoichiometry (list): List of integers representing the stoichiometry of each substrate.
    - product_smiles (list): List of SMILES strings representing the products in the reaction.
    - product_stoichiometry (list): List of integers representing the stoichiometry of each product.

    Output:
    - (str): A reaction string formatted in SMILES, with substrates and products separated by '>>'.
    """
    # Create the substrate part of the reaction string, ensuring individual molecules are separated
    substrate_part = '.'.join([s for s, sm in zip(substrate_smiles, substrate_stoichiometry) for _ in range(sm)])
    
    # Create the product part of the reaction string, ensuring individual molecules are separated
    product_part = '.'.join([s for s, sm in zip(product_smiles, product_stoichiometry) for _ in range(sm)])
    
    # Combine the substrate and product parts with '>>'
    reaction_string = f"{substrate_part}>>{product_part}"
    
    return reaction_string

def construct_reaction_rxnfile(substrate_mols, substrate_stoichiometry, product_mols, product_stoichiometry, substrates_names, products_names):
    """
    Constructs an RXN file from given substrates and products (RDKit MOL objects) with their stoichiometries.
    This function repeats molecules according to their stoichiometry for accurate representation.

    :param substrate_mols: List of RDKit MOL objects for substrates.
    :param substrate_stoichiometry: List of integers for substrate stoichiometries.
    :param product_mols: List of RDKit MOL objects for products.
    :param product_stoichiometry: List of integers for product stoichiometries.
    :param substrates_names: List of names (strings) for each substrate.
    :param products_names: List of names (strings) for each product.
    :return: Path to the generated RXN file.
    """
    

    rxn = Chem.rdChemReactions.ChemicalReaction()

    # Add substrates to the reaction, repeated according to stoichiometry
    for mol, stoich, name in zip(substrate_mols, substrate_stoichiometry, substrates_names):
        for _ in range(stoich):  # Repeat molecule addition based on stoichiometry
            mol_with_name = Chem.RWMol(mol)
            mol_with_name.SetProp("_Name", name)
            rxn.AddReactantTemplate(mol_with_name)
    
    # Add products to the reaction, repeated according to stoichiometry
    for mol, stoich, name in zip(product_mols, product_stoichiometry, products_names):
        for _ in range(stoich):  # Repeat molecule addition based on stoichiometry
            mol_with_name = Chem.RWMol(mol)
            mol_with_name.SetProp("_Name", name)
            rxn.AddProductTemplate(mol_with_name)
    
    # Convert the reaction to an RXN block string
    rxn_block = Chem.rdChemReactions.ReactionToRxnBlock(rxn)

    # Define the path for saving the RXN file
    rxn_file_path = 'temp_rxn_file.rxn'
    
    # Write the RXN block to the file
    with open(rxn_file_path, "w") as f:
        f.write(rxn_block)

    return rxn_file_path