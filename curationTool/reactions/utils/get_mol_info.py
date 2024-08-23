from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import rdkit.Chem as Chem
from reactions.reaction_info import calculate_total_charge
from reactions.utils.to_mol import smiles_with_explicit_hydrogens

def get_mol_info(mols):
    formulas, charges, mol_file_strings = [], [], []
    for mol in mols:
        mol.UpdatePropertyCache(strict=False)
        mol_formula = CalcMolFormula(mol)
        mol_charge = calculate_total_charge([mol])
        try:
            mol_smiles = Chem.MolToSmiles(mol, allHsExplicit=True)
            mol_smiles = smiles_with_explicit_hydrogens(mol_smiles)
            mol = Chem.MolFromSmiles(mol_smiles)
            mol = Chem.AddHs(mol,explicitOnly=True)
        except:
            pass
        mol_block = Chem.MolToMolBlock(mol)
        mol_file_strings.append(mol_block)
        formulas.append(mol_formula)
        charges.append(mol_charge)
    return formulas, charges, mol_file_strings