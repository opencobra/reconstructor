from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from rdkit.Chem.Descriptors import MolWt
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from reactions.reaction_info import calculate_total_charge
from reactions.utils.to_mol import smiles_with_explicit_hydrogens


def get_mol_info(mols):
    metabolite_formulas, metabolite_charges, mol_file_strings = [], [], []
    stereo_counts, stereo_locations_list = [], []
    smiles, inchi_keys, inchis, mol_weights = [], [], [], []

    for mol in mols:
        # Update property cache and calculate formula and charge
        mol.UpdatePropertyCache(strict=False)
        mol_formula = CalcMolFormula(mol)
        mol_charge = calculate_total_charge([mol])

        try:
            # Generate a SMILES string with explicit hydrogens
            mol_smiles = Chem.MolToSmiles(mol, allHsExplicit=True)
            mol_smiles = smiles_with_explicit_hydrogens(mol_smiles)
            mol = Chem.MolFromSmiles(mol_smiles)
            mol = Chem.AddHs(mol, explicitOnly=True)
        except Exception as e:
            print(f"Error processing molecule: {e}")
            continue
        smiles.append(mol_smiles)
        inchi_keys.append(Chem.MolToInchiKey(mol))
        inchis.append(Chem.MolToInchi(mol))
        mol_weights.append(MolWt(mol))
        # Embed 3D coordinates 
        try:
            if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                raise ValueError("Embedding failed for molecule.")
            AllChem.UFFOptimizeMolecule(mol)  # Optimize 3D structure
        except Exception as e:
            print(f"Error embedding 3D coordinates: {e}")
            continue

        # Convert to MolBlock for 3Dmol.js visualization
        mol_block = Chem.MolToMolBlock(mol)
        mol_file_strings.append(mol_block)
        metabolite_formulas.append(mol_formula)
        metabolite_charges.append(mol_charge)

        # Detect chiral centers (including unassigned centers)
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        unspecified_centers = [
            center for center in chiral_centers if center[1] == '?']

        # Count and store stereochemical information
        stereo_count = len(unspecified_centers)
        stereo_locations = [center[0] for center in unspecified_centers]

        stereo_counts.append(stereo_count)
        stereo_locations_list.append(stereo_locations)

    # Return gathered information
    return {
        'smiles': smiles,
        'inchi_keys': inchi_keys,
        'inchis': inchis,
        'mol_weights': mol_weights,
        'metabolite_formulas': metabolite_formulas,
        'metabolite_charges': metabolite_charges,
        'metabolite_mol_file_strings': mol_file_strings,
        'stereo_counts': stereo_counts,
        'stereo_locations_list': stereo_locations_list
    }
