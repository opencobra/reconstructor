from rdkit.Chem import rdChemReactions
from rdkit import Chem
from collections import defaultdict
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.rdchem import PeriodicTable
from rdkit.Chem.rdchem import GetPeriodicTable
import json 

def calculate_total_charge(molecules):
    """Calculate the total charge of a set of molecules."""
    total_charge = 0
    for molecule in molecules:
        for atom in molecule.GetAtoms():
            total_charge += atom.GetFormalCharge()
    return total_charge

def get_reaction_charges(reaction):
    """Return the total charge on the reactant and product side of a reaction."""
    # Parse the reaction

    # Calculate total charge for reactants and products
    reactant_charge = calculate_total_charge(reaction.GetReactants())
    product_charge = calculate_total_charge(reaction.GetProducts())

    return reactant_charge, product_charge

def count_elements(molecule):
    """Count elements in a molecule."""
    element_count = defaultdict(int)
    for atom in molecule.GetAtoms():
        element_count[atom.GetSymbol()] += 1
    return element_count

def is_reaction_balanced_count(file_path):
    """Check if a chemical reaction is balanced."""
    # Parse the reaction
    reaction = rdChemReactions.ReactionFromRxnFile(file_path)

    # Count elements in reactants and products
    reactant_elements = defaultdict(int)
    product_elements = defaultdict(int)

    for reactant in reaction.GetReactants():
        reactant.UpdatePropertyCache(strict=False)
        reactant = Chem.AddHs(reactant)
        for count in count_elements(reactant).items():
            reactant_elements[count[0]] += count[1]

    for product in reaction.GetProducts():
        product.UpdatePropertyCache(strict=False)
        product = Chem.AddHs(product)
        for count in count_elements(product).items():
            product_elements[count[0]] += count[1]
    return reactant_elements == product_elements, (reactant_elements,product_elements)

def is_reaction_balanced_charge(file_path):
    """Check if a chemical reaction is balanced."""
    # Parse the reaction
    reaction = rdChemReactions.ReactionFromRxnFile(file_path)

    # Calculate total charge for reactants and products
    reactant_charge, product_charge = get_reaction_charges(reaction)

    # Compare total charge
    return reactant_charge == product_charge, (reactant_charge, product_charge)

def get_molecular_formula(file_path,direction):
    """Construct the molecular formula of a reaction."""
    # Parse the reaction
    reaction = rdChemReactions.ReactionFromRxnFile(file_path)

    # Initialize a list to store molecular formulas
    reactant_formulas = []
    product_formulas = []
    symb_to_name = {}
    # Loop through reactants and products
    for reactant in reaction.GetReactants():
        try:
            Chem.SanitizeMol(reactant)  # Ensure molecule is properly sanitized
        except:
            pass
        try:
            Chem.AssignStereochemistry(reactant)  # Assign stereochemistry if needed
        except:
            pass
        try:
            reactant.UpdatePropertyCache()
        except:
            pass
        reactant_formulas.append(CalcMolFormula(reactant))
        symb_to_name.update({atom.GetSymbol(): PeriodicTable.GetElementName(GetPeriodicTable(),atom.GetAtomicNum()) for atom in reactant.GetAtoms()})


    for product in reaction.GetProducts():
        try:
            Chem.SanitizeMol(product)  # Ensure molecule is properly sanitized
        except:
            pass
        try:
            Chem.AssignStereochemistry(product)  # Assign stereochemistry if needed
        except:
            pass
        try:
            product.UpdatePropertyCache()
        except:
            pass
        # Chem.GetImplicitHs(product)  # Calculate implicit hydrogens
        product_formulas.append(CalcMolFormula(product))
        symb_to_name.update({atom.GetSymbol(): PeriodicTable.GetElementName(GetPeriodicTable(),atom.GetAtomicNum()) for atom in product.GetAtoms()})

    # Construct reaction formula
    reactant_side = ' + '.join(reactant_formulas)
    product_side = ' + '.join(product_formulas)
    reaction_formula = f"{reactant_side} -> {product_side}" if direction == 'forward' else f"{reactant_side} <=> {reactant_side}"

    return reaction_formula,symb_to_name

def get_reaction_info(rxn_file_path,direction):

    balanced_count,(subs_atoms,prods_atoms) = is_reaction_balanced_count(rxn_file_path)
    balanced_charge,(subs_charge, prods_charge) = is_reaction_balanced_charge(rxn_file_path)
    molc_formula,symb_to_name = get_molecular_formula(rxn_file_path,direction)
    return balanced_count,(subs_atoms,prods_atoms),balanced_charge, (subs_charge, prods_charge),molc_formula,symb_to_name

def construct_vmh_formula(reaction, subs_abbr, prods_abbr):
    """Construct the VMH formula of a reaction."""
    subs_stch = [float(x) for x in json.loads(reaction.subs_sch)]
    prods_stch = [float(x) for x in json.loads(reaction.prods_sch)]
    subs_comps = json.loads(reaction.subs_comps)
    prods_comps = json.loads(reaction.prods_comps)
    substrate_formula = ' + '.join([f"{subs_stch[i]} {subs_abbr[i]}[{subs_comps[i]}]" for i in range(len(subs_stch))])
    product_formula = ' + '.join([f"{prods_stch[i]} {prods_abbr[i]}[{prods_comps[i]}]" for i in range(len(prods_stch))])
    direction = '->' if reaction.direction == 'forward' else '<=>'
    formula = f"{substrate_formula} {direction} {product_formula}"
    return formula