"""Shared constants for the v1 programmatic surface."""

# Compartments offered by the web UI (ReactantsinputFieldManagement.js).
COMPARTMENTS = ['c', 'e', 'g', 'm', 'l', 'n', 'r', 'x', '-']

# Upper bound on stoichiometric coefficients. The rxn-file builder
# (construct_reaction_rxnfile) duplicates each molecule once per unit, so an
# unbounded coefficient lets a single request exhaust a worker's memory.
MAX_STOICHIOMETRY = 100

# API type values -> internal converter types (reactions/utils/to_mol.py).
# File-based types ('MDL Mol file', 'Draw') are intentionally not exposed.
API_TO_INTERNAL_TYPE = {
    'vmh': 'VMH',
    'chebi': 'ChEBI ID',
    'chebi_name': 'ChEBI Name',
    'pubchem': 'PubChem ID',
    'swisslipids': 'SwissLipids',
    'smiles': 'SMILES',
    'saved': 'Saved',
}
SUPPORTED_TYPES = sorted(API_TO_INTERNAL_TYPE)
