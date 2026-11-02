"""Pydantic schemas for API v1 (contracts/api-v1.md)."""
from typing import List, Optional

from ninja import Schema
from pydantic import field_validator, model_validator

from reactions.services.constants import (  # noqa: F401 (re-exported)
    API_TO_INTERNAL_TYPE,
    COMPARTMENTS,
    MAX_STOICHIOMETRY,
    SUPPORTED_TYPES,
)


class ErrorBody(Schema):
    code: str
    message: str
    detail: Optional[dict] = None


class ErrorEnvelope(Schema):
    error: ErrorBody


class MetaboliteInput(Schema):
    identifier: str
    type: str
    stoichiometry: int = 1
    compartment: str = 'c'
    name: Optional[str] = None

    @field_validator('identifier')
    @classmethod
    def check_identifier(cls, v):
        if not v or not v.strip():
            raise ValueError('identifier must not be empty')
        return v.strip()

    @field_validator('type')
    @classmethod
    def check_type(cls, v):
        if v not in API_TO_INTERNAL_TYPE:
            raise ValueError(
                f"unsupported type '{v}'; supported types: {', '.join(SUPPORTED_TYPES)}")
        return v

    @field_validator('stoichiometry')
    @classmethod
    def check_stoichiometry(cls, v):
        if v < 1 or v > MAX_STOICHIOMETRY:
            raise ValueError(
                f'stoichiometry must be a whole number between 1 and {MAX_STOICHIOMETRY}')
        return v

    @field_validator('compartment')
    @classmethod
    def check_compartment(cls, v):
        if v not in COMPARTMENTS:
            raise ValueError(
                f"unknown compartment '{v}'; accepted compartments: {', '.join(COMPARTMENTS)}")
        return v

    @model_validator(mode='after')
    def check_name_for_smiles(self):
        if self.type == 'smiles' and not (self.name and self.name.strip()):
            raise ValueError("'name' is required for metabolites of type 'smiles'")
        return self


class ResolveIn(Schema):
    identifier: str
    type: str

    @field_validator('identifier')
    @classmethod
    def check_identifier(cls, v):
        if not v or not v.strip():
            raise ValueError('identifier must not be empty')
        return v.strip()

    @field_validator('type')
    @classmethod
    def check_type(cls, v):
        if v not in API_TO_INTERNAL_TYPE:
            raise ValueError(
                f"unsupported type '{v}'; supported types: {', '.join(SUPPORTED_TYPES)}")
        return v


class VmhPresence(Schema):
    found: bool
    abbreviation: Optional[str] = None
    miriam: Optional[str] = None


class MetaboliteResolution(Schema):
    name: Optional[str] = None
    formula: Optional[str] = None
    charge: Optional[int] = None
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    mol_weight: Optional[float] = None
    vmh: VmhPresence


class ReactionSubmission(Schema):
    substrates: List[MetaboliteInput] = []
    products: List[MetaboliteInput] = []
    direction: str = 'forward'
    subsystem: Optional[str] = None
    organs: List[str] = []

    @field_validator('direction')
    @classmethod
    def check_direction(cls, v):
        if v not in ('forward', 'reversible'):
            raise ValueError("direction must be 'forward' or 'reversible'")
        return v

    @model_validator(mode='after')
    def check_sides(self):
        if not self.substrates and not self.products:
            raise ValueError(
                'enter at least one substrate or one product')
        return self


class VmhMatch(Schema):
    found: bool
    similar: bool
    url: Optional[object] = None
    formula: Optional[object] = None


class ProcessedReaction(Schema):
    reaction_id: int
    direction: str
    subsystem: Optional[str] = None
    molc_formula: Optional[object] = None
    balanced_count: Optional[object] = None
    balanced_charge: Optional[object] = None
    subs_atoms: Optional[object] = None
    prods_atoms: Optional[object] = None
    subs_charge: Optional[object] = None
    prods_charge: Optional[object] = None
    symb_to_name: Optional[object] = None
    metabolite_names: List[str] = []
    metabolite_formulas: List[str] = []
    metabolite_charges: List[object] = []
    metabolite_smiles: List[Optional[str]] = []
    metabolite_inchi_keys: List[Optional[str]] = []
    metabolite_mol_weights: List[Optional[float]] = []
    subs_found: List[bool] = []
    subs_miriams: List[Optional[object]] = []
    prod_found: List[bool] = []
    prod_miriams: List[Optional[object]] = []
    vmh: VmhMatch


class SaveReactionIn(Schema):
    name: str
    description: str = ''

    @field_validator('name')
    @classmethod
    def check_name(cls, v):
        if not v or not v.strip():
            raise ValueError('name must not be empty')
        return v.strip()


class SavedReactionOut(Schema):
    reaction_id: int
    name: str
    saved: bool = True


class ReactionListItem(Schema):
    reaction_id: int
    name: Optional[str] = None
    rxn_formula: Optional[object] = None
    balanced_mass: Optional[bool] = None
    balanced_charge: Optional[bool] = None
    vmh_found: bool = False


class ReactionList(Schema):
    count: int
    page: int
    page_size: int
    results: List[ReactionListItem]


class VmhFetchIn(Schema):
    abbreviation: str


class RheaFetchIn(Schema):
    rhea_id: str


class MeOut(Schema):
    user_id: int
    name: Optional[str] = None
    full_name: Optional[str] = None
    key_prefix: str
    key_created_at: Optional[str] = None
