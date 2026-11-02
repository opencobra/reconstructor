"""/api/v1/metabolites endpoints."""
from ninja import Router

from reactions.api.schemas import MetaboliteResolution, ResolveIn
from reactions.services.metabolite_resolution import resolve_metabolite

router = Router(tags=['metabolites'])


@router.post('/resolve', response=MetaboliteResolution,
             summary='Resolve a metabolite identifier')
def resolve(request, payload: ResolveIn):
    """Resolve a VMH abbreviation, ChEBI ID/name, SwissLipids ID, PubChem CID,
    SMILES string or saved-metabolite id to name, formula, charge, structure
    and VMH presence."""
    return resolve_metabolite(payload.identifier, payload.type, user=request.auth)
