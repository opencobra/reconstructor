"""/api/v1/sources endpoints: import reactions from VMH / Rhea."""
from ninja import Router

from reactions.api.schemas import ProcessedReaction, RheaFetchIn, VmhFetchIn
from reactions.services.external_sources import fetch_rhea_reaction, fetch_vmh_reaction

router = Router(tags=['sources'])


@router.post('/vmh/fetch', response={201: ProcessedReaction},
             summary='Import a reaction from VMH by abbreviation')
def vmh_fetch(request, payload: VmhFetchIn):
    _, data = fetch_vmh_reaction(payload.abbreviation.strip(), request.auth)
    return 201, data


@router.post('/rhea/fetch', response={201: ProcessedReaction},
             summary='Import a reaction from Rhea by id')
def rhea_fetch(request, payload: RheaFetchIn):
    _, data = fetch_rhea_reaction(payload.rhea_id, request.auth)
    return 201, data
