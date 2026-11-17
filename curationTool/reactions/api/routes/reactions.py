"""/api/v1/reactions endpoints."""
import json

from ninja import Router

from reactions.api.schemas import (
    ProcessedReaction,
    ReactionList,
    ReactionSubmission,
    SavedReactionOut,
    SaveReactionIn,
)
from reactions.services.reaction_processing import (
    get_user_reaction,
    process_reaction,
    save_reaction_for_user,
    serialize_reaction,
)
from reactions.utils.utils import safe_json_loads

router = Router(tags=['reactions'])


@router.post('', response={201: ProcessedReaction},
             summary='Create and balance-check a reaction')
def create(request, payload: ReactionSubmission):
    """Process a reaction (resolve metabolites, mass/charge balance, VMH match)
    and persist it. Save it to your collection with POST /reactions/{id}/save."""
    _, data = process_reaction(
        substrates=[m.dict() for m in payload.substrates],
        products=[m.dict() for m in payload.products],
        direction=payload.direction,
        subsystem=payload.subsystem,
        organs=payload.organs,
        user=request.auth,
    )
    return 201, data


@router.post('/{int:reaction_id}/save', response=SavedReactionOut,
             summary="Save a reaction to the caller's collection")
def save(request, reaction_id: int, payload: SaveReactionIn):
    reaction = save_reaction_for_user(
        request.auth, reaction_id, payload.name, payload.description)
    return SavedReactionOut(reaction_id=reaction.id, name=reaction.short_name)


@router.get('', response=ReactionList, summary="List the caller's saved reactions")
def list_reactions(request, page: int = 1, page_size: int = 25):
    page = max(1, page)
    page_size = max(1, min(page_size, 100))
    queryset = request.auth.saved_reactions.order_by('-id')
    count = queryset.count()
    start = (page - 1) * page_size
    results = []
    for reaction in queryset[start:start + page_size]:
        def first(field):
            values = safe_json_loads(field)
            return values[0] if isinstance(values, list) and values else values
        results.append({
            'reaction_id': reaction.id,
            'name': reaction.short_name,
            'rxn_formula': (safe_json_loads(reaction.rxn_formula)
                            if reaction.rxn_formula else first(reaction.molc_formula)),
            'balanced_mass': bool(first(reaction.balanced_count)),
            'balanced_charge': bool(first(reaction.balanced_charge)),
            'vmh_found': reaction.vmh_found,
        })
    return {'count': count, 'page': page, 'page_size': page_size,
            'results': results}


@router.get('/{int:reaction_id}', response=ProcessedReaction,
            summary='Get one of your reactions')
def get(request, reaction_id: int):
    reaction = get_user_reaction(request.auth, reaction_id)
    return serialize_reaction(reaction)
