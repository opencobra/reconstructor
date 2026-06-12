"""
JSON endpoints for the multi-reaction tabbed workspace (feature
001-multi-reaction-tabbed).

Each tab is a `ReactionGroup` (owner FK, name, is_active) holding an ordered set
of `Reaction` rows via the `ReactionGroupMembership` through-model. These views
mirror the request/response shape of the existing `send_to_workspace` /
`remove_from_workspace` endpoints (`vmh_views.py`) but target groups instead of
the single per-user `Workspace`.

Invariants enforced here:
- Always at least one group per user (FR-008): `get_default_group` seeds one.
- Deleting a group never deletes member reactions (FR-009).
- Removing a reaction from a group never deletes the saved reaction (FR-011).
- A reaction appears at most once per group (FR-012): duplicate-guarded on add.
"""
import ast
import json

from django.db import models
from django.http import JsonResponse
from django.shortcuts import get_object_or_404
from django.views.decorators.http import require_GET, require_POST

from reactions.models import (
    Reaction,
    ReactionGroup,
    ReactionGroupMembership,
    User,
)

DEFAULT_GROUP_NAME = 'My Reactions'


def get_default_group(user):
    """
    Return the user's active group, seeding a default one if none exist.

    Guarantees the FR-008 invariant (at least one group/tab per user) and that
    exactly one group is marked active. Used on first main-page load (FR-008,
    SC-005) and as a fallback target for add/remove.
    """
    active = user.reaction_groups.filter(is_active=True).first()
    if active:
        return active
    existing = user.reaction_groups.order_by('created_at').first()
    if existing:
        existing.is_active = True
        existing.save(update_fields=['is_active'])
        return existing
    return user.reaction_groups.create(
        name=DEFAULT_GROUP_NAME, is_active=True)


def _set_active(user, group):
    """Mark `group` as the user's single active group."""
    user.reaction_groups.exclude(pk=group.pk).filter(
        is_active=True).update(is_active=False)
    if not group.is_active:
        group.is_active = True
        group.save(update_fields=['is_active'])


def _clean_formula(value):
    """
    Several Reaction fields persist a stringified single-element list, e.g.
    `["atp[c] -> adp[c]"]`. Surface the human-readable formula without altering
    reaction semantics (Constitution IV) — return the inner string, trimmed.
    """
    if not value:
        return ''
    text = value.strip() if isinstance(value, str) else value
    if isinstance(text, str) and text.startswith('[') and text.endswith(']'):
        try:
            parsed = ast.literal_eval(text)
            if isinstance(parsed, (list, tuple)):
                return ' '.join(str(p) for p in parsed).strip()
            return str(parsed).strip()
        except (ValueError, SyntaxError):
            return text
    return str(text).strip()


def _reaction_card(reaction):
    """Compact summary used to render a reaction card in the group view."""
    return {
        'id': reaction.pk,
        'short_name': reaction.short_name or '',
        'rxn_formula': _clean_formula(reaction.rxn_formula),
        'molc_formula': _clean_formula(reaction.molc_formula),
        'substrates': reaction.substrates or '',
        'products': reaction.products or '',
        'mass_balanced': bool(reaction.balanced_count),
        'charge_balanced': bool(reaction.balanced_charge),
        'subsystem': reaction.subsystem or '',
        'direction': reaction.direction or '',
        'flags': [
            {
                'id': flag.pk,
                'name': flag.name_flag or '',
                'color': flag.color or '',
            }
            for flag in reaction.flags.all()
        ],
    }


def _group_summary(group):
    return {
        'id': group.pk,
        'name': group.name,
        'is_active': group.is_active,
        'count': group.memberships.count(),
    }


def _require_user(request, data=None):
    """Resolve the acting user from POST body or session; None if anonymous."""
    user_id = None
    if data is not None:
        user_id = data.get('userID') or data.get('user_id')
    if not user_id:
        user_id = request.GET.get('userID') or request.GET.get('user_id')
    if not user_id:
        user_id = request.session.get('userID')
    if not user_id:
        return None
    return get_object_or_404(User, pk=user_id)


@require_GET
def list_groups(request):
    """Return all of the user's groups plus the active group id (FR-005, FR-007)."""
    user = _require_user(request)
    if user is None:
        return JsonResponse({'status': 'anonymous', 'groups': []})
    get_default_group(user)  # ensure ≥1 group exists (FR-008)
    groups = [_group_summary(g) for g in user.reaction_groups.all()]
    active = next((g['id'] for g in groups if g['is_active']), None)
    return JsonResponse(
        {'status': 'success', 'groups': groups, 'active_group_id': active})


@require_GET
def group_contents(request):
    """
    Return one group's ordered members as reaction cards (FR-001).

    Query params: `groupId` (optional; defaults to the user's active group),
    `userID` (optional; falls back to session).
    """
    user = _require_user(request)
    if user is None:
        return JsonResponse({'status': 'anonymous', 'reactions': []})
    group_id = request.GET.get('groupId')
    if group_id:
        group = get_object_or_404(ReactionGroup, pk=group_id, owner=user)
        _set_active(user, group)
    else:
        group = get_default_group(user)
    cards = [
        _reaction_card(m.reaction)
        for m in group.memberships.select_related('reaction').prefetch_related(
            'reaction__flags')
    ]
    return JsonResponse({
        'status': 'success',
        'group': _group_summary(group),
        'reactions': cards,
    })


@require_POST
def create_group(request):
    """Create a new group and make it active (FR-009)."""
    data = json.loads(request.body)
    user = _require_user(request, data)
    if user is None:
        return JsonResponse(
            {'status': 'error', 'message': 'Login required'}, status=403)
    name = (data.get('name') or '').strip() or DEFAULT_GROUP_NAME
    group = user.reaction_groups.create(name=name)
    _set_active(user, group)
    return JsonResponse({'status': 'success', 'group': _group_summary(group)})


@require_POST
def rename_group(request):
    """Rename a group (FR-009)."""
    data = json.loads(request.body)
    user = _require_user(request, data)
    if user is None:
        return JsonResponse(
            {'status': 'error', 'message': 'Login required'}, status=403)
    group = get_object_or_404(
        ReactionGroup, pk=data.get('groupId'), owner=user)
    name = (data.get('name') or '').strip()
    if not name:
        return JsonResponse(
            {'status': 'error', 'message': 'Name required'}, status=400)
    group.name = name
    group.save(update_fields=['name', 'updated_at'])
    return JsonResponse({'status': 'success', 'group': _group_summary(group)})


@require_POST
def delete_group(request):
    """
    Delete a group. Member reactions are NOT deleted (FR-009): the M2M through
    rows cascade, but `Reaction` rows persist. Always leaves ≥1 group (FR-008).
    """
    data = json.loads(request.body)
    user = _require_user(request, data)
    if user is None:
        return JsonResponse(
            {'status': 'error', 'message': 'Login required'}, status=403)
    group = get_object_or_404(
        ReactionGroup, pk=data.get('groupId'), owner=user)
    was_active = group.is_active
    group.delete()
    # Re-seed / re-activate so the workspace is never tab-less (FR-008).
    fallback = get_default_group(user)
    if was_active:
        _set_active(user, fallback)
    return JsonResponse(
        {'status': 'success', 'active_group_id': fallback.pk})


@require_GET
def saved_reactions_for_group(request):
    """
    Return the user's saved reactions as picker cards for the add-from-saved
    flow (FR-010). Optionally flags which are already in a group so the UI can
    pre-empt duplicates (FR-012). Query params: `userID`, `groupId` (optional).
    """
    user = _require_user(request)
    if user is None:
        return JsonResponse({'status': 'anonymous', 'reactions': []})
    in_group = set()
    group_id = request.GET.get('groupId')
    if group_id:
        group = get_object_or_404(ReactionGroup, pk=group_id, owner=user)
        in_group = set(group.memberships.values_list('reaction_id', flat=True))
    cards = []
    for reaction in user.saved_reactions.prefetch_related('flags').all():
        card = _reaction_card(reaction)
        card['in_group'] = reaction.pk in in_group
        cards.append(card)
    return JsonResponse({'status': 'success', 'reactions': cards})


@require_POST
def add_to_group(request):
    """
    Add saved reactions into a group (FR-010), duplicate-guarded (FR-012).

    Mirrors `send_to_workspace`: body carries `userID`, `groupId` (optional,
    defaults to active), and `reactionIds` (list).
    """
    data = json.loads(request.body)
    user = _require_user(request, data)
    if user is None:
        return JsonResponse(
            {'status': 'error', 'message': 'Login required'}, status=403)
    group_id = data.get('groupId')
    group = (get_object_or_404(ReactionGroup, pk=group_id, owner=user)
             if group_id else get_default_group(user))
    reaction_ids = data.get('reactionIds')
    if reaction_ids is None and data.get('reactionId'):
        reaction_ids = [data.get('reactionId')]
    reaction_ids = reaction_ids or []

    next_pos = (group.memberships.aggregate(
        m=models.Max('position'))['m'] or 0)
    added, skipped = [], []
    for rid in reaction_ids:
        reaction = get_object_or_404(Reaction, pk=rid)
        _, created = ReactionGroupMembership.objects.get_or_create(
            group=group, reaction=reaction,
            defaults={'position': next_pos + 1})
        if created:
            next_pos += 1
            added.append(reaction.pk)
        else:
            skipped.append(reaction.pk)  # already a member (FR-012)
    return JsonResponse({
        'status': 'success', 'added': added, 'already_present': skipped})


@require_POST
def remove_from_group(request):
    """
    Remove a reaction from a group (FR-011). The saved `Reaction` is untouched;
    only the membership row is deleted.
    """
    data = json.loads(request.body)
    user = _require_user(request, data)
    if user is None:
        return JsonResponse(
            {'status': 'error', 'message': 'Login required'}, status=403)
    group_id = data.get('groupId')
    group = (get_object_or_404(ReactionGroup, pk=group_id, owner=user)
             if group_id else get_default_group(user))
    reaction_id = data.get('reactionId')
    ReactionGroupMembership.objects.filter(
        group=group, reaction_id=reaction_id).delete()
    return JsonResponse({'status': 'success'})


@require_POST
def discard_clone(request):
    """
    Discard an unsaved clone (US4 / FR-015): fully delete a clone reaction the
    user just created so nothing is persisted. Guards on ownership — only a
    reaction in the user's saved set is removed — then drops its group
    memberships and deletes the row. Never touches any other reaction.
    """
    data = json.loads(request.body)
    user = _require_user(request, data)
    if user is None:
        return JsonResponse(
            {'status': 'error', 'message': 'Login required'}, status=403)
    reaction_id = data.get('reactionId')
    reaction = get_object_or_404(Reaction, pk=reaction_id)
    if not user.saved_reactions.filter(pk=reaction.pk).exists():
        # Refuse to delete anything the user does not own (safety).
        return JsonResponse(
            {'status': 'error', 'message': 'Not your reaction'}, status=403)
    user.saved_reactions.remove(reaction)
    ReactionGroupMembership.objects.filter(reaction=reaction).delete()
    reaction.delete()
    return JsonResponse({'status': 'success'})
