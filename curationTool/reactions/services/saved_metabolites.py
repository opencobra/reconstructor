"""Access control for saved metabolites referenced by id (API type 'saved')."""
from django.db.models import Q

from reactions.models import SavedMetabolite


def visible_to_user(identifier, user):
    """True iff ``identifier`` is the pk of a SavedMetabolite the user owns or
    was shared. Foreign and nonexistent ids are indistinguishable to callers,
    so rejecting on this check leaks no existence information."""
    identifier = str(identifier).strip()
    if user is None or not identifier.isdigit():
        return False
    return SavedMetabolite.objects.filter(
        Q(owner=user) | Q(shared_with=user), pk=int(identifier)).exists()
