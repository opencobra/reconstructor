"""Generation, rotation and authentication of per-user API keys (FR-002/FR-003)."""
import hashlib
import secrets

from django.utils import timezone

from reactions.models import ApiKey

KEY_PREFIX = 'ctor_'
PREFIX_DISPLAY_LEN = 12


def _digest(raw_key):
    return hashlib.sha256(raw_key.encode('utf-8')).hexdigest()


def generate_key(user):
    """Create (or rotate) the user's API key. Returns the raw key — the only
    time it is ever available."""
    raw_key = KEY_PREFIX + secrets.token_urlsafe(32)
    ApiKey.objects.update_or_create(
        user=user,
        defaults={
            'key_hash': _digest(raw_key),
            'prefix': raw_key[:PREFIX_DISPLAY_LEN],
            'last_used_at': None,
        },
    )
    return raw_key


def authenticate(raw_key):
    """Return the owning User for a valid key, else None. Lookup is by SHA-256
    digest (indexed), so no secret comparison happens in Python."""
    if not raw_key or not raw_key.startswith(KEY_PREFIX):
        return None
    try:
        api_key = ApiKey.objects.select_related('user').get(
            key_hash=_digest(raw_key))
    except ApiKey.DoesNotExist:
        return None
    ApiKey.objects.filter(pk=api_key.pk).update(last_used_at=timezone.now())
    return api_key.user
