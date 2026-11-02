"""Constructor API v1 (feature 002-curation-api-mcp).

Mounted at /api/v1/ from reactions_project/urls.py. OpenAPI: /api/v1/openapi.json,
interactive docs: /api/v1/docs.
"""
from django.conf import settings
from ninja import NinjaAPI
from ninja.errors import AuthenticationError, Throttled
from ninja.errors import ValidationError as NinjaValidationError
from ninja.throttling import AuthRateThrottle

from reactions.api.auth import API_AUTH
from reactions.api.schemas import MeOut
from reactions.services.errors import ServiceError

api_throttle = AuthRateThrottle(getattr(settings, 'API_THROTTLE_RATE', '60/m'))

api = NinjaAPI(
    title='Constructor API',
    version='1.0.0',
    auth=API_AUTH,
    throttle=[api_throttle],
    urls_namespace='api_v1',
)


def _error(request, status, code, message, detail=None):
    body = {'error': {'code': code, 'message': message}}
    if detail:
        body['error']['detail'] = detail
    return api.create_response(request, body, status=status)


@api.exception_handler(ServiceError)
def service_error(request, exc):
    return _error(request, exc.http_status, exc.code, exc.message,
                  exc.detail or None)


@api.exception_handler(AuthenticationError)
def auth_error(request, exc):
    return _error(request, 401, 'unauthorized',
                  'Missing or invalid API key. Pass it as '
                  '"Authorization: Bearer <key>" or "X-API-Key: <key>".')


@api.exception_handler(Throttled)
def throttled(request, exc):
    return _error(request, 429, 'rate_limited',
                  'Rate limit exceeded; slow down and retry shortly.',
                  {'retry_after_seconds': exc.wait})


@api.exception_handler(NinjaValidationError)
def validation_error(request, exc):
    first = exc.errors[0] if exc.errors else {}
    loc = '.'.join(str(part) for part in first.get('loc', []) if part != 'body')
    msg = first.get('msg', 'Invalid input.')
    message = f'{loc}: {msg}' if loc else msg
    return _error(request, 422, 'validation_error', message,
                  {'errors': exc.errors})


@api.get('/health', auth=None, summary='Liveness check (no auth)')
def health(request):
    return {'status': 'ok', 'version': 'v1'}


@api.get('/me', response=MeOut, summary='Identify the caller')
def me(request):
    user = request.auth
    key = user.api_key
    return MeOut(
        user_id=user.pk,
        name=user.name,
        full_name=user.full_name,
        key_prefix=key.prefix,
        key_created_at=key.created_at.isoformat() if key.created_at else None,
    )


from reactions.api.routes.metabolites import router as metabolites_router  # noqa: E402
from reactions.api.routes.reactions import router as reactions_router  # noqa: E402
from reactions.api.routes.sources import router as sources_router  # noqa: E402

api.add_router('/metabolites', metabolites_router)
api.add_router('/reactions', reactions_router)
api.add_router('/sources', sources_router)
