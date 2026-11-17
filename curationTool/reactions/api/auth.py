"""API-key authentication for /api/v1/ (Bearer or X-API-Key header)."""
from ninja.security import APIKeyHeader, HttpBearer

from reactions.services import api_keys


class BearerKeyAuth(HttpBearer):
    def authenticate(self, request, token):
        return api_keys.authenticate(token)


class HeaderKeyAuth(APIKeyHeader):
    param_name = 'X-API-Key'

    def authenticate(self, request, key):
        return api_keys.authenticate(key)


API_AUTH = [BearerKeyAuth(), HeaderKeyAuth()]
