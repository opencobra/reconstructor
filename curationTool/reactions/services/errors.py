"""Typed errors raised by the service layer and mapped to API error envelopes.

Mapping to HTTP status codes lives in reactions/api/__init__.py; the services
themselves stay transport-agnostic.
"""


class ServiceError(Exception):
    """Base class: machine-readable code + human message (+ optional detail)."""

    code = 'processing_error'
    http_status = 500

    def __init__(self, message, detail=None):
        super().__init__(message)
        self.message = message
        self.detail = detail or {}


class ServiceValidationError(ServiceError):
    code = 'validation_error'
    http_status = 422


class UnsupportedTypeError(ServiceError):
    code = 'unsupported_type'
    http_status = 422


class UnknownIdentifierError(ServiceError):
    code = 'unknown_identifier'
    http_status = 422


class DuplicateNameError(ServiceError):
    code = 'duplicate_name'
    http_status = 409


class NotFoundError(ServiceError):
    code = 'not_found'
    http_status = 404


class UpstreamServiceError(ServiceError):
    code = 'upstream_unavailable'
    http_status = 502


class ProcessingError(ServiceError):
    code = 'processing_error'
    http_status = 500
