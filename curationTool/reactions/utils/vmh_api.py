"""Client for the VMH public API.

VMH (https://www.vmh.life) retired the legacy ``/_api/...`` endpoints and the
interim ``vmh.chatimd.org`` host; every read now goes through
``https://www.vmh.life/api/public`` (docs: /apiDocs, migration guide:
/apiMigration).  This module is the only place that talks to that API.

What callers can rely on:

* **"Not found" and "could not ask" are different outcomes.**  The strict
  helpers (``lookup_*``, ``*_exists``, :func:`match_metabolite_structure`,
  :func:`search_*`) return ``None``/``False``/``[]`` only when VMH answered
  that there is no such entry; when VMH could not be asked (network error,
  timeout, HTTP 429/5xx, unexpected payload) they raise
  :class:`VMHUnavailableError`.  The older ``find_*``/``gene_*`` helpers keep
  their historical never-raise contract for read paths where an unanswered
  question may be treated as "no"; they log a warning instead.
* **Rate limit.**  Reads need no key, but VMH then allows 100 requests per hour
  per client IP, shared by all endpoints.  ``VMH_API_KEY`` (a VMH
  ``public_unlimited`` key, sent as ``x-api-key``) lifts the limit.  After an
  HTTP 429 every call fails fast until VMH's advertised reset time.
* **Caching.**  Answers, including "not found", are cached in the ``vmh``
  cache (``VMH_CACHE_TTL`` / ``VMH_NEGATIVE_CACHE_TTL`` seconds); failures are
  never cached.
* **Structure lookups** (InChIKey / SMILES) use the SPARQL endpoint (~1 s).
  The ``advancedSearch`` regex filters scan the whole table (~40 s) and are
  only used as a fallback.  About half of VMH's metabolites carry a SMILES but
  no InChIKey; with a VMH API key their InChIKeys are derived locally with
  RDKit into an index that is cached and rebuilt daily
  (:func:`_derived_inchikey_index`).
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import re
import tempfile
import threading
import time
from contextlib import contextmanager
from typing import Any, Callable, Dict, Iterator, List, Optional, Sequence
from urllib.parse import quote

import requests
from django.conf import settings
from django.core.cache import caches
from django.core.cache.backends.base import InvalidCacheBackendError

try:  # POSIX only; without it the derived index is built without a cross-process lock
    import fcntl
except ImportError:  # pragma: no cover - Windows
    fcntl = None

logger = logging.getLogger(__name__)

VMH_SITE_URL = "https://www.vmh.life"
VMH_API_URL = VMH_SITE_URL + "/api/public"

DEFAULT_TIMEOUT = 20        # seconds; the handler endpoints answer in well under 1 s
SPARQL_TIMEOUT = 30         # literal-pattern SPARQL queries take ~1-2 s (more when VMH is cold)
SLOW_SEARCH_TIMEOUT = 90    # advancedSearch regex filters scan the whole table (~40 s)
BULK_TIMEOUT = 180          # whole-table SPARQL export used for the derived index (~25 s)

DEFAULT_CACHE_TTL = 24 * 3600
DEFAULT_NEGATIVE_CACHE_TTL = 3600
INDEX_MAX_AGE = 24 * 3600           # rebuild the derived InChIKey index daily...
INDEX_CACHE_TTL = 7 * 24 * 3600     # ...but keep serving a stale one meanwhile
INDEX_RETRY_AFTER_FAILURE = 600     # do not retry a failed index build for 10 minutes

_CACHE_VERSION = 1
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_ENTREZ_RE = re.compile(r"^\d+(?:\.\d+)?$")
_ARROW = r"(?:->|<=>)"
# Characters with a meaning in JavaScript regular expressions (incl. the 'u' flag).
_JS_REGEX_SPECIAL = re.compile(r"([\\^$.|?*+()\[\]{}/])")


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

def _setting(name: str, default: Any = None) -> Any:
    """A Django setting (so tests can override it), else the environment."""
    if hasattr(settings, name):
        value = getattr(settings, name)
    else:
        value = os.getenv(name)
    return default if value in (None, "") else value


def _int_setting(name: str, default: int) -> int:
    try:
        return int(_setting(name, default))
    except (TypeError, ValueError):
        return default


def _bool_setting(name: str, default: bool) -> bool:
    value = _setting(name, default)
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


# Hosts VMH no longer serves from (vmh.chatimd.org was the staging host of the
# revamped site; it now presents another service's TLS certificate).
_RETIRED_HOSTS = ("vmh.chatimd.org",)
_retired_host_warned: set = set()


def _url_setting(name: str, default: str) -> str:
    value = str(_setting(name, default)).strip().rstrip("/")
    if any(host in value.lower() for host in _RETIRED_HOSTS):
        if name not in _retired_host_warned:
            _retired_host_warned.add(name)
            logger.warning("%s=%s points at a retired VMH host; using %s instead. "
                           "Update the setting.", name, value, default)
        return default
    return value


def vmh_api_base_url() -> str:
    """Base URL of the VMH public API (``VMH_API_BASE_URL`` overrides it)."""
    return _url_setting("VMH_API_BASE_URL", VMH_API_URL)


def _vmh_api_key() -> str:
    return str(_setting("VMH_API_KEY", "") or "").strip()


def vmh_api_key_configured() -> bool:
    """True when a VMH API key is configured (it lifts VMH's anonymous rate limit)."""
    return bool(_vmh_api_key())


def vmh_public_base_url() -> str:
    """Prefix used to build clickable VMH links such as /metabolite/{abbr}."""
    return _url_setting("VMH_BASE_URL", VMH_SITE_URL)


def vmh_metabolite_url(abbr: str) -> str:
    return f"{vmh_public_base_url()}/metabolite/{quote(str(abbr), safe='')}"


def vmh_reaction_url(abbr: str) -> str:
    return f"{vmh_public_base_url()}/reaction/{quote(str(abbr), safe='')}"


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class VMHUnavailableError(Exception):
    """VMH could not answer: network error, timeout, rate limit, server error or
    an unexpected response.  Never raised for "no such entry"."""

    def __init__(self, reason: str, *, retry_after: Optional[int] = None,
                 status: Optional[int] = None):
        super().__init__(reason)
        self.reason = reason
        self.retry_after = retry_after
        self.status = status

    @property
    def rate_limited(self) -> bool:
        return self.status == 429

    def user_message(self, action: str = "reach VMH") -> str:
        """A complete sentence suitable for showing to a curator or API client."""
        message = f"Could not {action}: {self.reason}."
        if self.retry_after:
            minutes = max(1, math.ceil(self.retry_after / 60))
            message += (f" Please try again in about {minutes} "
                        f"minute{'s' if minutes != 1 else ''}.")
        else:
            message += " Please try again shortly."
        if self.rate_limited and not vmh_api_key_configured():
            message += (" (Administrators: set VMH_API_KEY to a VMH 'public_unlimited'"
                        " API key to lift VMH's limit of 100 requests per hour.)")
        return message

    def __str__(self) -> str:
        return self.user_message()


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

_MISS = object()


def _cache():
    try:
        return caches["vmh"]
    except InvalidCacheBackendError:
        return caches["default"]


def _cache_key(*parts: Any) -> str:
    raw = json.dumps(parts, sort_keys=True, default=str).encode("utf-8")
    return f"vmh:{_CACHE_VERSION}:{hashlib.sha256(raw).hexdigest()}"


def _state_key(name: str) -> str:
    return f"vmh:{_CACHE_VERSION}:state:{name}"


def _cache_get(key: str, default: Any = None) -> Any:
    try:
        return _cache().get(key, default)
    except Exception as exc:  # a broken cache must never break VMH lookups
        logger.warning("VMH cache read failed (%s); continuing without the cache", exc)
        return default


def _cache_set(key: str, value: Any, timeout: int) -> None:
    try:
        _cache().set(key, value, timeout)
    except Exception as exc:
        logger.warning("VMH cache write failed (%s)", exc)


def _is_empty(value: Any) -> bool:
    return value is None or value is False or value == [] or value == {}


def _remember(key: str, value: Any) -> None:
    ttl = (_int_setting("VMH_NEGATIVE_CACHE_TTL", DEFAULT_NEGATIVE_CACHE_TTL)
           if _is_empty(value) else _int_setting("VMH_CACHE_TTL", DEFAULT_CACHE_TTL))
    if ttl > 0:
        _cache_set(key, {"value": value}, ttl)


def _memoize(kind: str, ident: Any, compute: Callable[[], Any]) -> Any:
    """Cache ``compute()``'s answer under (API base, kind, ident).  Exceptions
    propagate and are not cached."""
    key = _cache_key(vmh_api_base_url(), kind, ident)
    hit = _cache_get(key, _MISS)
    if isinstance(hit, dict) and "value" in hit:
        return hit["value"]
    value = compute()
    _remember(key, value)
    return value


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------

_thread_state = threading.local()


def _session() -> requests.Session:
    session = getattr(_thread_state, "session", None)
    if session is None:
        session = requests.Session()
        _thread_state.session = session
    return session


def _send(method: str, url: str, *, params=None, data=None, headers=None,
          timeout: float = DEFAULT_TIMEOUT) -> requests.Response:
    """The single network seam (tests patch this)."""
    return _session().request(method, url, params=params, data=data,
                              headers=headers, timeout=timeout)


def _clamp_seconds(seconds: float) -> int:
    return max(1, min(int(math.ceil(seconds)), 3600))


def _retry_after(resp: requests.Response) -> int:
    header = (resp.headers.get("Retry-After") or "").strip()
    if header.isdigit():
        return _clamp_seconds(int(header))
    try:
        reset = float(resp.headers.get("x-ratelimit-reset"))
    except (TypeError, ValueError):
        return 60
    if reset > 1e12:        # epoch milliseconds
        reset /= 1000.0
    if reset > 1e9:         # epoch seconds
        return _clamp_seconds(reset - time.time())
    return _clamp_seconds(reset)  # seconds from now


def _blocked_for() -> int:
    until = _cache_get(_state_key("blocked-until")) or 0
    try:
        remaining = float(until) - time.time()
    except (TypeError, ValueError):
        return 0
    return int(math.ceil(remaining)) if remaining > 0 else 0


def _block(seconds: int) -> None:
    _cache_set(_state_key("blocked-until"), time.time() + seconds, int(seconds) + 5)


def _note_rate_limit_headers(resp: requests.Response, sent_key: bool) -> None:
    limit = resp.headers.get("x-ratelimit-limit")
    if limit is None:
        return
    if sent_key and not _cache_get(_state_key("key-ignored-logged")):
        logger.warning(
            "VMH treated the configured VMH_API_KEY as anonymous (limit %s requests/hour "
            "applies). Check that the key is valid and of the 'public_unlimited' type.",
            limit)
        _cache_set(_state_key("key-ignored-logged"), True, 3600)
    remaining = str(resp.headers.get("x-ratelimit-remaining") or "")
    if remaining.isdigit() and int(remaining) <= 10:
        logger.warning(
            "VMH API: %s of %s anonymous requests left this hour. Set VMH_API_KEY to "
            "lift the limit.", remaining, limit)


def _json_body(resp: requests.Response) -> Any:
    try:
        return resp.json()
    except ValueError:
        return None


def _request(path: str, *, params: Optional[Dict[str, Any]] = None,
             data: Optional[Dict[str, Any]] = None, method: str = "GET",
             timeout: Optional[float] = None) -> Any:
    """Call the VMH API.

    Returns the parsed JSON for a 2xx answer and ``None`` for VMH's JSON
    "not found" (HTTP 404).  Raises VMHUnavailableError for everything else.
    """
    wait = _blocked_for()
    if wait:
        raise VMHUnavailableError("VMH's request limit has been reached",
                                  retry_after=wait, status=429)

    url = vmh_api_base_url() + path
    base_headers = {"Accept": "application/json"}
    key = _vmh_api_key()
    sent_key = bool(key) and not _cache_get(_state_key("key-rejected"))
    timeout = timeout or _int_setting("VMH_API_TIMEOUT", DEFAULT_TIMEOUT)
    try:
        headers = dict(base_headers, **({"x-api-key": key} if sent_key else {}))
        resp = _send(method, url, params=params, data=data, headers=headers, timeout=timeout)
        if sent_key and resp.status_code in (401, 403):
            logger.error(
                "VMH rejected VMH_API_KEY (HTTP %s); retrying without it. Replace the key "
                "in the server environment.", resp.status_code)
            _cache_set(_state_key("key-rejected"), True, 3600)
            sent_key = False
            resp = _send(method, url, params=params, data=data, headers=base_headers,
                         timeout=timeout)
    except requests.Timeout as exc:
        logger.warning("VMH request timed out after %ss: %s", timeout, url)
        raise VMHUnavailableError("VMH did not respond in time") from exc
    except requests.RequestException as exc:
        logger.warning("VMH request failed: %s (%s)", url, exc)
        raise VMHUnavailableError("could not connect to VMH") from exc

    _note_rate_limit_headers(resp, sent_key)
    status = resp.status_code
    if status == 429:
        retry = _retry_after(resp)
        _block(retry)
        logger.warning("VMH rate limit reached (HTTP 429); pausing VMH requests for %ss", retry)
        raise VMHUnavailableError("VMH's request limit has been reached",
                                  retry_after=retry, status=429)

    body = _json_body(resp)
    if status == 404:
        if isinstance(body, dict):
            return None  # VMH's JSON "not found" envelope
        logger.error("VMH API endpoint not found (non-JSON 404): %s - has the VMH API moved?",
                     url)
        raise VMHUnavailableError(
            "the VMH API endpoint was not found (the VMH API may have changed)", status=404)
    if status >= 500:
        logger.warning("VMH server error HTTP %s for %s", status, url)
        raise VMHUnavailableError(f"VMH returned a server error (HTTP {status})",
                                  status=status)
    if status >= 400 or status < 200 or status >= 300:
        logger.error("VMH rejected the request (HTTP %s) %s: %s", status, url,
                     body if body is not None else resp.text[:300])
        raise VMHUnavailableError(f"VMH rejected the request (HTTP {status})", status=status)
    if body is None:
        logger.error("VMH returned a non-JSON response for %s (content-type %s)", url,
                     resp.headers.get("content-type"))
        raise VMHUnavailableError("VMH returned an unexpected response", status=status)
    return body


def _search_all(path: str, list_key: str, params: Dict[str, Any],
                normalize: Callable[[Dict[str, Any]], Dict[str, Any]], *,
                page_length: int = 1000, max_pages: int = 10,
                timeout: Optional[float] = None) -> List[Dict[str, Any]]:
    """Collect every page of an advancedSearch listing."""
    rows: List[Dict[str, Any]] = []
    for page in range(1, max_pages + 1):
        data = _request(path, params=dict(params, page=page, pageLength=page_length),
                        timeout=timeout)
        if not isinstance(data, dict):
            break
        batch = [row for row in (data.get(list_key) or []) if isinstance(row, dict)]
        rows.extend(normalize(row) for row in batch)
        total = data.get("totalRecords")
        if not batch or not isinstance(total, int) or page * page_length >= total:
            break
    return rows


# ---------------------------------------------------------------------------
# SPARQL
# ---------------------------------------------------------------------------

def _sparql_literal(value: str) -> str:
    escaped = (str(value).replace("\\", "\\\\").replace('"', '\\"')
               .replace("\n", "\\n").replace("\r", "\\r").replace("\t", "\\t"))
    return f'"{escaped}"'


def _sparql_select(query: str, *, timeout: float = SPARQL_TIMEOUT) -> List[Dict[str, str]]:
    """Run a SPARQL SELECT (POSTed, so long SMILES literals are fine) and return
    its rows as plain ``{variable: value}`` dicts."""
    data = _request("/sparql", method="POST", data={"query": query, "format": "json"},
                    timeout=timeout)
    bindings = ((data or {}).get("results") or {}).get("bindings") if isinstance(data, dict) else None
    if not isinstance(bindings, list):
        raise VMHUnavailableError("VMH's SPARQL endpoint returned an unexpected response")
    return [
        {name: str(term.get("value", "")) for name, term in binding.items()
         if isinstance(term, dict)}
        for binding in bindings if isinstance(binding, dict)
    ]


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

def _normalize_metabolite(row: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "abbreviation": row.get("abbreviation") or "",
        "fullName": row.get("fullName") or row.get("name") or "",
        "smile": row.get("smile") or row.get("smiles") or "",
        "inchiString": row.get("inchiString") or "",
        "inchiKey": row.get("inchiKey") or "",
        "charge": row.get("charge"),
        "chargedFormula": row.get("chargedFormula") or "",
        "pubChemId": row.get("pubChemId") or "",
        "synonyms": row.get("synonyms") or "",
    }


def _normalize_reaction(row: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "abbreviation": row.get("abbreviation") or "",
        "description": row.get("description") or row.get("name") or "",
        "formula": row.get("formula") or "",
        "subsystem": row.get("subsystem") or "",
    }


def _normalize_gene_row(row: Dict[str, Any]) -> Dict[str, Any]:
    normalized = dict(row)
    entrez = row.get("entrez_id") or row.get("gene_number") or ""
    normalized["entrez_id"] = str(entrez) if entrez else ""
    # Historical name of the Entrez ID in the retired /_api/genes endpoint.
    normalized["gene_number"] = normalized["entrez_id"]
    return normalized


def _metabolite_rows_from_sparql(rows: List[Dict[str, str]], **extra: Any) -> List[Dict[str, Any]]:
    """One normalised record per abbreviation, sorted by abbreviation (the order
    the previous API returned matches in)."""
    by_abbr: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        abbr = (row.get("abbr") or "").strip()
        if abbr and abbr not in by_abbr:
            by_abbr[abbr] = _normalize_metabolite(
                dict(extra, abbreviation=abbr, fullName=row.get("name") or ""))
    return [by_abbr[abbr] for abbr in sorted(by_abbr)]


# ---------------------------------------------------------------------------
# Metabolites (strict: raise VMHUnavailableError when VMH cannot be asked)
# ---------------------------------------------------------------------------

def lookup_metabolite(abbr: str) -> Optional[Dict[str, Any]]:
    """Metabolite by abbreviation (``/getMetHandler``; case-insensitive).

    A *former* abbreviation resolves to the current record, flagged with
    ``is_former=True``; ``abbreviation`` is then the current one.
    """
    abbr = (abbr or "").strip()
    if not abbr:
        return None
    return _memoize("met-abbr", abbr, lambda: _fetch_metabolite(abbr))


def _fetch_metabolite(abbr: str) -> Optional[Dict[str, Any]]:
    data = _request("/getMetHandler", params={"id": abbr})
    if not isinstance(data, dict):
        return None
    rows = [row for row in (data.get("metData") or data.get("metaboliteData") or [])
            if isinstance(row, dict) and row.get("abbreviation")]
    if not rows:
        return None
    wanted = abbr.lower()
    row = (next((r for r in rows if r["abbreviation"] == abbr), None)
           or next((r for r in rows if r["abbreviation"].lower() == wanted), None))
    if row is not None:
        record = _normalize_metabolite(row)
        record["is_former"] = False
        return record
    former = {str(item).lower() for item in (data.get("formerAbbrs") or [])}
    if wanted in former:
        record = _normalize_metabolite(rows[0])
        record["is_former"] = True
        return record
    return None


def metabolite_abbreviation_exists(abbr: str) -> bool:
    """True if VMH knows ``abbr`` as a current *or former* metabolite abbreviation."""
    return bool(lookup_metabolite(abbr))


def search_metabolites_by_name(term: str) -> List[Dict[str, Any]]:
    """Metabolites whose abbreviation or full name contains ``term``
    (case-insensitive substring search, ``advancedSearch?searchTerm=``)."""
    term = (term or "").strip()
    if not term:
        return []
    return _memoize("met-search", term.lower(), lambda: _search_all(
        "/advancedSearch/metabolite", "metabolites", {"searchTerm": term},
        _normalize_metabolite, page_length=10000, max_pages=3))


def lookup_metabolite_by_name(name: str) -> Optional[Dict[str, Any]]:
    """Metabolite whose full name equals ``name`` (case-insensitive)."""
    name = (name or "").strip()
    if not name:
        return None
    wanted = name.lower()
    for row in search_metabolites_by_name(name):
        if (row.get("fullName") or "").strip().lower() == wanted:
            return row
    return None


def metabolite_name_exists(name: str) -> bool:
    return bool(lookup_metabolite_by_name(name))


def _advanced_metabolite_regex_search(field: str, value: str, *,
                                      case_sensitive: bool) -> List[Dict[str, Any]]:
    """Slow fallback for structure lookups when SPARQL is unavailable."""
    param = {"inchiKey": "inchiKeyRegex", "smile": "smileRegex"}[field]
    data = _request("/advancedSearch/metabolite",
                    params={param: f"^{_js_regex_escape(value)}$", "page": 1, "pageLength": 100},
                    timeout=SLOW_SEARCH_TIMEOUT)
    rows = (data.get("metabolites") or []) if isinstance(data, dict) else []

    def same(candidate: str) -> bool:
        candidate = (candidate or "").strip()
        return candidate == value if case_sensitive else candidate.upper() == value.upper()

    matches = [_normalize_metabolite(row) for row in rows
               if isinstance(row, dict) and row.get("abbreviation") and same(row.get(field))]
    return sorted(matches, key=lambda row: row["abbreviation"])


def _metabolites_with(prop: str, value: str, *, field: str, case_sensitive: bool,
                      extra: Dict[str, Any]) -> List[Dict[str, Any]]:
    query = ("SELECT DISTINCT ?abbr ?name WHERE { ?m a vmh:Metabolite ; "
             f"vmh:{prop} {_sparql_literal(value)} ; vmh:abbreviation ?abbr . "
             "OPTIONAL { ?m rdfs:label ?name } }")
    try:
        return _metabolite_rows_from_sparql(_sparql_select(query), **extra)
    except VMHUnavailableError as exc:
        if exc.rate_limited:
            raise
        logger.warning("VMH SPARQL lookup by %s failed (%s); falling back to the slow "
                       "advancedSearch filter", prop, exc.reason)
        return _advanced_metabolite_regex_search(field, value, case_sensitive=case_sensitive)


def lookup_metabolite_by_inchikey(inchi_key: str) -> Optional[Dict[str, Any]]:
    """Metabolite with this standard InChIKey: VMH's own InChIKeys first, then
    keys derived from the SMILES of metabolites VMH stores without one."""
    key = (inchi_key or "").strip().upper()
    if not _INCHIKEY_RE.match(key):
        return None
    cache_key = _cache_key(vmh_api_base_url(), "met-inchikey", key)
    hit = _cache_get(cache_key, _MISS)
    if isinstance(hit, dict) and "value" in hit:
        return hit["value"]

    rows = _metabolites_with("inchiKey", key, field="inchiKey", case_sensitive=False,
                             extra={"inchiKey": key})
    if rows:
        _remember(cache_key, rows[0])
        return rows[0]

    index = _derived_inchikey_index()
    if index is None:
        return None  # incomplete answer: do not cache "not found"
    for abbr in index.get(key) or []:
        record = lookup_metabolite(abbr)
        if record:
            record = dict(record, inchiKey=record.get("inchiKey") or key)
            _remember(cache_key, record)
            return record
    _remember(cache_key, None)
    return None


def lookup_metabolite_by_smiles(smiles: str) -> Optional[Dict[str, Any]]:
    """Metabolite whose stored SMILES is exactly ``smiles``."""
    smiles = (smiles or "").strip()
    if not smiles:
        return None

    def compute():
        rows = _metabolites_with("smiles", smiles, field="smile", case_sensitive=True,
                                 extra={"smile": smiles})
        return rows[0] if rows else None

    return _memoize("met-smiles", smiles, compute)


def match_metabolite_structure(*, inchi_key: str = "", inchi_string: str = "",
                               smiles: str = "") -> Optional[Dict[str, Any]]:
    """Find a VMH metabolite by structure: InChIKey (computed from the InChI
    when no key is given), then exact SMILES."""
    key = (inchi_key or "").strip().upper()
    if not key and inchi_string:
        key = _inchikey_from_inchi(inchi_string)
    if key:
        row = lookup_metabolite_by_inchikey(key)
        if row:
            return row
    if smiles:
        return lookup_metabolite_by_smiles(smiles)
    return None


# ---------------------------------------------------------------------------
# Derived InChIKey index (metabolites VMH stores with a SMILES but no InChIKey)
# ---------------------------------------------------------------------------

# The index lives in its own JSON file (not the Django cache, whose culling
# could evict it) next to a lock file, in VMH_CACHE_DIR.
_index_memory: Dict[str, Any] = {"by_key": None, "built_at": 0.0, "base": None}
_index_refreshing = threading.Event()


def _rdkit():
    from rdkit import Chem, rdBase  # local import: RDKit is heavy and only needed here
    return Chem, rdBase


def _inchikey_from_inchi(inchi_string: str) -> str:
    Chem, rdBase = _rdkit()
    try:
        with rdBase.BlockLogs():
            return Chem.InchiToInchiKey(inchi_string) or ""
    except Exception:
        return ""


def _inchikey_from_smiles(smiles: str) -> str:
    Chem, _ = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return ""
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            return ""
    try:
        return Chem.MolToInchiKey(mol) or ""
    except Exception:
        return ""


def _index_dir() -> str:
    directory = str(_setting("VMH_CACHE_DIR", "")
                    or os.path.join(tempfile.gettempdir(), "constructor-cache", "vmh"))
    os.makedirs(directory, exist_ok=True)
    return directory


def _index_path() -> str:
    return os.path.join(_index_dir(), f"derived-inchikey-index-v{_CACHE_VERSION}.json")


def _lock_path() -> str:
    return os.path.join(_index_dir(), "derived-inchikey-index.lock")


def _read_index_file():
    """(by_key, built_at) from disk, or (None, 0.0) if absent, unreadable,
    expired or built for another VMH API base URL."""
    try:
        with open(_index_path(), encoding="utf-8") as handle:
            data = json.load(handle)
        by_key, built_at = data["by_key"], float(data["built_at"])
        if (isinstance(by_key, dict) and data.get("base") == vmh_api_base_url()
                and time.time() - built_at < INDEX_CACHE_TTL):
            return by_key, built_at
    except (OSError, ValueError, KeyError, TypeError):
        pass
    return None, 0.0


def _write_index_file(by_key: Dict[str, List[str]], built_at: float) -> None:
    path = _index_path()
    temporary = f"{path}.{os.getpid()}.{threading.get_ident()}.tmp"
    try:
        with open(temporary, "w", encoding="utf-8") as handle:
            json.dump({"base": vmh_api_base_url(), "built_at": built_at, "by_key": by_key},
                      handle)
        os.replace(temporary, path)  # atomic: readers never see a partial file
    except OSError as exc:
        logger.warning("Could not save the derived VMH InChIKey index (%s)", exc)
        try:
            os.remove(temporary)
        except OSError:
            pass


def _remember_index(by_key: Dict[str, List[str]], built_at: float) -> Dict[str, List[str]]:
    _index_memory.update(by_key=by_key, built_at=built_at, base=vmh_api_base_url())
    return by_key


@contextmanager
def _index_build_lock(blocking: bool, timeout: float = 150.0) -> Iterator[bool]:
    """Cross-process lock so only one gunicorn worker builds the index."""
    if fcntl is None:
        yield True
        return
    fd = os.open(_lock_path(), os.O_RDWR | os.O_CREAT, 0o644)
    acquired = False
    try:
        deadline = time.time() + timeout
        while True:
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                acquired = True
                break
            except OSError:
                if not blocking or time.time() >= deadline:
                    break
                time.sleep(0.5)
        yield acquired
    finally:
        if acquired:
            fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)


_KEYLESS_METABOLITES_QUERY = (
    "SELECT DISTINCT ?abbr ?smiles WHERE { ?m a vmh:Metabolite ; vmh:abbreviation ?abbr ; "
    "vmh:smiles ?smiles . OPTIONAL { ?m vmh:inchiKey ?key } FILTER(!BOUND(?key)) }")


def _download_derived_index() -> Dict[str, List[str]]:
    # JSON results: VMH's CSV/TSV writers misalign columns when a variable is unbound.
    rows = _sparql_select(_KEYLESS_METABOLITES_QUERY, timeout=BULK_TIMEOUT)
    smiles_by_abbr: Dict[str, set] = {}
    for row in rows:
        abbr = (row.get("abbr") or "").strip()
        smiles = (row.get("smiles") or "").strip()
        if abbr and smiles:
            smiles_by_abbr.setdefault(abbr, set()).add(smiles)

    _, rdBase = _rdkit()
    by_key: Dict[str, List[str]] = {}
    with rdBase.BlockLogs():  # ~10% of VMH SMILES are R-group structures RDKit rejects
        for abbr in sorted(smiles_by_abbr):
            for smiles in sorted(smiles_by_abbr[abbr]):
                key = _inchikey_from_smiles(smiles)
                if key and abbr not in by_key.setdefault(key, []):
                    by_key[key].append(abbr)
    if not by_key:
        logger.warning("The derived VMH InChIKey index is empty (%d keyless metabolites "
                       "exported)", len(smiles_by_abbr))
    return by_key


def _build_derived_index(blocking: bool) -> Optional[Dict[str, List[str]]]:
    if _cache_get(_state_key("index-build-failed")):
        return None
    with _index_build_lock(blocking) as acquired:
        if not acquired:
            return None
        # Another worker may have finished a build while we waited for the lock.
        by_key, built_at = _read_index_file()
        if by_key is not None and time.time() - built_at < INDEX_MAX_AGE:
            return _remember_index(by_key, built_at)
        started = time.time()
        try:
            by_key = _download_derived_index()
        except VMHUnavailableError as exc:
            logger.warning("Could not build the derived VMH InChIKey index: %s", exc.reason)
            _cache_set(_state_key("index-build-failed"), True, INDEX_RETRY_AFTER_FAILURE)
            return None
        built_at = time.time()
        _write_index_file(by_key, built_at)
        logger.info("Built the derived VMH InChIKey index (%d keys) in %.1fs",
                    len(by_key), built_at - started)
        return _remember_index(by_key, built_at)


def _refresh_derived_index_in_background() -> None:
    if _index_refreshing.is_set():
        return
    _index_refreshing.set()

    def run():
        try:
            _build_derived_index(blocking=False)
        except Exception:  # pragma: no cover - defensive: never kill the worker
            logger.exception("Background refresh of the derived VMH InChIKey index failed")
        finally:
            _index_refreshing.clear()

    threading.Thread(target=run, name="vmh-derived-index", daemon=True).start()


def _derived_inchikey_index() -> Optional[Dict[str, List[str]]]:
    """InChIKey -> [abbreviation] for VMH metabolites without a stored InChIKey,
    or ``None`` when the index is unavailable (so "not found" is not final).

    Needs a VMH API key: anonymous SPARQL answers are silently truncated
    (observed: a few dozen to ~900 of ~11,000 rows), which would yield a
    partial index that looks complete."""
    if not _bool_setting("VMH_DERIVED_INDEX", True) or not vmh_api_key_configured():
        return {}
    now = time.time()
    if (_index_memory["by_key"] is not None and _index_memory["base"] == vmh_api_base_url()
            and now - _index_memory["built_at"] < INDEX_MAX_AGE):
        return _index_memory["by_key"]
    by_key, built_at = _read_index_file()
    if by_key is not None:
        _remember_index(by_key, built_at)
        if now - built_at >= INDEX_MAX_AGE:
            _refresh_derived_index_in_background()  # serve the stale index meanwhile
        return by_key
    return _build_derived_index(blocking=True)


def warm_derived_inchikey_index() -> int:
    """Build (or load) the derived index now; returns its number of keys."""
    index = _derived_inchikey_index()
    if index is None:
        raise VMHUnavailableError("the derived InChIKey index could not be built")
    return len(index)


# ---------------------------------------------------------------------------
# Reactions (strict)
# ---------------------------------------------------------------------------

def lookup_reaction(abbr: str) -> Optional[Dict[str, Any]]:
    """Reaction by abbreviation via ``/getRxnHandler``; former abbreviations
    resolve and are flagged ``is_former=True``.

    VMH matches abbreviations case-sensitively (the retired API did not), so the
    upper- and lower-case spellings are tried too ('hex1' -> 'HEX1')."""
    abbr = (abbr or "").strip()
    if not abbr:
        return None

    def compute():
        for variant in dict.fromkeys((abbr, abbr.upper(), abbr.lower())):
            record = _fetch_reaction(variant)
            if record:
                return record
        return None

    return _memoize("rxn-abbr", abbr, compute)


def _fetch_reaction(abbr: str) -> Optional[Dict[str, Any]]:
    wanted = abbr.lower()
    data = _request("/getRxnHandler", params={"id": abbr})
    if not isinstance(data, dict):
        return None
    rows = [row for row in (data.get("rxnData") or data.get("reactionData") or [])
            if isinstance(row, dict) and row.get("abbreviation")]
    if not rows:
        return None
    record = _normalize_reaction(rows[0])
    if not record["formula"] and isinstance(data.get("formula"), str):
        record["formula"] = data["formula"]
    if record["abbreviation"].lower() == wanted:
        record["is_former"] = False
        return record
    former = {str(item).lower() for item in (data.get("formerRxnAbbrs") or [])}
    if wanted in former or data.get("redirectAbbr"):
        record["is_former"] = True
        return record
    return None


def reaction_abbreviation_exists(abbr: str) -> bool:
    """True if VMH knows ``abbr`` as a reaction abbreviation: as given, in upper
    or lower case, or as a former abbreviation."""
    return bool(lookup_reaction(abbr))


def search_reactions_by_text(term: str) -> List[Dict[str, Any]]:
    """Reactions whose abbreviation or description contains ``term``."""
    term = (term or "").strip()
    if not term:
        return []
    return _memoize("rxn-search", term.lower(), lambda: _search_all(
        "/advancedSearch/reaction", "reactions", {"searchTerm": term}, _normalize_reaction,
        page_length=5000, max_pages=3))


def reaction_description_exists(description: str) -> bool:
    """True if a VMH reaction has exactly this description (case-insensitive)."""
    wanted = (description or "").strip().lower()
    if not wanted:
        return False
    return any((row.get("description") or "").strip().lower() == wanted
               for row in search_reactions_by_text(description))


def _js_regex_escape(text: str) -> str:
    return _JS_REGEX_SPECIAL.sub(r"\\\1", str(text))


def build_formula_regex(substrates: Sequence[str], products: Sequence[str]) -> str:
    """JavaScript regex for VMH formulas such as
    ``atp[c] + glc_D[c] -> adp[c] + g6p[c] + h[c]`` that have every substrate
    before the arrow, every product after it, and exactly as many terms on each
    side.  Stoichiometry and compartments are left to the caller to compare."""
    def side_shape(count: int, *, left: bool) -> str:
        if count == 0:
            return r"\s*" if left else r"\s*$"
        shape = r"[^+]*" + (r"(?:\+[^+]*){%d}" % (count - 1) if count > 1 else "")
        return shape if left else shape + "$"

    parts = ["^(?=" + side_shape(len(substrates), left=True) + _ARROW
             + side_shape(len(products), left=False) + ")"]
    for abbr in dict.fromkeys(substrates):
        parts.append(r"(?=(?:(?!%s).)*(?:^|\s)%s\[)" % (_ARROW, _js_regex_escape(abbr)))
    for abbr in dict.fromkeys(products):
        parts.append(r"(?=.*%s.*\s%s\[)" % (_ARROW, _js_regex_escape(abbr)))
    return "".join(parts)


def search_reactions_by_participants(substrates: Sequence[str],
                                     products: Sequence[str]) -> List[Dict[str, Any]]:
    """VMH reactions with exactly these metabolites (abbreviations) on each side."""
    pattern = build_formula_regex(list(substrates), list(products))
    return _memoize("rxn-participants", pattern, lambda: _search_all(
        "/advancedSearch/reaction", "reactions", {"formulaRegex": pattern},
        _normalize_reaction, page_length=500, max_pages=10))


# ---------------------------------------------------------------------------
# Genes (strict)
# ---------------------------------------------------------------------------

def lookup_gene(identifier: str) -> Optional[Dict[str, Any]]:
    """Gene by HGNC symbol or Entrez ID (VMH's transcript suffix such as
    ``2645.1`` is accepted).  Returns ``{'rows': [...], 'reactions': [...]}``."""
    identifier = (identifier or "").strip()
    if not identifier:
        return None
    return _memoize("gene", identifier, lambda: _fetch_gene(identifier))


def _gene_symbols_for_entrez(entrez_id: str) -> List[str]:
    query = ("SELECT DISTINCT ?symbol WHERE { ?g a vmh:Gene ; "
             f"vmh:entrezId {_sparql_literal(entrez_id)} ; vmh:symbol ?symbol }}")
    return sorted({row["symbol"] for row in _sparql_select(query) if row.get("symbol")})


def _fetch_gene(identifier: str) -> Optional[Dict[str, Any]]:
    # /getGeneHandler takes HGNC symbols only; map Entrez IDs through SPARQL.
    symbols = ([identifier] if not _ENTREZ_RE.match(identifier)
               else _gene_symbols_for_entrez(identifier.split(".", 1)[0]))
    for symbol in symbols:
        data = _request("/getGeneHandler", params={"id": symbol})
        if not isinstance(data, dict):
            continue
        rows = [_normalize_gene_row(row) for row in (data.get("geneData") or [])
                if isinstance(row, dict)]
        if rows:
            reactions = [_normalize_reaction(row) for row in (data.get("rxnData") or [])
                         if isinstance(row, dict)]
            return {"rows": rows, "reactions": reactions}
    return None


# ---------------------------------------------------------------------------
# Never-raise helpers (historical contract of this module)
# ---------------------------------------------------------------------------

def _quietly(compute: Callable[[], Any], default: Any, what: str) -> Any:
    try:
        return compute()
    except VMHUnavailableError as exc:
        logger.warning("VMH unavailable during %s: %s", what, exc.reason)
        return default


def find_metabolite_by_abbreviation(abbr: str) -> Optional[Dict[str, Any]]:
    """Current metabolite abbreviation (case-insensitive), else ``None``."""
    row = _quietly(lambda: lookup_metabolite(abbr), None, f"metabolite lookup {abbr!r}")
    return row if row and not row.get("is_former") else None


def find_metabolite_by_full_name(name: str) -> Optional[Dict[str, Any]]:
    return _quietly(lambda: lookup_metabolite_by_name(name), None,
                    f"metabolite name lookup {name!r}")


def find_metabolite_by_smiles(smiles: str) -> Optional[Dict[str, Any]]:
    return _quietly(lambda: lookup_metabolite_by_smiles(smiles), None, "SMILES lookup")


def find_metabolite_by_inchikey(inchi_key: str, *, inchi_string: str = "",
                                smiles: str = "") -> Optional[Dict[str, Any]]:
    return _quietly(lambda: match_metabolite_structure(
        inchi_key=inchi_key, inchi_string=inchi_string, smiles=smiles), None,
        "structure lookup")


def find_reaction_by_abbreviation(abbr: str, *, allow_former: bool = False
                                  ) -> Optional[Dict[str, Any]]:
    row = _quietly(lambda: lookup_reaction(abbr), None, f"reaction lookup {abbr!r}")
    if row and row.get("is_former") and not allow_former:
        return None
    return row


def reaction_name_exists(name: str) -> bool:
    return bool(_quietly(lambda: reaction_description_exists(name), False,
                         f"reaction name lookup {name!r}"))


def gene_data_new(gene_id: str) -> List[Dict[str, Any]]:
    """VMH gene rows for an HGNC symbol or Entrez ID (``[]`` when unknown)."""
    gene = _quietly(lambda: lookup_gene(str(gene_id)), None, f"gene lookup {gene_id!r}")
    return list(gene["rows"]) if gene else []


def gene_exists_new(gene_id: str) -> bool:
    return bool(gene_data_new(gene_id))


def gene_reaction_formulas(gene_id: str) -> List[str]:
    """Formulas of the VMH reactions associated with a gene."""
    gene = _quietly(lambda: lookup_gene(str(gene_id)), None, f"gene lookup {gene_id!r}")
    return [row["formula"] for row in (gene or {}).get("reactions", []) if row.get("formula")]


def gene_symbol_from_row(row: Dict[str, Any]) -> str:
    for key in ("symbol", "gene_symbol", "geneSymbol", "name"):
        value = row.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ""
