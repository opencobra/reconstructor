"""Client for the ChEBI 2.0 REST API.

EBI retired the ChEBI SOAP web service
(https://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl now answers with an
HTML error page); compounds are served by
https://www.ebi.ac.uk/chebi/backend/api/public (docs: /chebi/backend/api/docs/).

``get_compound`` / ``find_compound_by_name`` return ``None`` when ChEBI has no
such entry and raise a ``requests.RequestException`` when ChEBI cannot be asked.
Each ChEBI request takes ~15 s (sometimes much longer) however many compounds it
asks for, so answers are cached and callers resolving several compounds should
first call ``prefetch_compounds`` (one batched request).
"""

import html
import logging
import re
from typing import Any, Dict, Optional, Tuple

import requests
from django.conf import settings
from django.core.cache import cache

logger = logging.getLogger(__name__)

CHEBI_API_URL = "https://www.ebi.ac.uk/chebi/backend/api/public"
TIMEOUT = 60
BATCH_TIMEOUT = 90
BATCH_SIZE = 50
CACHE_TTL = 24 * 3600
NEGATIVE_CACHE_TTL = 3600

_ID_RE = re.compile(r"^(?:chebi:?\s*)?(\d+)$", re.IGNORECASE)
_TAG_RE = re.compile(r"<[^>]+>")
_MISS = object()


def _base_url() -> str:
    return (getattr(settings, "CHEBI_API_BASE_URL", "") or CHEBI_API_URL).rstrip("/")


def normalize_chebi_id(value: Any) -> str:
    """'CHEBI:15377', 'chebi:15377', 'CHEBI15377' and '15377' -> 'CHEBI:15377'
    (anything else is returned stripped)."""
    text = str(value or "").strip()
    match = _ID_RE.match(text)
    return f"CHEBI:{match.group(1)}" if match else text


def _plain(name: Any) -> str:
    """Lower-cased name without ChEBI's HTML markup ('<small>D</small>-glucose')."""
    return _TAG_RE.sub("", html.unescape(str(name or ""))).strip().lower()


def _get(path: str, params: Optional[Dict[str, Any]] = None,
         timeout: float = TIMEOUT) -> Optional[Any]:
    response = requests.get(_base_url() + path, params=params, timeout=timeout,
                            headers={"Accept": "application/json"})
    if response.status_code == 404:
        return None
    response.raise_for_status()
    return response.json()


def _cache_lookup(key: str):
    try:
        hit = cache.get(key, _MISS)
    except Exception:  # a broken cache must not break lookups
        return _MISS
    return hit["value"] if isinstance(hit, dict) and "value" in hit else _MISS


def _cache_store(key: str, value) -> None:
    try:
        cache.set(key, {"value": value}, CACHE_TTL if value is not None else NEGATIVE_CACHE_TTL)
    except Exception as exc:
        logger.warning("ChEBI cache write failed (%s)", exc)


def _cached(key: str, compute):
    hit = _cache_lookup(key)
    if hit is not _MISS:
        return hit
    value = compute()
    _cache_store(key, value)
    return value


def _compound_key(accession: str) -> str:
    return f"chebi:1:compound:{accession}"


def _reduce_compound(data: Any) -> Optional[Dict[str, Any]]:
    if not isinstance(data, dict) or not data.get("chebi_accession"):
        return None
    structure = data.get("default_structure") or {}
    return {
        "chebi_accession": data["chebi_accession"],
        "name": data.get("ascii_name") or _plain(data.get("name")) or data["chebi_accession"],
        "inchi": structure.get("standard_inchi") or "",
        "inchi_key": structure.get("standard_inchi_key") or "",
        "smiles": structure.get("smiles") or "",
    }


def prefetch_compounds(chebi_ids) -> None:
    """Fetch the uncached compounds among ``chebi_ids`` with batched requests so
    later ``get_compound`` calls are answered from the cache.  Best effort:
    on failure the per-compound lookups simply ask ChEBI themselves."""
    wanted = []
    for value in chebi_ids or []:
        accession = normalize_chebi_id(value)
        if (accession.startswith("CHEBI:") and accession not in wanted
                and _cache_lookup(_compound_key(accession)) is _MISS):
            wanted.append(accession)
    for start in range(0, len(wanted), BATCH_SIZE):
        chunk = wanted[start:start + BATCH_SIZE]
        try:
            data = _get("/compounds/", params={"chebi_ids": ",".join(chunk)},
                        timeout=BATCH_TIMEOUT)
        except (requests.RequestException, ValueError) as exc:
            logger.warning("ChEBI batch lookup of %d compounds failed (%s); "
                           "looking them up one by one", len(chunk), exc)
            return
        entries = {}
        for key, entry in (data or {}).items() if isinstance(data, dict) else []:
            if isinstance(entry, dict):
                entries[normalize_chebi_id(key)] = entry
                entries[normalize_chebi_id(entry.get("standardized_chebi_id"))] = entry
        for accession in chunk:
            entry = entries.get(accession)
            if entry is None:
                continue  # not answered: leave it to the single lookup
            if entry.get("exists") is False:
                _cache_store(_compound_key(accession), None)
                continue
            compound = _reduce_compound(entry.get("data"))
            if compound:
                _cache_store(_compound_key(accession), compound)


def get_compound(chebi_id: Any) -> Optional[Dict[str, Any]]:
    """Compound by ChEBI ID (primary or secondary), reduced to what Constructor
    uses: accession, name, InChI, InChIKey, SMILES."""
    accession = normalize_chebi_id(chebi_id)
    if not accession.startswith("CHEBI:"):
        return None

    return _cached(_compound_key(accession),
                   lambda: _reduce_compound(_get(f"/compound/{accession}/")))


def find_compound_by_name(name: str) -> Optional[Dict[str, Any]]:
    """Best ChEBI hit for a name: an exact (case-insensitive) name match among
    the search results, else the single result of an exact synonym match, else
    the top-ranked result (ChEBI's ranking does not always put the exact name
    first, e.g. 'dioxygen').  Returns ``{'chebi_accession', 'name'}``."""
    term = (name or "").strip()
    if not term:
        return None

    def fetch():
        data = _get("/es_search/", params={"term": term, "size": 50})
        results = [item.get("_source", item) for item in ((data or {}).get("results") or [])
                   if isinstance(item, dict)]
        results = [item for item in results if item.get("chebi_accession")]
        if not results:
            return None
        wanted = term.lower()
        exact = [item for item in results
                 if _plain(item.get("ascii_name")) == wanted or _plain(item.get("name")) == wanted]
        chosen = (sorted(exact, key=lambda item: -(item.get("stars") or 0))[0] if exact
                  else results[0])
        return {"chebi_accession": chosen["chebi_accession"],
                "name": chosen.get("ascii_name") or _plain(chosen.get("name"))}

    return _cached(f"chebi:1:name:{term.lower()}", fetch)


def resolve(identifier: str, *, by_name: bool) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """(compound, error message) for a ChEBI ID or name."""
    if by_name:
        hit = find_compound_by_name(identifier)
        if not hit:
            return None, f"Metabolite `{identifier}` does not exist in ChEBI"
        identifier = hit["chebi_accession"]
    compound = get_compound(identifier)
    if not compound:
        return None, f"Metabolite {identifier} does not exist in ChEBI"
    return compound, None
