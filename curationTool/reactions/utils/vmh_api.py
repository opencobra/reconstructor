"""Helpers for VMH API access (new API first, old API fallback).

This module centralizes:
- Hardcoded VMH API endpoints (new + old)
- x-api-key handling for new API
- URL builders for VMH web links
- Common search/get wrappers used across the codebase
"""

from __future__ import annotations

import logging
import os
import re
from typing import Any, Dict, List, Optional
from urllib.parse import quote

import requests
from django.conf import settings

logger = logging.getLogger(__name__)

# Hardcoded in one place, as requested.
VMH_NEW_API_BASE = "https://vmh2.life"
VMH_OLD_API_BASE = "https://www.vmh.life"

_DEFAULT_TIMEOUT = 10


def _vmh_api_key() -> str:
    return os.getenv("VMH_API_KEY", "") or getattr(settings, "VMH_API_KEY", "") or ""


def vmh_public_base_url() -> str:
    """Prefix used to build clickable VMH links such as /metabolite/{abbr}."""
    base = os.getenv("VMH_BASE_URL", "") or getattr(settings, "VMH_BASE_URL", "")
    base = (base or VMH_OLD_API_BASE).strip()
    return base.rstrip("/")


def vmh_metabolite_url(abbr: str) -> str:
    return f"{vmh_public_base_url()}/metabolite/{quote(str(abbr), safe='')}"


def vmh_reaction_url(abbr: str) -> str:
    return f"{vmh_public_base_url()}/reaction/{quote(str(abbr), safe='')}"


def _safe_get(
    url: str,
    *,
    params: Optional[Dict[str, Any]] = None,
    headers: Optional[Dict[str, str]] = None,
    timeout: int = _DEFAULT_TIMEOUT,
    verify: bool = True,
) -> Optional[requests.Response]:
    try:
        return requests.get(url, params=params, headers=headers, timeout=timeout, verify=verify)
    except Exception as exc:
        logger.warning("VMH API request failed: %s — %s", url, exc)
        return None


def vmh_new_get(path: str, *, params: Optional[Dict[str, Any]] = None, timeout: int = _DEFAULT_TIMEOUT) -> Optional[requests.Response]:
    headers = {"Accept": "application/json"}
    key = _vmh_api_key()
    if key:
        headers["x-api-key"] = key
    return _safe_get(
        f"{VMH_NEW_API_BASE}{path}",
        params=params,
        headers=headers,
        timeout=timeout,
        verify=True,
    )


def vmh_old_get(path: str, *, params: Optional[Dict[str, Any]] = None, timeout: int = _DEFAULT_TIMEOUT) -> Optional[requests.Response]:
    # Preserve legacy behavior (verify=False) used throughout this codebase.
    return _safe_get(
        f"{VMH_OLD_API_BASE}{path}",
        params=params,
        timeout=timeout,
        verify=False,
    )


def _json_or_empty(resp: Optional[requests.Response]) -> Dict[str, Any]:
    if not resp:
        return {}
    if resp.status_code != 200:
        logger.warning("VMH API returned %s for %s", resp.status_code, resp.url)
        return {}
    try:
        return resp.json()
    except Exception as exc:
        logger.warning("VMH API response parse error for %s: %s", resp.url, exc)
        return {}


def _normalize_new_metabolite(row: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "abbreviation": row.get("abbreviation") or row.get("id") or "",
        "fullName": row.get("fullName") or row.get("name") or "",
        "smile": row.get("smile") or "",
        "inchiString": row.get("inchiString") or "",
        "inchiKey": row.get("inchiKey") or "",
        "charge": row.get("charge"),
        "pubChemId": row.get("pubChemId") or "",
        "synonyms": row.get("synonyms") or "",
    }


def _normalize_old_metabolite(row: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "abbreviation": row.get("abbreviation") or "",
        "fullName": row.get("fullName") or "",
        "smile": row.get("smile") or "",
        "inchiString": row.get("inchiString") or "",
        "inchiKey": row.get("inchiKey") or "",
        "charge": row.get("charge"),
        "pubChemId": row.get("pubChemId") or "",
        "synonyms": row.get("synonyms") or "",
    }


def _normalize_new_reaction(row: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "abbreviation": row.get("abbreviation") or row.get("id") or "",
        "description": row.get("description") or row.get("name") or "",
        "formula": row.get("formula") or "",
        "subsystem": row.get("subsystem") or "",
    }


def _normalize_old_reaction(row: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "abbreviation": row.get("abbreviation") or "",
        "description": row.get("description") or "",
        "formula": row.get("formula") or "",
        "subsystem": row.get("subsystem") or "",
    }


def search_metabolites_new(params: Dict[str, Any], *, timeout: int = _DEFAULT_TIMEOUT) -> List[Dict[str, Any]]:
    query = {"page": 1, "pageLength": 10000, "colSort": "abbreviation", "sortOrder": "ASC"}
    query.update(params)
    resp = vmh_new_get("/advancedSearch/metabolite", params=query, timeout=timeout)
    data = _json_or_empty(resp)
    rows = data.get("metabolites") or []
    return [_normalize_new_metabolite(row) for row in rows if isinstance(row, dict)]


def search_metabolites_old(params: Dict[str, Any], *, timeout: int = _DEFAULT_TIMEOUT) -> List[Dict[str, Any]]:
    resp = vmh_old_get("/_api/metabolites/", params=params, timeout=timeout)
    data = _json_or_empty(resp)
    rows = data.get("results") or []
    return [_normalize_old_metabolite(row) for row in rows if isinstance(row, dict)]


def search_reactions_new(params: Dict[str, Any], *, timeout: int = _DEFAULT_TIMEOUT) -> List[Dict[str, Any]]:
    query = {"page": 1, "pageLength": 30000, "colSort": "abbreviation", "sortOrder": "ASC"}
    query.update(params)
    resp = vmh_new_get("/advancedSearch/reaction", params=query, timeout=timeout)
    data = _json_or_empty(resp)
    rows = data.get("reactions") or []
    return [_normalize_new_reaction(row) for row in rows if isinstance(row, dict)]


def search_reactions_old(params: Dict[str, Any], *, timeout: int = _DEFAULT_TIMEOUT) -> List[Dict[str, Any]]:
    resp = vmh_old_get("/_api/reactions/", params=params, timeout=timeout)
    data = _json_or_empty(resp)
    rows = data.get("results") or []
    return [_normalize_old_reaction(row) for row in rows if isinstance(row, dict)]


def find_metabolite_by_abbreviation(abbr: str) -> Optional[Dict[str, Any]]:
    rows = search_metabolites_new({"abbreviation": abbr})
    exact = [r for r in rows if (r.get("abbreviation") or "").lower() == abbr.lower()]
    if exact:
        return exact[0]

    rows_old = search_metabolites_old({"abbreviation": abbr})
    exact_old = [r for r in rows_old if (r.get("abbreviation") or "").lower() == abbr.lower()]
    return exact_old[0] if exact_old else None


def find_metabolite_by_full_name(name: str) -> Optional[Dict[str, Any]]:
    rows = search_metabolites_new({"fullName": name})
    exact = [r for r in rows if (r.get("fullName") or "").lower() == name.lower()]
    if exact:
        return exact[0]

    rows_old = search_metabolites_old({"fullName": name})
    exact_old = [r for r in rows_old if (r.get("fullName") or "").lower() == name.lower()]
    return exact_old[0] if exact_old else None


def find_metabolite_by_inchi_string_old(inchi_string: str) -> Optional[Dict[str, Any]]:
    if not inchi_string:
        return None
    rows_old = search_metabolites_old({"inchiString": inchi_string})
    for row in rows_old:
        if (row.get("inchiString") or "").strip() == inchi_string.strip():
            return row
    return rows_old[0] if rows_old else None


def find_metabolite_by_smiles(smiles: str) -> Optional[Dict[str, Any]]:
    if not smiles:
        return None
    safe_smiles = re.escape(smiles)
    rows = search_metabolites_new({"smileRegex": f"^{safe_smiles}$"})
    exact = [r for r in rows if (r.get("smile") or "") == smiles]
    return exact[0] if exact else None


def find_metabolite_by_inchikey(
    inchi_key: str,
    *,
    inchi_string: str = "",
    smiles: str = "",
) -> Optional[Dict[str, Any]]:
    if not inchi_key:
        if inchi_string:
            old_row = find_metabolite_by_inchi_string_old(inchi_string)
            if old_row:
                return old_row
        if smiles:
            return find_metabolite_by_smiles(smiles)
        return None

    safe_key = re.escape(inchi_key)
    rows = search_metabolites_new({"inchiKeyRegex": f"^{safe_key}$"})
    exact = [r for r in rows if (r.get("inchiKey") or "").upper() == inchi_key.upper()]
    if exact:
        return exact[0]

    if inchi_string:
        old_row = find_metabolite_by_inchi_string_old(inchi_string)
        if old_row:
            return old_row
    if smiles:
        smiles_row = find_metabolite_by_smiles(smiles)
        if smiles_row:
            return smiles_row
    return None


def find_reaction_by_abbreviation(abbr: str) -> Optional[Dict[str, Any]]:
    rows = search_reactions_new({"abbreviation": abbr})
    exact = [r for r in rows if (r.get("abbreviation") or "").lower() == abbr.lower()]
    if exact:
        return exact[0]

    rows_old = search_reactions_old({"abbreviation": abbr})
    exact_old = [r for r in rows_old if (r.get("abbreviation") or "").lower() == abbr.lower()]
    return exact_old[0] if exact_old else None


def reaction_name_exists(name: str) -> bool:
    rows = search_reactions_new({"searchTerm": name})
    if any((r.get("description") or "").lower() == name.lower() for r in rows):
        return True

    rows_old = search_reactions_old({"description": name})
    return any((r.get("description") or "").lower() == name.lower() for r in rows_old)


def gene_data_new(gene_id: str, *, timeout: int = _DEFAULT_TIMEOUT) -> List[Dict[str, Any]]:
    resp = vmh_new_get("/getGeneHandler", params={"id": gene_id}, timeout=timeout)
    data = _json_or_empty(resp)
    rows = data.get("geneData") or []
    return [row for row in rows if isinstance(row, dict)]


def gene_exists_new(gene_id: str) -> bool:
    return bool(gene_data_new(gene_id))


def gene_rows_old(*, symbol: str = "", gene_number: str = "", timeout: int = _DEFAULT_TIMEOUT) -> List[Dict[str, Any]]:
    params: Dict[str, Any] = {}
    if symbol:
        params["symbol"] = symbol
    if gene_number:
        params["gene_number"] = gene_number
    if not params:
        return []

    resp = vmh_old_get("/_api/genes/", params=params, timeout=timeout)
    data = _json_or_empty(resp)
    rows = data.get("results") or []
    return [row for row in rows if isinstance(row, dict)]


def gene_symbol_from_row(row: Dict[str, Any]) -> str:
    for key in ("symbol", "gene_symbol", "geneSymbol", "name"):
        value = row.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ""
