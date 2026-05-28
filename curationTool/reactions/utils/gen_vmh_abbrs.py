import random
import json
import os
import sys
import traceback
from pathlib import Path
from django.conf import settings
from reactions.utils.vmh_api import find_metabolite_by_abbreviation, find_reaction_by_abbreviation

# Use HTTP client to communicate with MATLAB container
from reactions.utils.MatlabHTTPClient import MatlabSessionManager

# Cache file path - stored in media folder which is mounted in docker
ABBR_CACHE_FILE = Path(settings.MEDIA_ROOT) / 'metabolite_abbr_cache.json'


def _load_abbr_cache():
    """Load the abbreviation cache from disk."""
    if ABBR_CACHE_FILE.exists():
        try:
            with open(ABBR_CACHE_FILE, 'r', encoding='utf-8') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"[WARN gen_vmh_abbrs] Failed to load cache: {e}", flush=True)
            return {}
    return {}


def _save_abbr_cache(cache):
    """Save the abbreviation cache to disk."""
    try:
        # Ensure directory exists
        ABBR_CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
        with open(ABBR_CACHE_FILE, 'w', encoding='utf-8') as f:
            json.dump(cache, f, indent=2, ensure_ascii=False)
    except IOError as e:
        print(f"[WARN gen_vmh_abbrs] Failed to save cache: {e}", flush=True)

def check_reaction_abbr_exists(abbr):
    row = find_reaction_by_abbreviation(abbr)
    return bool(row and (row.get("abbreviation") or "").lower() == abbr.lower())


def check_met_abbr_exists(abbr):
    row = find_metabolite_by_abbreviation(abbr)
    return bool(row and (row.get("abbreviation") or "").lower() == abbr.lower())


def gen_reaction_abbr(sub_abbr, prod_abbr, reaction):
    subs_comps = json.loads(reaction.subs_comps or '[]')
    prod_comps = json.loads(reaction.prods_comps or '[]')
    comp = (subs_comps or prod_comps or ['c'])[0]
    candidates = sub_abbr + prod_abbr
    if not candidates:
        raise ValueError('Cannot generate a reaction abbreviation without metabolites.')
    if len(candidates) == 1:
        abbr = candidates[0].upper() + comp
    else:
        abbr = random.choice(candidates).upper() + comp
    exists = check_reaction_abbr_exists(abbr)
    while exists:
        abbr = abbr + '_'
        exists = check_reaction_abbr_exists(abbr)
    return abbr


def gen_metabolite_abbr(
        metabolite,
        mtype,
        metabolite_name,
        search_func):
    if mtype == 'VMH':
        return metabolite

    # Check cache first (key by metabolite_name since that's what MATLAB uses)
    cache = _load_abbr_cache()
    cache_key = metabolite_name.strip().lower()
    
    if cache_key in cache:
        cached_abbr = cache[cache_key]
        print(f"[DEBUG gen_vmh_abbrs] Cache hit for '{metabolite_name}': {cached_abbr}", flush=True)
        return cached_abbr

    found, abbr = search_func(
        [metabolite], [mtype], None, side='substrates', nofile=True, return_abbr=True)
    found, abbr = found[0], abbr[0]
    if found:
        print(f"[DEBUG gen_vmh_abbrs] Metabolite '{metabolite_name}' found in VMH with abbr: {abbr}", flush=True)
        # Cache the VMH abbreviation too
        cache[cache_key] = abbr
        _save_abbr_cache(cache)
        return abbr
    else:
        print(f"[DEBUG gen_vmh_abbrs] Metabolite '{metabolite_name}' NOT found in VMH, calling MATLAB...", flush=True)
        matlab_session = MatlabSessionManager()
        print(f"[DEBUG gen_vmh_abbrs] MatlabSessionManager created, calling generateVMHMetAbbr...", flush=True)
        result = matlab_session.execute('generateVMHMetAbbr', metabolite_name)
        print(f"[DEBUG gen_vmh_abbrs] MATLAB result: {result}", flush=True)
        abbr = result['result'] if result['status'] == 'success' else metabolite_name
        abbr = abbr[-1] if isinstance(abbr, list) else abbr
        exists = check_met_abbr_exists(abbr)
        while exists:
            abbr = abbr + '_'
            exists = check_met_abbr_exists(abbr)
        print(f"[DEBUG gen_vmh_abbrs] Final abbreviation: {abbr}", flush=True)
        
        # Save to cache
        cache[cache_key] = abbr
        _save_abbr_cache(cache)
        
        return abbr
