import random
import requests
import json
import os
import sys
import traceback
from django.conf import settings

# AGGRESSIVE DEBUG LOGGING
print(f"=" * 80, flush=True)
print(f"[DEBUG gen_vmh_abbrs] MODULE LOADING STARTED", flush=True)
print(f"[DEBUG gen_vmh_abbrs] Python executable: {sys.executable}", flush=True)
print(f"[DEBUG gen_vmh_abbrs] Python version: {sys.version}", flush=True)
print(f"[DEBUG gen_vmh_abbrs] sys.path: {sys.path}", flush=True)

# Check environment variables
matlab_remote_enabled = os.getenv('MATLAB_REMOTE_ENABLED', 'NOT_SET')
print(f"[DEBUG gen_vmh_abbrs] MATLAB_REMOTE_ENABLED = '{matlab_remote_enabled}'", flush=True)

# Check if matlabengine package exists
try:
    import matlab
    print(f"[DEBUG gen_vmh_abbrs] ✓ 'matlab' package found at: {matlab.__file__}", flush=True)
    try:
        import matlab.engine
        print(f"[DEBUG gen_vmh_abbrs] ✓ 'matlab.engine' module found at: {matlab.engine.__file__}", flush=True)
    except ImportError as e:
        print(f"[DEBUG gen_vmh_abbrs] ✗ 'matlab.engine' NOT found: {e}", flush=True)
except ImportError as e:
    print(f"[DEBUG gen_vmh_abbrs] ✗ 'matlab' package NOT found: {e}", flush=True)

# Check if MatlabSessionManagerRemote exists
import os.path
remote_manager_path = os.path.join(os.path.dirname(__file__), 'MatlabSessionManagerRemote.py')
print(f"[DEBUG gen_vmh_abbrs] Looking for MatlabSessionManagerRemote at: {remote_manager_path}", flush=True)
print(f"[DEBUG gen_vmh_abbrs] File exists? {os.path.exists(remote_manager_path)}", flush=True)

# Use remote MATLAB session manager if enabled, otherwise use local
skip = False
MatlabSessionManager = None

try:
    if matlab_remote_enabled.lower() == 'true':
        print(f"[DEBUG gen_vmh_abbrs] >>> Importing MatlabSessionManagerRemote...", flush=True)
        from reactions.utils.MatlabSessionManagerRemote import MatlabSessionManager
        print(f"[DEBUG gen_vmh_abbrs] >>> SUCCESS! MatlabSessionManager = {MatlabSessionManager}", flush=True)
    else:
        print(f"[DEBUG gen_vmh_abbrs] >>> Importing MatlabSessionManager (local)...", flush=True)
        from reactions.utils.MatlabSessionManager import MatlabSessionManager
        print(f"[DEBUG gen_vmh_abbrs] >>> SUCCESS! MatlabSessionManager = {MatlabSessionManager}", flush=True)
except Exception as e:
    print(f"[ERROR gen_vmh_abbrs] ✗✗✗ IMPORT FAILED ✗✗✗", flush=True)
    print(f"[ERROR gen_vmh_abbrs] Exception: {e}", flush=True)
    print(f"[ERROR gen_vmh_abbrs] Exception type: {type(e).__name__}", flush=True)
    print(f"[ERROR gen_vmh_abbrs] Full traceback:", flush=True)
    traceback.print_exc(file=sys.stdout)
    sys.stdout.flush()
    skip = True

print(f"[DEBUG gen_vmh_abbrs] Final state: MatlabSessionManager = {MatlabSessionManager}, skip = {skip}", flush=True)
print(f"[DEBUG gen_vmh_abbrs] MODULE LOADING COMPLETED", flush=True)
print(f"=" * 80, flush=True)
sys.stdout.flush()

def check_reaction_abbr_exists(abbr):
    BASE_URL = settings.OLD_VMH_BASE_URL
    endpoint = f"{BASE_URL}_api/reactions/?abbreviation={abbr}"
    # Make the GET request
    response = requests.get(endpoint, verify=False)
    if response.json().get('count', 0) == 0:
        return False
    else:
        return True


def check_met_abbr_exists(abbr):
    BASE_URL = settings.OLD_VMH_BASE_URL
    endpoint = f"{BASE_URL}_api/metabolites/?abbreviation={abbr}"
    # Make the GET request
    response = requests.get(endpoint, verify=False)
    if response.json().get('count', 0) == 0:
        return False
    else:
        return True


def gen_reaction_abbr(sub_abbr, prod_abbr, reaction):
    comp = json.loads(reaction.subs_comps)[0]
    if len(sub_abbr) == 1 and len(prod_abbr) == 1:
        abbr = (sub_abbr[0]).upper() + comp
    else:
        if random.random() > 0.5:
            abbr = sub_abbr[random.choice(range(len(sub_abbr)))].upper() + comp
        else:
            abbr = prod_abbr[random.choice(
                range(len(prod_abbr)))].upper() + comp
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

    found, abbr = search_func(
        [metabolite], [mtype], None, side='substrates', nofile=True, return_abbr=True)
    found, abbr = found[0], abbr[0]
    if found:
        return abbr
    else:
        if MatlabSessionManager is None:
            raise RuntimeError("MatlabSessionManager is not available. MATLAB integration is not configured.")
        
        matlab_session = MatlabSessionManager()
        result = matlab_session.execute('generateVMHMetAbbr', metabolite_name)
        abbr = result['result'] if result['status'] == 'success' else metabolite_name
        # abbr = abbr[-1] if isinstance(abbr, list) else abbr
        exists = check_met_abbr_exists(abbr)
        while exists:
            abbr = abbr + '_'
            exists = check_met_abbr_exists(abbr)
        return abbr
