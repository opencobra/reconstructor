import random
import requests
import json
from reactions.utils.MatlabSessionManager import MatlabSessionManager


def check_reaction_abbr_exists(abbr):
    BASE_URL = 'https://www.vmh.life/'
    endpoint = f"{BASE_URL}_api/reactions/?abbreviation={abbr}"
    # Make the GET request
    response = requests.get(endpoint, verify=False)
    if response.json().get('count', 0) == 0:
        return False
    else:
        return True


def check_met_abbr_exists(abbr):
    BASE_URL = 'https://www.vmh.life/'
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
        matlab_session = MatlabSessionManager()
        result = matlab_session.execute('generateVMHMetAbbr', metabolite_name)
        abbr = result['result'] if result['status'] == 'success' else metabolite_name
        # abbr = abbr[-1] if isinstance(abbr, list) else abbr
        exists = check_met_abbr_exists(abbr)
        while exists:
            abbr = abbr + '_'
            exists = check_met_abbr_exists(abbr)
        return abbr
