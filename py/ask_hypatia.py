import numpy as np
import requests
from astropy.io import ascii
import random

key = '136611fd615f8a90aafb1da408f0e5b3'

def retrieve_star_names(exo_hosts=True, API_KEY=key, writeto='host_names.txt', exclude_blank=False):

    catalog = requests.get("https://hypatiacatalog.com/hypatia/api/v2/data", auth=(API_KEY, "api_token"))
    names_all = []
    print(catalog.json()['values'])
    for row in catalog.json()['values']:  # list of dicts
        names_all.append(row['name'])

    if not exo_hosts and not exclude_blank:
        names_final = names_all
    else:
        params = {}
        if exo_hosts:
            params['name'] = names_all
            address = "https://hypatiacatalog.com/hypatia/api/v2/star"
        # if (not not exclude_blank):
        #     address = "https://hypatiacatalog.com/hypatia/api/v2/composition"  # works for just
        #     params['element'] = []
        #     for el in exclude_blank:
        #         params['element'].append(el)

        entry = requests.get(address, auth=(API_KEY, "api_token"),
                             params=params)

        names_final = []
        for ii, r in enumerate(entry.json()):  # list of dicts
            if not exo_hosts or (exo_hosts and len(r['planets']) > 0):
                # only add if more than 0 planets
                names_final.append(names_all[ii])
            # if not exclude_blank or (or not not exclude_blank and 'mg' in exclude_blank):
            #     # only add if not missing this measurement - todo

    if writeto is not None:
        with open(writeto, 'w') as filehandle:
            filehandle.writelines("%s\n" % a for a in names_final)
    return names_final


def random_star(n=1, names_file='host_names.txt', **kwargs):
    # get a random star from database - N.B. there are 1331 entries for exo hosts
    if names_file is not None:
        with open(names_file, 'r') as filehandle:
            names = [a.rstrip() for a in filehandle.readlines()]
    else:
        names = retrieve_star_names(**kwargs)

    if n == -1:  # return all stars
        return names  # in same order as file

    idx = random.sample(range(0, len(names)), n)
    sample_names = []
    for i in idx:
        sample_names.append(names[i])
    return sample_names


def star_composition(oxide_list=None, star_name='sun', API_KEY=key, verbose=False,
                ca_sol=6.33 - 12, al_sol=6.47 - 12, fe_sol=7.45 - 12, si_sol=7.52 - 12, mg_sol=7.54 - 12,
                na_sol=6.30 - 12, **kwargs):
    """ star id is same as hypatia catalog with spaces e.g. "HIP 12345" """

    els = []
    for ox in oxide_list:
        els.append(ox[:2].lower())

    if star_name != 'sun':  # not needed for solar values
        params = {"name": [star_name] * len(els), "element": els, "solarnorm": ["lod09"] * len(els)}
        entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", auth=(API_KEY, "api_token"),
                             params=params)

        if np.size(entry.json()) == 0:
            raise Exception('No entry found in Hypatia Catalog:', star_name)
        # print('loaded json', entry.json())
    nH_star = []
    for ii, el in enumerate(els):
        # absolute not working for some reason so get difference from solar via lodders norm
        if star_name == 'sun':
            nH = 0
        else:
            try:
                nH = entry.json()[ii]['mean']
                if verbose:
                    print('retrieving from hypatia...', entry.json()[ii]['element'], 'mean =', entry.json()[ii]['mean'])
                sol_val = eval(el + '_sol')
                nH_star.append(nH + sol_val)
            except TypeError:
                print(star_name, 'does not have measured', el, 'and should be removed from stellar sample')
                return None

    return nH_star  # will always be in same order as oxides list