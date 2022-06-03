import numpy as np
import requests
# from astropy.io import ascii
import random
import pickle as pkl
import os

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
    print(os.getcwd())
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


def star_composition(oxide_list=None, star='sun', API_KEY=key, verbose=False, use_local_composition=False, **kwargs):
    """ star id is same as hypatia catalog with spaces e.g. "HIP 12345" """
    import parameters as p

    els = []
    for ox in oxide_list:
        els.append(ox[:2].lower())

    def do_local():
        try:
            path = get_directory(star)
            with open(path + '/dat.pkl', "rb") as pfile:
                dat = pkl.load(pfile)
            nH_star = dat.nH_star
        except FileNotFoundError:
            return None
        return nH_star

    params = {"name": [star] * len(els), "element": els, "solarnorm": ["lod09"] * len(els)}
    if not use_local_composition:
        try:
            entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", auth=(API_KEY, "api_token"),
                             params=params)

            if np.size(entry.json()) == 0:
                raise Exception('No entry found in Hypatia Catalog:', star)
            # print('loaded json', entry.json())
            nH_star = []
            for ii, el in enumerate(els):
                # absolute not working for some reason so get difference from solar via lodders norm
                try:
                    nH = entry.json()[ii]['mean']
                    if verbose:
                        print('retrieving from hypatia...', entry.json()[ii]['element'], 'mean =', entry.json()[ii]['mean'])
                    sol_val = eval('p.' + el + '_sol')
                    nH_star.append(nH + sol_val)
                except TypeError:
                    print(star, 'does not have measured', el, 'and should be removed from stellar sample')
                    return None
        except (ConnectionError, requests.exceptions.ConnectionError) as e:
            # try loading from file
            nH_star = do_local()
    else:
        nH_star = do_local()

    return nH_star  # will always be in same order as oxides list


def get_directory(star_name, existing_dir='hypatia1M_1600K_80Fe/'):
    from perplexdata import output_parent_default
    """retrieve directory of existing run for a given star
    note this requires default use of paths as given in perplexdata.py"""
    import glob
    sn = star_name.replace(' ', '')
    # print('searching for', output_parent_default + existing_dir + '*' + sn + '*')
    path = glob.glob(output_parent_default + existing_dir + '*' + sn + '*')
    # print('path', path)
    return path[0]  # should only be one but glob.glob returns list
