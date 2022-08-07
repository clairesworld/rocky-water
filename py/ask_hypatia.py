import numpy as np
import requests
import random
import pickle as pkl
import os

# Enter API key generated from https://www.hypatiacatalog.com/api
key = '39728dfb7021ad38c12aa84d66a7f6fd'


def retrieve_star_names(exo_hosts=True, API_KEY=key, writeto='host_names.txt', exclude_blank=False):
    """
    Get list of star names from the Hypatia catalogue. Note: this can take a few minutes, only need to run if names
    file does not exist.

    Parameters
    ----------
    exo_hosts : bool
        If True, return only names of exoplanet host stars
    API_KEY : str
        API key generated from https://www.hypatiacatalog.com/api
    writeto : str or None
        Path + filename to write list of stars, if None then write nothing
    exclude_blank : bool
        [todo] Exclude star names with no measured mg, si, fe etc.

    Returns
    -------
    names : list
        List of star names
    """
    catalog = requests.get("https://hypatiacatalog.com/hypatia/api/v2/data", auth=(API_KEY, "api_token"))
    names_all = []
    # print(catalog.json()['values'])
    for row in catalog.json()['values']:  # loop through list of dicts
        names_all.append(row['name'])  # add to list - todo there is a faster way to do this but cba

    if not exo_hosts and not exclude_blank:
        names = names_all  # simply return all star names in the catalogue

    else:
        params = {}
        if exo_hosts:
            params['name'] = names_all
            address = "https://hypatiacatalog.com/hypatia/api/v2/star"
        # if (not not exclude_blank):
        #     address = "https://hypatiacatalog.com/hypatia/api/v2/composition"
        #     params['element'] = []
        #     for el in exclude_blank:
        #         params['element'].append(el)

        entry = requests.get(address, auth=(API_KEY, "api_token"),
                             params=params)

        names = []  # initialise blank list
        for ii, r in enumerate(entry.json()):  # loop through list of dicts
            if not exo_hosts or (exo_hosts and len(r['planets']) > 0):
                names.append(names_all[ii])  # only add if more than 0 planets
            # if not exclude_blank or (or not not exclude_blank and 'mg' in exclude_blank):
            #     # only add if not missing this measurement - todo

    if writeto is not None:
        with open(writeto, 'w') as filehandle:
            filehandle.writelines("%s\n" % a for a in names)  # write names to file line by line
    return names


def random_star(n=1, names_file='host_names.txt', **kwargs):
    """
    Get random star(s) from an input list of star names - N.B. there are ~1331 entries for exoplanet hosts and ~8304
    entries for all stars in the Hypatia catalogue.

    Parameters
    ----------
    n : int
        Number of random stars to return, if -1 then return all names. Note if n = len(names) then all stars will be
        returned but in random order.
    names_file : str or None
        Path + filename to read from; i.e., generated from retrieve_star_names(), if None then call this function
    kwargs : dict
        Keyword arguments to pass to retrieve_star_names() if applicable

    Returns
    -------
    sample_names : list
        List of star names
    """
    from inspect import getsourcefile
    from os.path import abspath
    wd = abspath(getsourcefile(lambda: 0))  # get absolute path of this file

    if names_file is not None:
        with open('/' + os.path.join(*(wd.split(os.path.sep)[:-1])) + '/' + names_file, 'r') as filehandle:
            names = [a.rstrip() for a in filehandle.readlines()]  # read file line by line into list
    else:
        names = retrieve_star_names(**kwargs)

    if n == -1:  # return all stars
        return names  # in same order as file

    idx = random.sample(range(0, len(names)), n)  # get random indices
    sample_names = []
    for i in idx:
        sample_names.append(names[i])  # todo faster way to do this but cba
    return sample_names


def star_composition(oxide_list=None, star='sun', API_KEY=key, use_local_composition=False,
                     output_parent_path=None, get_hypatia_min=None, get_hypatia_max=None, verbose=True,
                     **kwargs):
    """
    Retrieve the elemental abundances of an individual star. Requires file parameters.py definining solar abundances for
    each element X in log10(N_X/N_H); e.g., mg_sol = 7.54 - 12

    Parameters
    ----------
    oxide_list : list
        Names of elements to retrieve, items are case-insensitive and can be oxide; e.g., 'Mg', 'mg', or 'MgO'
    star : str
        Name of the star in the same format as Hypatia Catalogue including spaces; e.g. "HIP 12345"
    API_KEY : str
        API key generated from https://www.hypatiacatalog.com/api
    use_local_composition : bool
        If True, load nH_star from a pickled object or textfile 'nH_star.txt' (one abundance ber line), rather than
        query Hypatia
    output_parent_path : str or Path object
        Parent directory in which to search for and load pickled object; only used if use_local_composition is True
    get_hypatia_min : list
        Names of elements for which the minimum abundance should be retrieived, given star exists in >1 catalogue
    get_hypatia_max : list
        Names of elements for which the maximum abundance should be retrieived, given star exists in >1 catalogue
    verbose : bool
        More detailed explanation in output if True

    Returns
    _______
    nH_star : list
        Absolute stellar abundances for each item in oxide_list X, as log10(N_X/N_H)
    """

    import parameters as p
    if get_hypatia_min is None:  # elements for which you want the minimum value
        get_hypatia_min = []
    if get_hypatia_max is None:  # elements for which you want the maximum value
        get_hypatia_max = []

    els = []
    for ox in oxide_list:
        els.append(ox[:2].lower())

    def do_local():
        try:
            path = find_existing_directory(star, output_parent=output_parent_path, **kwargs)
            try:
                nH_star = np.loadtxt(path + '/nH_star.txt')
            except OSError:
                # no data
                return None
            except FileNotFoundError:  # try pickled file
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
                raise Exception('No entry found in Hypatia Catalogue:', star)
            # print('loaded json', entry.json())
            nH_star = []
            for ii, el in enumerate(els):
                # absolute not working for some reason so get difference from solar via lodders norm
                try:
                    if (el not in get_hypatia_max) and (el not in get_hypatia_min):
                        nH = entry.json()[ii]['mean']
                        if verbose:
                            print('retrieving from Hypatia...', entry.json()[ii]['element'], 'mean =',
                                  entry.json()[ii]['mean'])
                    else:
                        all_values = [dd['value'] for dd in entry.json()[ii]['all_values']]  # this is a list of dicts
                        print('retieving from Hypatia...', entry.json()[ii]['element'], 'all values =', all_values)
                        if (el in get_hypatia_min) and (el in get_hypatia_max):
                            raise NotImplementedError(
                                'Hypatia query: cannot retrive both get min and max simultaneously for', el)
                        elif el in get_hypatia_min:
                            nH = min(all_values)
                        elif el in get_hypatia_max:
                            nH = max(all_values)
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


def find_existing_directory(star_name, existing_dir='hypatia0,1M_1600K_88Fe/', output_parent=None):
    """
    Retrieve directory of already-existing run for a given star. Best with default use of paths as given in
    perplexdata.py; i.e., output_path_default / existing_dir / *star_name*

    Parameters
    ----------
    star_name : str
        Name/ID of star as retrieved from Hypatia Catalogue
    existing_dir : str or Path
        Directory containing existing runs
    output_parent : str or None
        Shared parent directory of existing_dir and

    Returns
    -------
    path[0] : str
        Directory containing output data for desired existing run
    """
    from perplexdata import output_parent_default
    import glob

    if output_parent is None:
        output_parent = output_parent_default
    sn = star_name.replace(' ', '')
    # print('searching for', output_parent_default + existing_dir + '*' + sn + '*')
    try:
        path = glob.glob(output_parent + existing_dir + '*' + sn + '*')
        return path[0]  # should only find one item; glob.glob returns list
    except IndexError as e:
        # try one more directory up
        try:
            output_parent = '/' + os.path.join(*(output_parent.split(os.path.sep)[:-2])) + '/'
            path = glob.glob(output_parent + existing_dir + '*' + sn + '*')
            return path[0]  # should only be one but glob.glob returns list
        except IndexError as e:
            print('path', path)
            print('output_parent', output_parent)
            print('output_parent + existing', output_parent + existing_dir)
            raise e
