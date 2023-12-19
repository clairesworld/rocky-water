import numpy as np
import py.parameters as p


def mass_ratio_to_mol_ratio(m_i, m_j, s_i='', s_j=''):
    """ convert mass ratio i/j of oxides to mol ratio where s_i is a string of the oxide chemical formula """
    try:
        n_i = int(s_i[2])  # get oxide stoichiometry
    except ValueError:
        n_i = 1
    try:
        n_j = int(s_j[2])
    except ValueError:
        n_j = 1
    return (m_i / eval('p.M_' + s_i) / n_i) / (m_j / eval('p.M_' + s_j) / n_j)


NaTi_sol = (10 ** p.na_sol) / (10 ** p.ti_sol)
NaTi_bse = mass_ratio_to_mol_ratio(0.36, 0.201, 'Na2O', 'TiO2')
# print('NaTi_bse', NaTi_bse, 'NaTi_sol', NaTi_sol, 'ratio', NaTi_bse / NaTi_sol)


def get_stellar_percent_abundance(nH_star, which='moles', oxide_list=None):
    # convert list of nH_star into relative abundance in mol% (with respect to other elements in list)
    rel = [10 ** x for x in nH_star]

    if which == 'mass':  # n = m / M
        rel = [n * eval('p.M_' + ox[:2]) for (n, ox) in zip(rel, oxide_list)]  # will be renormalised
    elif which == 'moles':
        pass

    # renormalise
    print('-----------\nstellar composition (%)\n-----------')
    factor = 100 / sum(rel)
    rel = [x * factor for x in rel]
    [print("{:5.2f}%".format(x)) for x in rel]
    return rel


def mole_ratio_uncertainty(err_Q, err_R):
    """
    from Hinklel Young & Wheeler (2022), AJ, eqn. (20)
    returns error on mole ratio Q/R
    err_Q and err_R is half the total uncertainty (+/-) on [Q/H]_* and [R/H}_* in dex
    """
    return np.log(10) * np.sqrt(err_Q ** 2 + err_R ** 2)


# print('delta_Mg/Si',mole_ratio_uncertainty(0.07, 0.05))


def stellar_mantle_corelightelements(oxide_list, nH_star, core_eff, depletion_NaTi=None, core_Si_wtpt=None, verbose=True, **kwargs):
    """
    Convert stellar abundances directly to bulk mantle composition in terms of wt% oxides. Requires file parameters.py
    definining molar masses of each desired element in g/mol; e.g. M_Mg = 24.305

    Parameters
    ----------
    oxide_list : list
        Names of oxides to include, items are case-insensitive; e.g., 'mgo', 'Al2O3'
    nH_star : list
        Absolute log10(N_X/N_H) for the star; must be in same order as oxide_list
    core_eff : float
        Mole fraction of Fe in core (instead of FeO in mantle) with respect to total bulk planet Fe
    depletion_NaTi : float
        Depletion of Na/Ti in Earth BSE with respect to solar

    Returns
    -------
    wt_oxides_dict : dict
        Dictionary where keys are from oxide_list and entries are the mantle's proportion of that oxide in weight
        percent
    """

    if depletion_NaTi is None:
        depletion_NaTi = NaTi_bse / NaTi_sol

    # get molar masses
    try:
        M_oxides = [eval('p.M_' + ox) for ox in oxide_list]
    except Exception as e:
        print('oxide list', oxide_list)
        for ox in oxide_list:
            print(eval('p.M_' + ox))
        raise Exception('missing molar mass in parameters.py :', e)

    wt_oxides = [1]  # give denomenator 1 for now
    X_ratio_mol = [1]  # molar ratios of oxides
    mol_Fe_core = None  # for tracking moles of Fe in core wrt moles of bulk planet

    for ii, ox in enumerate(oxide_list):
        # print('calculating oxide', ox)
        if ii > 0:
            if ox[2] == '2':  # 2 moles cation per 1 mol oxide
                # print(ox)
                X_ratio_mol.append(0.5 * 10 ** nH_star[ii] / 10 ** nH_star[0])  # 2 mols Na per 1 mol Na2O e.g.
            else:
                X_ratio_mol.append(10 ** nH_star[ii] / 10 ** nH_star[0])  # cancel out H abundance

            # now do special mods
            if ox == 'FEO' or ox == 'FeO':
                # some molar percentage of Fe goes to core instead of to mantle FeO
                # this is the stage where core Fe is extracted from wt_oxides (being the bulk mantle composition)
                mol_Fe_core = X_ratio_mol[ii] * core_eff  # how many moles of Fe in core wrt bulk planet
                X_ratio_mol[ii] = X_ratio_mol[ii] * (1 - core_eff)  # update how many moles of Fe (as FeO) in mantle
                # print('mol Fe mantle', X_ratio_mol[ii])
                # print('mol_Fe_core', mol_Fe_core)

            if ox == 'NA2O' or ox == 'Na2O':
                # account for depletion: Na2O
                try:
                    idx_Ti = oxide_list.index('TiO2')
                    NaTi_star = 10 ** nH_star[ii] / 10 ** nH_star[idx_Ti]  # molar ratio
                    # print('Na/Ti_star', NaTi_star)
                    X_ratio_mol[ii] = depletion_NaTi * NaTi_star * X_ratio_mol[idx_Ti] * 2
                except ValueError:
                    raise NotImplementedError(
                        'ERROR: Can only define Na depletion with respect to Ti, must enter Ti first in oxide_list')

    # adjust according to core light elements
    if core_Si_wtpt is not None:
        # calculate moles of Si needed in core
        sife_molratio_core = (core_Si_wtpt / 100) / (1 - (core_Si_wtpt / 100)) * (p.M_Fe / p.M_Si)
        # print('core_Si_wtpt', core_Si_wtpt, 'wt%')
        # print('core_si_molefrac', sife_molratio_core)
        mol_Si_core = mol_Fe_core * sife_molratio_core

        # print('Si/Fe molar ratio core', sife_molratio_core)
        # print('mol_Si_core', mol_Si_core)
        # print('(p.M_Fe / p.M_Si)', (p.M_Fe / p.M_Si))

        # calculate how many moles must be extracted from mantle
        idx_Si = oxide_list.index('SiO2')
        if verbose:
            print('Mg/Si_bulk =', X_ratio_mol[oxide_list.index('MgO')] / X_ratio_mol[idx_Si])
        mol_Si_mantle_postcore = X_ratio_mol[idx_Si] - mol_Si_core
        X_ratio_mol[idx_Si] = mol_Si_mantle_postcore
        if verbose:
            print('Mg/Si_mantle =', X_ratio_mol[oxide_list.index('MgO')] / X_ratio_mol[idx_Si])

    # print('mol ratio oxides\n-----------')
    # for chem, val in zip(oxide_list, X_ratio_mol):
    #     try:
    #         print("{0:<7}".format(chem[:2]), "{:5.2f}".format(val*100))
    #     except TypeError:
    #         print('oxide_list', oxide_list, '\n wt_oxides', wt_oxides)

    # adjust weight percents of mantle oxides
    for ii, ox in enumerate(oxide_list):
        if ii > 0:
            # print('ii, ox', ii, ox)
            M_ratio = M_oxides[ii] / M_oxides[0]
            # print('M_ratio', M_ratio)
            m_ratio = X_ratio_mol[ii] * M_ratio  # convert from mole ratio to mass ratio assuming all metals in oxides
            # print('m_ratio', m_ratio)
            wt_oxides.append(m_ratio)
            # print('-')

    mass_core_check = (mol_Fe_core * p.M_Fe / M_oxides[0]) + (mol_Si_core * p.M_Si / M_oxides[0])
    # print('mass core', mass_core_check)
    mass_mantle_check = sum(wt_oxides)
    # print('mass mantle', mass_mantle_check)
    # print('CMF', mass_core_check / (mass_mantle_check + mass_core_check))

    # now have ratios in wt% 1:X1:X2:X3
    wt_oxides = np.array(wt_oxides)
    wt_oxides = wt_oxides / sum(wt_oxides) * 100  # normalise so total wt% is 100

    if verbose:
        print('wt.% oxides\n-----------')
        for chem, val in zip(oxide_list, wt_oxides):
            try:
                print("{0:<7}".format(chem), "{:5.2f}%".format(val))
            except TypeError:
                print('oxide_list', oxide_list, '\n wt_oxides', wt_oxides)

    # put in dict
    wt_oxides_dict = {k: v for k, v in zip(oxide_list, wt_oxides)}
    return wt_oxides_dict


def core_mass_fraction(wt_oxides, core_eff, core_light_wtpt=None):
    """
    Calculates core mass fraction given bulk FeO in mantle and the wt fraction of total FeO that ends up in core.
    Doesn't require knowing the planet mass
    core_eff = n_Fe_core / (n_FeO_mantle + n_Fe_core) = m_Fe_core / (m_Fe_mantle + m_Fe_core) as n_Fe = n_FeO in mtl
    core_light_wt: dict of core light element compositon in wt%
    """
    # todo: still might want to triple check this is consistent with stellar composition constraints

    if core_light_wtpt is None:
        core_light_wtpt = {'none': 0}

    try:
        x_Fe_mm = wt_oxides['FeO'] * (p.M_Fe / p.M_FeO)  # mass fraction Fe in mantle wrt total mantle mass (excluding O)
        # print('x_Fe_m wrt mantle', x_Fe_mm)
    except KeyError:
        x_Fe_mm = 0

    if core_eff == 1:
        # TODO
        x_core_m = None  # mass fraction Fe in core wrt total mantle mass
    else:
        x_Fe_cm = x_Fe_mm * core_eff / (1 - core_eff)  # mass fraction Fe in core wrt total mantle mass
        x_core_m = x_Fe_cm * (100 + sum(core_light_wtpt.values())) / 100  # mass fraction total core wrt total mantle mass
        # print('x_Fe_c wrt mantle', x_Fe_cm)

    mantle_mass_fraction = 100 / (100 + x_core_m)
    # print('mantle mass fraction', mantle_mass_fraction)

    # scale by mass_mtl/M_p --> should be smaller than x_Fe_cm
    x_core = x_core_m * mantle_mass_fraction
    # print('x_core wrt planet', x_core)
    CMF = x_core / 100

    print('\ncore mass fraction = {:5.3f}\n'.format(CMF))
    return CMF

# core_mass_fraction(wt_oxides, core_eff, core_light_wtpt=None)


def stellar_mantle(oxide_list, nH_star, core_eff, depletion_NaTi=None, **kwargs):
    """
    Convert stellar abundances directly to bulk mantle composition in terms of wt% oxides. Requires file parameters.py
    definining molar masses of each desired element in g/mol; e.g. M_Mg = 24.305

    Parameters
    ----------
    oxide_list : list
        Names of oxides to include, items are case-insensitive; e.g., 'mgo', 'Al2O3'
    nH_star : list
        Absolute log10(N_X/N_H) for the star; must be in same order as oxide_list
    core_eff : float
        Mole fraction of Fe in core (instead of FeO in mantle) with respect to total bulk planet Fe
    depletion_NaTi : float
        Depletion of Na/Ti in Earth BSE with respect to solar

    Returns
    -------
    wt_oxides_dict : dict
        Dictionary where keys are from oxide_list and entries are the mantle's proportion of that oxide in weight
        percent
    """

    if depletion_NaTi is None:
        depletion_NaTi = NaTi_bse / NaTi_sol

    # get molar masses
    try:
        M_oxides = [eval('p.M_' + ox) for ox in oxide_list]
    except Exception as e:
        print('oxide list', oxide_list)
        for ox in oxide_list:
            print(eval('p.M_' + ox))
        raise Exception('missing molar mass in parameters.py :', e)

    wt_oxides = [1]  # give denomenator 1 for now
    X_ratio_mol = [1]  # molar ratios of oxides
    for ii, ox in enumerate(oxide_list):
        # print('calculating oxide', ox)
        if ii > 0:
            if ox[2] == '2':  # 2 moles cation per 1 mol oxide
                # print(ox)
                X_ratio_mol.append(0.5 * 10 ** nH_star[ii] / 10 ** nH_star[0])  # 2 mols Na per 1 mol Na2O e.g.
            else:
                X_ratio_mol.append(10 ** nH_star[ii] / 10 ** nH_star[0])  # cancel out H abundance

            # now do special mods
            if ox == 'FEO' or ox == 'FeO':
                # some molar percentage of Fe goes to core instead of to mantle FeO
                # this is the stage where core Fe is extracted from wt_oxides (being the bulk mantle composition)
                X_ratio_mol[ii] = X_ratio_mol[ii] * (1 - core_eff)

            if ox == 'NA2O' or ox == 'Na2O':
                # account for depletion: Na2O
                try:
                    idx_Ti = oxide_list.index('TiO2')
                    NaTi_star = 10 ** nH_star[ii] / 10 ** nH_star[idx_Ti]  # molar ratio
                    # print('Na/Ti_star', NaTi_star)
                    X_ratio_mol[ii] = depletion_NaTi * NaTi_star * X_ratio_mol[idx_Ti] * 2
                except ValueError:
                    raise NotImplementedError(
                        'ERROR: Can only define Na depletion with respect to Ti, must enter Ti first in oxide_list')

            # convert from mole ratio to mass ratio assuming all metals in oxides
            M_ratio = M_oxides[ii] / M_oxides[0]
            m_ratio = X_ratio_mol[ii] * M_ratio
            wt_oxides.append(m_ratio)

    # now have ratios in wt% 1:X1:X2:X3
    wt_oxides = np.array(wt_oxides)
    # print('len wt_oxides', len(wt_oxides))
    wt_oxides = wt_oxides / sum(wt_oxides) * 100  # normalise so total wt% is 100

    print('wt.% oxides\n-----------')
    for chem, val in zip(oxide_list, wt_oxides):
        try:
            print("{0:<7}".format(chem), "{:5.2f}%".format(val))
        except TypeError:
            print('oxide_list', oxide_list, '\n wt_oxides', wt_oxides)

    # put in dict
    wt_oxides_dict = {k: v for k, v in zip(oxide_list, wt_oxides)}
    # print('Na depletion factor:', (mass_ratio_to_mol_ratio(wt_oxides_dict['Na2O'], wt_oxides_dict['TiO2'], 'Na2O', 'TiO2')) / NaTi_star, '(want', depletion_NaTi, ')')
    # print('Na/Ti planet', mass_ratio_to_mol_ratio(wt_oxides_dict['Na2O'], wt_oxides_dict['TiO2'], 'Na2O', 'TiO2'))
    return wt_oxides_dict

def o2_molar_ratio(n_FeO, X_ferric=0.05, **kwargs):
    """ Given a ratio of Fe3+/total Fe, and the number of moles of Fe(II and III)O (=n_FeO), what is the molar abundance of O2?
    4FeO + O2 <--> 2Fe2O3  ===  FeO + 1/4O2 <--> FeO1.5
    X_ferric = X_FeO1.5 / (X_FeO + X_FeO1.5)
    """

    n_ferric = X_ferric * n_FeO
    n_ferrous = n_FeO - n_ferric
    n_O2 = X_ferric / 4 * n_FeO  # 2 mol Fe2O3 for every O2 and 2 mol Fe3+ for every mol Fe2O3

    # print('n_Fe2+', n_ferrous, 'n_O2', n_O2, 'n_Fe3+', n_ferric, 'r', X_ferric)
    return n_O2


def test_ferric_ratio_from_O2(m_FeOstar, m_O2):
    """ 4FeO + O2 <--> 2Fe2O3, where O2 is limiting """

    n_FeO = m_FeOstar / p.M_FeO
    n_O2 = m_O2 / p.M_O2
    # n_ferrous_eq = n_O2 * 4  # number of mol Fe(2+)O being used = 2
    n_Fe2O3_eq = n_O2 * 2  # number of mol Fe(3+)2O3 taking up oxygen = 1
    n_ferric_eq = n_Fe2O3_eq * 2  # number of mol Fe(3+) from above - i.e. 2 mol Fe(3+) for every mol Fe2O3
    # n_ferrous_extra = n_FeO - n_ferric_eq  # number of mol Fe2+ left = 8
    X_ferric = n_ferric_eq / n_FeO  # molar ratio of Fe(3+) to total Fe
    # print('n_Fe2+', n_ferrous_extra, 'n_O2', n_O2, 'n_Fe3+', n_ferric_eq, 'r', X_ferric)
    return X_ferric


# test_ferric_ratio_from_O2(n_FeO=10, n_O2=0.12)


def insert_o2_from_wt_comp(wt_oxides_dict, **kwargs):
    """ Get O2 abundance in wt% given existing oxide composition in wt% """

    # convert to mole fraction
    molar_dict = {}
    for key, m in wt_oxides_dict.items():
        if key != 'O2':
            M = eval('p.M_' + key)
            molar_dict[key] = m / M

    # enter O2 molar ratio from FeO
    molar_dict['O2'] = o2_molar_ratio(molar_dict['FeO'], **kwargs)

    # # print
    # print('mol.% oxides\n-----------')
    # factor = 100 / sum(molar_dict.values())
    # for k in molar_dict:
    #     print("{0:<7}".format(k), "{:5.2f}%".format(molar_dict[k] * factor))

    # reconvert to wt percent
    new_wt_dict = {}
    for key, n in molar_dict.items():
        M = eval('p.M_' + key)
        new_wt_dict[key] = n * M  # convert to mole fraction

    # renormalise
    print('-----------\nwt.% oxides with O2 added\n-----------')
    factor = 100 / sum(new_wt_dict.values())
    for k in new_wt_dict:
        new_wt_dict[k] = new_wt_dict[k] * factor
        print("{0:<7}".format(k), "{:5.2f}%".format(new_wt_dict[k]))
    print('-----------')
    return new_wt_dict


def insert_fe2o3(wt_oxides_dict, X_ferric, mass_tot=100, verbose=True):
    m_FeOstar = wt_oxides_dict['FeO']
    m_FeO, m_Fe2O3 = partition_FeOstar(m_FeOstar, X_ferric)

    # renormalise
    wt_oxides_dict['FeO'] = m_FeO
    wt_oxides_dict['Fe2O3'] = m_Fe2O3

    if verbose:
        print('-----------\nwt.% oxides with Fe2O3 added\n-----------')
    factor = mass_tot / sum(wt_oxides_dict.values())
    for k in wt_oxides_dict:
        wt_oxides_dict[k] = wt_oxides_dict[k] * factor
        if verbose:
            print("{0:<7}".format(k), "{:5.2f}%".format(wt_oxides_dict[k]))
    print('-----------')
    return wt_oxides_dict


def partition_FeOstar(m_FeOstar, X_ferric):
    # A = (2 * (1 - X_ferric) * p.M_FeO) / (X_ferric * p.M_Fe2O3)
    # m_FeO = A / (1 + A)
    # m_Fe2O3 = m_FeOstar - m_FeO

    m_FeO = 2 * p.M_FeO * m_FeOstar * (1 - X_ferric) / (p.M_Fe2O3 * X_ferric - 2 * p.M_FeO * X_ferric + 2 * p.M_FeO)
    m_Fe2O3 = p.M_Fe2O3 * X_ferric * m_FeOstar / (p.M_Fe2O3 * X_ferric - 2 * p.M_FeO * X_ferric + 2 * p.M_FeO)

    # print('m_FeO', m_FeO, 'g')
    # print('m_Fe2O3', m_Fe2O3, 'g')
    return m_FeO, m_Fe2O3


def test_ferric_ratio(m_FeO, m_Fe2O3):
    n_Fe2O3 = m_Fe2O3 / p.M_Fe2O3
    n_FeIII = 2 * n_Fe2O3
    n_FeO = m_FeO / p.M_FeO
    n_FeII = n_FeO
    Xfer = n_FeIII / n_FeII
    return Xfer


def get_element_ratio(ratio_str, wt_oxides):
    """

    Parameters
    ----------
    ratio_str : string of elements e.g. 'Mg/Si'
    wt_oxides : dict of wt oxide composition

    Returns
    -------
    float of element ratio
    """

    # turn wt oxides keys into elements
    k_els = [k[:2] for k in wt_oxides.keys()]

    els = ratio_str.split('/')  # should get 2
    n = np.zeros(2)

    for ii, el in enumerate(els):
        idx = k_els.index(el)
        ox = list(wt_oxides.keys())[idx]
        M = eval('p.M_' + ox)

        # get moles of oxide
        n_ox = wt_oxides[ox] / M

        # parse stoichiometry to get moles of element - if 3rd char is a digit
        try:
            a = int(ox[2])
        except ValueError:
            a = 1
        n[ii] = a * n_ox

        # check for Fe2O3
        if (el == 'Fe') and ('Fe2O3' in wt_oxides):
            idx = list(wt_oxides.keys()).index('Fe2O3')
            n_Fe2O3 = wt_oxides['Fe2O3'] / p.M_Fe2O3 * 2
            n[ii] = n[ii] + n_Fe2O3

    # print(ratio_str, n[0] / n[1])
    return n[0] / n[1]


def total_refractory_O(wt_oxides_dict):
    # renormalise
    factor = 100 / sum(wt_oxides_dict.values())
    for k in wt_oxides_dict:
        wt_oxides_dict[k] = wt_oxides_dict[k] * factor

    # get mass of refractory oxygen in bulk silicate planet (out of 100 g)
    m_O_tot = 0
    # if 'Fe2O3' in wt_oxides_dict:  # this has no O2..?
    for k, mass in wt_oxides_dict.items():
        cat, n_O = k.split('O')
        # print('\n', k, 'cat', cat, 'n_o', n_O)

        # get oxygen mass
        if not n_O:
            n_O = 1  # 1 oxygen per molecule
        else:
            n_O = int(n_O)
        m_O = n_O * p.M_O
        # print('m_O', m_O)

        # get cation mass
        if not cat:  # just O2
            m_cat = 0
        else:
            try:
                n_cat = int(cat[2])
            except (IndexError, TypeError):
                n_cat = 1
            m_cat = n_cat * eval('p.M_' + cat[:2])
        # print('m_cat', m_cat)

        # fraction of mass of molecule that is O
        frac_m_O = m_O / (m_O + m_cat)
        # print('frac_m_O', frac_m_O)

        # accumulate total
        m_O_tot += frac_m_O * mass

    return m_O_tot


def Fecore_mass_fraction(wt_oxides, core_eff):
    """
    Calculates core mass fraction given bulk FeO in mantle and the wt fraction of total FeO that ends up in core.
    Doesn't require knowing the planet mass
    core_eff = n_Fe_core / (n_FeO_mantle + n_Fe_core) = m_Fe_core / (m_Fe_mantle + m_Fe_core) as n_Fe = n_FeO in mtl
    core_light_wt: dict of core light element compositon in wt%
    """
    # todo: still might want to triple check this is consistent with stellar composition constraints

    try:
        x_Fe_mm = wt_oxides['FeO'] * (p.M_Fe / p.M_FeO)  # mass fraction Fe in mantle wrt total mantle mass
    except KeyError:
        raise Exception(
            'ERROR: cmf is undefined for a given molar core efficiency if no FeO in bulk mantle. try fixing input cmf instead.')

    if core_eff == 1:
        # TODO
        x_Fe_cm = 1 - x_Fe_mm  # mass fraction Fe in core wrt total mantle mass
    else:
        x_Fe_cm = -x_Fe_mm * core_eff / (core_eff - 1)  # mass fraction Fe in core wrt total mantle mass
    mantle_mass_fraction = 100 / (100 + x_Fe_cm)

    # scale by mass_mtl/M_p --> should be smaller than x_Fe_cm
    x_Fe_c = x_Fe_cm * mantle_mass_fraction
    # print('x_Fe_c wrt planet', x_Fe_c)
    CMF = x_Fe_c / 100

    Fe_mass_fraction = x_Fe_c + x_Fe_mm * mantle_mass_fraction  # store bulk planet Fe fraction for posterity
    # print('x_Fe_pl', x_Fe_pl)

    print('\ncore mass fraction = {:5.3f}\n'.format(CMF))
    return CMF


def core_eff_from_cmf(CMF, M_p, wt_oxides, core_light_wt=None):
    M_c = CMF * M_p  # core mass in kg
    M_m = M_p - M_c  # mantle mass in kg
    n_Fe_mtl = wt_oxides['FeO'] / 100 * M_m / p.M_FeO  # n_Fe = n_FeO
    n_Fe_core = M_c / p.M_Fe  # pure Fe core
    core_eff = n_Fe_core / (n_Fe_core + n_Fe_mtl)
    return core_eff


# wt_oxides_MD95 = {'SiO2': 45.0, 'MgO': 37.8, 'CaO': 3.55, 'Al2O3': 4.45, 'FeO': 8.05}
# wt_oxides1 = update_MgSi(1.6, wt_oxides_MD95)
#
# wt_oxides1 = insert_o2_from_wt_comp(wt_oxides1, X_ferric=0.03)
# wt_oxides2 = update_MgSi(0.8, wt_oxides1)
#
# print('Mg/Si = 0.8, O (wt%) =', total_refractory_O(wt_oxides2))
# print('Mg/Si = 1.5, O (wt%) =', total_refractory_O(wt_oxides1))
#
#
# M_Mg2SiO4 = (2 * p.M_Mg) + (p.M_Si) + (4 * p.M_O)
# M_MgSiO3 = (p.M_Mg) + (p.M_Si) + (3 * p.M_O)
# print('M Mg2SiO4', M_Mg2SiO4, 'g', ', frac O', (4 * p.M_O) / M_Mg2SiO4)
# print('M MgSiO3', M_MgSiO3, 'g', ', frac O', (3 * p.M_O) / M_MgSiO3)

# dmm = {'SiO2': 44.71, 'Al2O3':3.98, 'FeO': 8.008 + 0.191, 'MgO': 38.73, 'CaO': 3.17}
# test_ferric_ratio(8.008, 0.191)
# insert_fe2o3(dmm, 0.03, mass_tot=sum(dmm.values()))
#
#
# import sympy as sym
#
# m_FeOstar, X_ferric, m_FeO, m_Fe2O3, M_FeO, M_Fe2O3 = sym.symbols('m_FeOstar,X_ferric,m_FeO,m_Fe2O3,M_FeO,M_Fe2O3')
# eq1 = sym.Eq(m_FeO+m_Fe2O3,m_FeOstar)
# eq2 = sym.Eq((2*m_Fe2O3/M_Fe2O3) / (m_FeO / M_FeO + 2*m_Fe2O3/M_Fe2O3),X_ferric)
# result = sym.solve([eq1,eq2],(m_FeO,m_Fe2O3))
# print(result)


