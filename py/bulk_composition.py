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


NaTi_sol = 10 ** p.na_sol / 10 ** p.ti_sol
NaTi_bse = mass_ratio_to_mol_ratio(0.36, 0.2, 'Na2O', 'TiO2')


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
    X_ratio_mol = [1]
    for ii, ox in enumerate(oxide_list):
        # print('calculating oxide', ox)
        if ii > 0:
            if ox == 'AL2O3' or ox == 'Al2O3':
                X_ratio_mol.append(0.5 * 10 ** nH_star[ii] / 10 ** nH_star[0])  # 2 mols Al per 1 mol Al2O3
            elif ox == 'NA2O' or ox == 'Na2O':
                X_ratio_mol.append(0.5 * 10 ** nH_star[ii] / 10 ** nH_star[0])  # 2 mols Na per 1 mol Na2O3
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
                    NaTi_star = X_ratio_mol[ii] / X_ratio_mol[idx_Ti]
                    X_ratio_mol[ii] = depletion_NaTi * NaTi_star * X_ratio_mol[idx_Ti]
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
        print("{0:<7}".format(chem), "{:5.2f}%".format(val))

    # put in dict
    wt_oxides_dict = {k: v for k, v in zip(oxide_list, wt_oxides)}
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


# o2_molar_ratio(10, 0.05)


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
    print('wt.% oxides with O2 added\n-----------')
    factor = 100 / sum(new_wt_dict.values())
    for k in new_wt_dict:
        new_wt_dict[k] = new_wt_dict[k] * factor
        print("{0:<7}".format(k), "{:5.2f}%".format(new_wt_dict[k]))
    return new_wt_dict


def insert_fe2o3(wt_oxides_dict, X_ferric, mass_tot=100):
    m_FeOstar = wt_oxides_dict['FeO']
    m_FeO, m_Fe2O3 = partition_FeOstar(m_FeOstar, X_ferric)

    # renormalise
    wt_oxides_dict['FeO'] = m_FeO
    wt_oxides_dict['Fe2O3'] = m_Fe2O3

    print('wt.% oxides with Fe2O3 added\n-----------')
    factor = mass_tot / sum(wt_oxides_dict.values())
    for k in wt_oxides_dict:
        wt_oxides_dict[k] = wt_oxides_dict[k] * factor
        print("{0:<7}".format(k), "{:5.2f}%".format(wt_oxides_dict[k]))

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


