import numpy as np
import parameters as p


def stellar_mantle(oxide_list, nH_star, core_eff, **kwargs):
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

    Returns
    -------
    wt_oxides_dict : dict
        Dictionary where keys are from oxide_list and entries are the mantle's proportion of that oxide in weight
        percent
    """

    # get molar masses
    try:
        M_oxides = [eval('p.M_' + ox) for ox in oxide_list]
    except Exception as e:
        print('oxide list', oxide_list)
        for ox in oxide_list:
            print(eval('p.M_' + ox))
        raise Exception('missing molar mass in parameters.py :', e)

    wt_oxides = [1]  # give denomenator 1 for now
    for ii, ox in enumerate(oxide_list):
        if ii > 0:
            if ox == 'AL2O3' or ox == 'Al2O3':
                X_ratio_mol = 0.5 * 10 ** nH_star[ii] / 10 ** nH_star[0]  # 2 mols Al per 1 mol Al2O3
            elif ox == 'NA2O' or ox == 'Na2O':
                X_ratio_mol = 0.5 * 10 ** nH_star[ii] / 10 ** nH_star[0]  # 2 mols Na per 1 mol Na2O3
            else:
                X_ratio_mol = 10 ** nH_star[ii] / 10 ** nH_star[0]  # cancel out H abundance

            if ox == 'FEO' or ox == 'FeO':
                # some molar percentage of Fe goes to core instead of to mantle FeO
                # this is the stage where core Fe is extracted from wt_oxides (being the bulk mantle composition)
                X_ratio_mol = X_ratio_mol * (1 - core_eff)
                idx_FeO = ii

            if ox == 'O2':
                # get oxygen abundance from Fe3+/Fe2+ ratio
                try:
                    X_FeO = X_ratio_mol[idx_FeO]
                except:
                    raise Exception("Can't get O2 abundance without FeO in system")
                X_ratio_mol = o2_molar_ratio(X_FeO, **kwargs)

            # convert from mole ratio to mass ratio assuming all metals in oxides
            M_ratio = M_oxides[ii] / M_oxides[0]
            m_ratio = X_ratio_mol * M_ratio
            wt_oxides.append(m_ratio)

    # now have ratios in wt% 1:X1:X2:X3
    wt_oxides = np.array(wt_oxides)
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

    print('n_Fe2+', n_ferrous, 'n_O2', n_O2, 'n_Fe3+', n_ferric, 'r', X_ferric)
    return n_O2
# o2_molar_ratio(10, 0.05)


def ferric_ratio(n_FeO, n_O2):
    """ 4FeO + O2 <--> 2Fe2O3, where O2 is limiting """

    n_ferrous_eq = n_O2 * 4  # number of mol Fe(2+)O being used = 2
    n_Fe2O3_eq = n_O2 * 2  # number of mol Fe(3+)2O3 taking up oxygen = 1
    n_ferric_eq = n_Fe2O3_eq * 2  # number of mol Fe(3+) from above - i.e. 2 mol Fe(3+) for every mol Fe2O3
    n_ferrous_extra = n_FeO - n_ferric_eq  # number of mol Fe2+ left = 8
    X_ferric = n_ferric_eq / n_FeO  # molar ratio of Fe(3+) to total Fe
    print('n_Fe2+', n_ferrous_extra, 'n_O2', n_O2, 'n_Fe3+', n_ferric_eq, 'r', X_ferric)
# ferric_ratio(n_FeO=10, n_O2=0.12)


def o2_from_wt_comp(wt_oxides_dict, **kwargs):
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
    print('wt.% oxides\n-----------')
    factor = 100 / sum(new_wt_dict.values())
    for k in new_wt_dict:
        new_wt_dict[k] = new_wt_dict[k] * factor
        print("{0:<7}".format(k), "{:5.2f}%".format(new_wt_dict[k]))
    return new_wt_dict
