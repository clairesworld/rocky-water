import numpy as np
import parameters as p


def stellar_mantle(oxide_list, nH_star, core_eff):
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

