import numpy as np
import parameters as p

def bulk_composition(oxide_list, nH_star, core_eff):
    """ nH_star is list of the absolute log(N_X/N_H) for the star, in same order
         would need to convert from wrt solar if necessary
         by convention base element (denomenator for oxide wt ratios) is first in list
         core_efficiency is fraction of moles Fe that go into core instead of mantle FeO wrt total moles Fe
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
                # poss confusing but this is the stage where core Fe is extracted from wt_oxides (which refers to mtl)
                X_ratio_mol = X_ratio_mol * (1 - core_eff)

            # convert from mole ratio to mass ratio assuming all metals in oxides
            M_ratio = M_oxides[ii] / M_oxides[0]
            m_ratio = X_ratio_mol * M_ratio
            wt_oxides.append(m_ratio)

    # now have ratios in wt% 1:X1:X2:X3... normalise so total wt% is 100
    wt_oxides = np.array(wt_oxides)
    wt_oxides = wt_oxides / sum(wt_oxides) * 100

    print('wt.% oxides\n-----------')
    for chem, val in zip(oxide_list, wt_oxides):
        print("{0:<7}".format(chem), "{:5.2f}%".format(val))

    # put in dict
    wt_oxides_dict = {k: v for k, v in zip(oxide_list, wt_oxides)}
    return wt_oxides_dict


# from parameters import ca_sol, fe_sol, al_sol, mg_sol, si_sol, na_sol
# oxide_list_default = oxide_list_default = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO']  # 'NaO
# solar = {'ca_sol': ca_sol, 'fe_sol': fe_sol, 'al_sol': al_sol, 'mg_sol': mg_sol, 'si_sol': si_sol, 'na_sol': na_sol}
# nH_sun = [solar[ox[:2].lower() + '_sol'] for ox in oxide_list_default]
#
# FeH_range = (-0.5 + fe_sol, 0.5 + fe_sol)
# for f in FeH_range:
#     nH_sun[-1] = f
#     print(bulk_composition(oxide_list_default, nH_sun, 0.88))
