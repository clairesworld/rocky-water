import numpy as np
import bulk_composition as bulk
import parameters as p
from parameters import ca_sol, fe_sol, al_sol, mg_sol, si_sol, na_sol


def get_mgsi(wt_oxides):
    n_MgO = wt_oxides['MgO'] / p.M_MgO  # convert to moles
    n_SiO2 = wt_oxides['SiO2'] / p.M_SiO2
    mgsi = n_MgO / n_SiO2  # n_Mg = n_MgO etc
    return mgsi

M_p = 100  # assume 100 g
CMF = 0.325
core_eff = 0.84
core_Si_mol = 0.0646  # tronnes 2019 mole fraction of Si on core
# core_Si_wt = 0.036
core_Si_wt = 0.1
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO']
wt_oxides_PM = {'SiO2': 45, 'MgO': 37.8, 'CaO': 3.55, 'Al2O3': 4.45, 'FeO': 8.05}  # Tronnes 2019 primitive mantle

solar = {'ca_sol': ca_sol, 'fe_sol': fe_sol, 'al_sol': al_sol, 'mg_sol': mg_sol, 'si_sol': si_sol, 'na_sol': na_sol}
nH_star = [solar[ox[:2].lower() + '_sol'] for ox in oxide_list if ox != 'O2']
wt_oxides_sol = bulk.stellar_mantle(oxide_list, nH_star, core_eff)

""" calculate core mass fraction, mantle wt% composition given a mole fraction of Si in core """

wt_oxides = wt_oxides_sol

print('earth', get_mgsi(wt_oxides_PM))

print('Mg/Si before core partitioning', get_mgsi(wt_oxides))  # bulk.mass_ratio_to_mol_ratio(wt_oxides['MgO'], wt_oxides['SiO2'], s_i='MgO', s_j='SiO2'))

if CMF is None:
    # calculate CMF

    try:
        x_Fe_mm = wt_oxides['FeO'] * (p.M_Fe / p.M_FeO)  # mass fraction Fe in mantle wrt total mantle mass
    except KeyError:
        raise Exception(
            'ERROR: cmf is undefined for a given molar core efficiency if no FeO in bulk mantle. try fixing input cmf instead.')

    if core_eff == 1:
        x_Fe_cm = 1 - x_Fe_mm  # mass fraction Fe in core wrt total mantle mass
    else:
        x_Fe_cm = -x_Fe_mm * core_eff / (core_eff - 1)  # mass fraction Fe in core wrt total mantle mass
    mantle_mass_fraction = 100 / (100 + x_Fe_cm)

    # scale by mass_mtl/M_p --> should be smaller than x_Fe_cm
    x_Fe_c = x_Fe_cm * mantle_mass_fraction
    CMF = x_Fe_c / 100
    Fe_mass_fraction = x_Fe_c + x_Fe_mm * mantle_mass_fraction  # store bulk planet Fe fraction for posterity

    print('\ncore mass fraction = {:5.3f}\n'.format(CMF))
else:
    mantle_mass_fraction = 1 - CMF

# remove Si into core

m_Si_core = CMF * core_Si_wt * M_p
m_Si_mantle_old = wt_oxides['SiO2']/100 * mantle_mass_fraction * M_p * p.M_Si / (p.M_O * 2)
m_Si_mantle_new = m_Si_mantle_old - m_Si_core
print('m Si core', m_Si_core)
print('m Si mantle old', m_Si_mantle_old, 'new', m_Si_mantle_new)

m_Mg_old = wt_oxides['MgO']/100 * mantle_mass_fraction * M_p * p.M_Mg / p.M_O
print('m Mg mantle', m_Mg_old)

n_Mg = m_Mg_old / p.M_Mg  # convert to moles
n_Si = m_Si_mantle_new / p.M_Si
mgsi_new = n_Mg / n_Si
print('MgSi new', mgsi_new)


# renormalise
new_wt_dict = wt_oxides.copy()
new_wt_dict['SiO2'] = m_Si_mantle_new / (mantle_mass_fraction * M_p) / (p.M_Si / (p.M_O * 2)) * 100

print('-----------\nwt.% oxides\n-----------')
factor = 100 / sum(new_wt_dict.values())
for k in new_wt_dict:
    new_wt_dict[k] = new_wt_dict[k] * factor
    print("{0:<7}".format(k), "{:5.2f}%".format(new_wt_dict[k]))
print('-----------')
print('MgSi retry', get_mgsi(new_wt_dict))