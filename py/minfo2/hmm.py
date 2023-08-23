import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
import numpy as np
import meltsfugacitydata as mfug
import perplexfugacitydata as pfug
import py.bulk_composition as bulk

wt_oxides_DMM_ext = {'SiO2': 44.71, 'MgO': 38.73, 'CaO': 3.17, 'Al2O3': 3.98, 'FeO': 8.17,
                     'Na2O': 0.28, 'TiO2': 0.13, 'Cr2O3': 0.57}  # Workman & Hart depleted mantle, extended

opp_mlt = mfug.output_parent_default  # fo2plt.output_parent_mlt + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'
opp_px = pfug.output_parent_default  # fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'

""" check on certain cases """

T_iso = 1373
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1374.5  # endpoint can't equal T_of_interest
pressures_of_interest = np.linspace(p_min, p_max, 2)  # for alphaMELTS
# core_eff = 0.88
Xf = 0.03
star = None
# name = '1M_' + str(core_eff) + 'Ceff_' + star + '_999K_' + str(Xf).replace('.', ',') + 'fer'
wt_ox = wt_oxides_DMM_ext

for (name, core_eff) in zip(['Earth_ce7'], [0.70]):
    wt_ox = {'SiO2': 44.71, 'MgO': 38.73, 'CaO': 3.17, 'Al2O3': 3.98, 'FeO': 1.0,
     'Na2O': 0.28, 'TiO2': 0.13, 'Cr2O3': 0.57}  # Earth but took out FeO
    mfug.fo2_from_oxides(name=name, pressures_of_interest=pressures_of_interest, T_final=T_iso,
                         test_oxides=wt_ox, X_ferric=Xf, star=star, core_efficiency=core_eff,
                         oxide_list=list(wt_ox.keys()), oxides=list(wt_ox.keys()),
                         compare_buffer='qfm', output_parent_path=opp_mlt, )
    pfug.fo2_from_oxides(name=name, T_iso=T_iso, p_min=p_min, p_max=p_max,
                         test_oxides=wt_ox, X_ferric=Xf, star=star, core_efficiency=core_eff,
                         oxide_list=list(wt_ox.keys()), oxides=list(wt_ox.keys()),
                         compare_buffer='qfm', output_parent_path=opp_px,
                         run_on_grid=False,  # not fully tested if this works
                         points_file='/home/claire/Works/perple_x/points_files/1373K_isotherm_endpoints.xy')


""" look @ melts stuff """

# print('\n\n-----------------------------------------------------------\nmelts data')
# mdat = mfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp_mlt, verbose=False)
# mdat.print_comp()


""" look @ px stuff """

# print('\n\n-----------------------------------------------------------\nplerplex data')
# pdat = pfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp_px, load_results_csv=True, verbose=False)
# pdat.print_comp()
# # pdat.ferric_composition_calc(T_iso=1373,  verbose=True)
# # df = pdat.read_ferric_phase_comp(phase='Opx', T_of_interest=1373, p_of_interest=1)
# # print(df.head())
#
# d = pdat.read_phase_main_components(p_of_interest=1, T_of_interest=1373, component='Fe2O3',
#                                     phases=pfug.solution_phases_default, absolute_abundance=True, verbose=False)
# print(d)
#
# bulk.get_element_ratio('Mg/Si', pdat.wt_oxides)
