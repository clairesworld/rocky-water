import numpy as np
import perplexfugacitydata as pf
import meltsfugacitydata as mf

T_iso = 1373
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
pressures_of_interest = np.linspace(p_min, p_max, 15)  # bar, for alphaMELTS
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O']
px_melt_phases = ['ctjL', 'dijL', 'enL', 'geik']

X_ferric = 0.03
core_eff = 0.88
output_sub = 'hypatia_' + str(int(core_eff * 100)) + 'coreeff_' + str(int(X_ferric * 100)) + 'ferric/'

# output_parent_path = pf.output_parent_apollo + output_sub
# perplex_path = '/raid1/cmg76/perple_x/'
output_parent_path = pf.output_parent_default
perplex_path = '/home/claire/Works/perple_x/'

# # only need to run this once to get build files (i.e., bulk composition) and vertex output files
# note this will first download the stellar composisions for the whole sample and then run perple_x
pf.fo2_from_hypatia(p_min, p_max, n_sample=1, T_min=T_min, T_max=T_max, isotherm=T_iso,
                    X_ferric=X_ferric, core_efficiency=core_eff,
                    planet_kwargs={'Tp': 999, 'oxides': oxide_list, 'excluded_phases': px_melt_phases},
                    #solve_interior=False, --> already a parameter
                    check_comp=True, suppress_output=False, run=True, verbose=True,
                    output_parent_path=output_parent_path, perplex_path=perplex_path,
                    mu0_file='data_tables/mu_o2_standard.tab', compare_buffer='qfm',
                    use_local_compositon=False, names_file='/home/claire/Works/rocky-water/py/host_names.txt',
                    # use_local_composition=True, existing_dir='hypatia88Fe/',  # try local first
                    # existing_output_parent=pf.output_parent_apollo,  <== existing kwarg
                    # restart='2MASS 23155829+3127462'
                    )

# # can (re)do fo2 calculations with this using existing vertex files
# pf.fo2_from_local(output_parent_path=output_parent_path, mu0_file='data_tables/mu_o2_standard.tab',
#                   compare_buffer='qfm', check_comp=True, suppress_output=False, isotherm=T_iso)
#

# do solar and Earth
# for (name, test_oxides, star) in zip(['dmm', 'sun_88Fe'], [wt_oxides_DMM, None], [None, 'sun']):
#     pf.fo2_from_oxides(name=name, p_min=p_min, p_max=p_max, T_min=T_min, T_max=T_max, test_oxides=test_oxides,
#                     X_ferric=X_ferric, isotherm=T_iso, core_efficiency=core_eff, star=star,
#                     run=True, compare_buffer='qfm',
#                     suppress_output=False, check_comp=True, verbose=True,
#                     mu0_file='data_tables/mu_o2_standard.tab', output_parent_path=pf.output_parent_default,
#                     )


# # # like Stolper
# pf.fo2_from_oxides(name='Stolper', p_min=p_min, p_max=p_max, T_min=T_min, T_max=T_max, test_oxides=pf.wt_oxides_DMM_ext,
#                 X_ferric=0.031, isotherm=T_iso,
#                 run=True, compare_buffer='qfm',
#                 suppress_output=False, check_comp=True, verbose=True,
#                 mu0_file='data_tables/mu_o2_standard.tab', output_parent_path=pf.output_parent_default,
#                 )


# # do pMELTS calc
# oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O']
# output_parent_path_melts = '/home/claire/Works/min-fo2/alphamelts_output/' + output_sub
# mf.fo2_from_hypatia(pressures_of_interest, n_sample=1, T_final=T_iso, oxides=oxide_list,
#                     X_ferric=X_ferric, core_efficiency=core_eff, planet_kwargs={'Tp': 999},
#                     #solve_interior=False, --> already a parameter
#                     suppress_output=False, verbose=True,
#                     output_parent_path=output_parent_path_melts,
#                     compare_buffer='qfm',
#                     use_local_compositon=False,
#                     # use_local_composition=True, existing_dir='hypatia88Fe/',  # try local first
#                     # existing_output_parent=pf.output_parent_apollo,  <== existing kwarg
#                     # restart='2MASS 23155829+3127462'
#                     )


"""direct compare alphamelts and perplex"""

# name = 'dmm_ext'
# star = None
# test_oxides = pf.wt_oxides_DMM_ext
# mf.fo2_from_oxides(name=name, pressures_of_interest=pressures_of_interest, T_final=T_iso, test_oxides=test_oxides,
#                    core_efficiency=core_eff, X_ferric=X_ferric, star=star, verbose=True, planet_kwargs={},
#                    output_parent_path=output_parent_path_melts, oxide_list=oxide_list,
#                    compare_buffer='qfm',)  # need oxide_list for MELTS?

# pf.fo2_from_oxides(name=name, p_min=p_min, p_max=p_max, T_min=T_min, T_max=T_max, test_oxides=test_oxides,
#                    X_ferric=X_ferric, isotherm=T_iso, core_efficiency=core_eff, star=star,
#                    run=True, compare_buffer='qfm',
#                    suppress_output=True, check_comp=True, verbose=True,
#                    mu0_file='data_tables/mu_o2_standard.tab', output_parent_path=pf.output_parent_default,
#                    )
