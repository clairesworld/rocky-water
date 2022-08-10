import numpy as np
import perplexfugacity as pf

# match Stolper fig 3 @ 1373 K?
T_iso = 1373
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest

X_ferric = 0.03
core_eff = 0.88
output_sub = 'hypatia' + str(int(core_eff * 100)) + 'Fe/'
output_parent_path = pf.output_parent_apollo + output_sub
perplex_path = '/raid1/cmg76/perple_x/'

# # only need to run this once to get build files (i.e., bulk composition) and vertex output files
# note this will first download the stellar composisions for the whole sample and then run perple_x
pf.fo2_from_hypatia(p_min, p_max, n_sample=-1, T_min=T_min, T_max=T_max, isotherm=T_iso,
                    X_ferric=X_ferric, core_efficiency=core_eff, planet_kwargs={'Tp': 999},
                    check_comp=True, suppress_output=False, run=True, verbose=True,
                    output_parent_path=output_parent_path, perplex_path=perplex_path,
                    mu0_file='data_tables/mu_o2_standard.tab', compare_buffer='qfm',
                    use_local_composition=True, existing_dir=output_sub,  # try local first
                    restart='2MASS 23155829+3127462'
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
