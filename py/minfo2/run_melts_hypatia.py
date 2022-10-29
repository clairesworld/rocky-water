import numpy as np
import meltsfugacitydata as mf

T_iso = 1373
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
pressures_of_interest = np.linspace(p_min, p_max, 15)  # bar, for alphaMELTS
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O']
px_melt_phases = ['ctjL', 'dijL', 'enL', 'geik']

X_ferric = 0.03
core_eff = 0.88
output_sub = 'hypatia_' + str(int(core_eff * 100)) + 'coreeff_' + str(int(X_ferric * 100)) + 'ferric_ext/'

output_parent_path = pf.output_parent_apollo + output_sub
perplex_path = '/raid1/cmg76/perple_x/'


fo2_from_hypatia(pressures_of_interest, n_sample=-1, core_efficiency=0.88, planet_kwargs={},
                     output_parent_path=output_parent_default)
mf.fo2_from_oxides(name=name, pressures_of_interest=pressures_of_interest, T_final=T_iso, test_oxides=test_oxides,
                   core_efficiency=core_eff, X_ferric=X_ferric, star=star, verbose=True, planet_kwargs={},
                   output_parent_path=output_parent_path_melts, oxide_list=oxide_list,
                   compare_buffer='qfm',)  # need oxide_list for MELTS?