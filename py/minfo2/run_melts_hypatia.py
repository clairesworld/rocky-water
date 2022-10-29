import numpy as np
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(os.path.dirname(PARENT_DIR))


import meltsfugacitydata as mf

alphamelts_path_apollo = '/raid1/cmg76/Works/alphamelts/'
output_parent_apollo = '/raid1/cmg76/alphamelts/output/rocky-fo2/'

T_iso = 1373
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
pressures_of_interest = np.linspace(p_min, p_max, 15)  # bar, for alphaMELTS
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O']

X_ferric = 0.03
core_eff = 0.88
output_sub = 'hypatia_' + str(int(core_eff * 100)) + 'coreeff_' + str(int(X_ferric * 100)) + 'ferric_ext/'
output_parent_path = output_parent_apollo + output_sub

mf.fo2_from_hypatia(pressures_of_interest, n_sample=-1, core_efficiency=core_eff, X_ferric=X_ferric,
                    T_final=T_iso, verbose=True, oxide_list=oxide_list,
                    planet_kwargs={}, compare_buffer='qfm',
                    output_parent_path=output_parent_path,
                    alphamelts_path=alphamelts_path_apollo)

