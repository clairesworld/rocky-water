import numpy as np
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(os.path.dirname(PARENT_DIR))

import perplexfugacitydata as pf

"""
source /raid1/cmg76/venv/bin/activate
cd ~/Works/rocky-water/
"""

perplex_path = '/raid1/cmg76/perple_x/'
T_iso = 1373
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
pressures_of_interest = np.linspace(p_min, p_max, 15)  # bar, for alphaMELTS
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O']
px_melt_phases = ['ctjL', 'dijL', 'enL', 'geik']

X_ferric = 0.07
core_eff = [0.88, 0.85,  0.8, 0.7, 0.95, 0.99, 0.75, 0.9,  0.65]

for ce in core_eff:
    output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(X_ferric * 100)) + 'ferric_ext/'
    output_parent_path = pf.output_parent_apollo + output_sub

    # only need to run this once to get build files (i.e., bulk composition) and vertex output files
    # note this will first download the stellar composisions for the whole sample and then run perple_x
    pf.fo2_from_hypatia(p_min, p_max, n_sample=-1, T_min=T_min, T_max=T_max, isotherm=T_iso,
                        X_ferric=X_ferric, core_efficiency=ce,
                        planet_kwargs={'Tp': 999, 'oxides': oxide_list, 'excluded_phases': px_melt_phases},
                        #solve_interior=False, --> already a parameter
                        check_comp=True, suppress_output=False, run=True, verbose=True,
                        output_parent_path=output_parent_path, perplex_path=perplex_path,
                        mu0_file='data_tables/mu_o2_standard.tab', compare_buffer='qfm',
                        names_file='/home/cmg76/Works/rocky-water/py/host_names.txt',
                        # use_local_compositon=False,
                        use_local_composition=True, existing_dir='hypatia_88coreeff_3ferric_ext/',  # try local first
                        # existing_output_parent=pf.output_parent_apollo,  # <== existing kwarg
                        # restart='2MASS 23155829+3127462'
                        run_vertex=True,
                        )
