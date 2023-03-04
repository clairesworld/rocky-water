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
python3 py/minfo2/process_perplex_local.py
"""

# user input
T_of_interests = [float(x) for x in input('Enter T_of_interest, separated by spaces (default 1373): ').split()]
X_ferric = [float(x) for x in input('Enter X_ferric, separated by spaces (default 0.03): ').split()]
core_eff = [float(x) for x in input('Enter core_eff, separated by spaces (default 0.88): ').split()]
if not T_of_interests:
    T_of_interests = [1373]
if not X_ferric:
    X_ferric = [0.03]
if not core_eff:
    core_eff = [0.88]

perplex_path = '/raid1/cmg76/perple_x/'
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'Na2O']#, 'Cr2O3']
px_melt_phases = ['ctjL', 'dijL', 'enL', 'geik']

for T_of_interest in T_of_interests:
    for ce in core_eff:
        for Xf in X_ferric:
            output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(Xf * 100)) + 'ferric_ext/'
            output_parent_path = pf.output_parent_apollo + output_sub

            # recalculate fo2 or at a different temperature
            pf.fo2_from_local(output_parent_path=output_parent_path, T_iso=T_of_interest, run_werami=True,
                              do_ferric_comp=True, do_system_comp=True, p_min=p_min, p_max=p_max, X_ferric=Xf,
                              perplex_path=perplex_path, mu0_file='data_tables/mu_o2_standard.tab',
                              compare_buffer='qfm', suppress_output=True, verbose=True)
