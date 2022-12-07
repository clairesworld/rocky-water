import numpy as np
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(os.path.dirname(PARENT_DIR))

import perplexfugacitydata as pf
import py.main as rw

"""
source /raid1/cmg76/venv/bin/activate
cd ~/Works/rocky-water/
"""

perplex_path = '/raid1/cmg76/perple_x/'
T_iso = 1373
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
pressures_of_interest = np.linspace(p_min, p_max, 15)  # bar, for alphaMELTS
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O'] #, 'Cr2O3']
px_melt_phases = ['ctjL', 'dijL', 'enL', 'geik']

X_ferric = 0.03
# core_eff = [0.7, 0.65]  #, 0.85,  0.8, 0.7, 0.95, 0.99, 0.75, 0.9,  0.65]
core_eff = [0.88]

for ce in core_eff:
    output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(X_ferric * 100)) + 'ferric_ext/'
    output_parent_path = pf.output_parent_apollo + output_sub

    subfolders = rw.get_run_dirs(output_path=output_parent_path)
    if subfolders:
        for sub in subfolders:
            name = os.path.basename(sub)
            # star = name.split('_')[2]

            dat = pf.init_from_results(name, X_ferric, load_results_csv=True, output_parent_path=output_parent_path)
            if dat:
                dat.ferric_composition_calc(T_iso=T_iso, p_min=p_min, p_max=p_max, verbose=True)
                print('finished', name)
            else:
                print('error loading')
    else:
        print('no local output found in', output_parent_path)