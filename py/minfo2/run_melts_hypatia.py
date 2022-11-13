import numpy as np
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(os.path.dirname(PARENT_DIR))

import py.minfo2.meltsfugacitydata as mf

"""
source /raid1/cmg76/venv/bin/activate
cd ~/Works/rocky-water/
"""

# set paths
alphamelts_path_apollo = '/raid1/cmg76/alphamelts/'
opp_apollo = '/raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/'
alphamelts_path_starlite = '/home/claire/Works/alphamelts/'
opp_starlite = '/home/claire/Works/min-fo2/alphamelts_output/hypatia_local2/'
opp_galson = '/home/claire/Works/min-fo2/alphamelts_output/earth-tea23/'

# set these
T_iso = 1373
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
pressures_of_interest = np.linspace(p_min, p_max, 15)  # bar, for alphaMELTS
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O']
X_ferric = [0.01, 0.03, 0.05, 0.07, 0.09]
core_eff = [0.88, 0.85,  0.8, 0.7, 0.95, 0.99, 0.75, 0.9,  0.65]
skip_stars = ['HIP 522', 'HIP 801', 'HIP 102409']
location = 'apollo'  # 'starlite'

# run
if location == 'apollo':
    source = opp_apollo
    alphamelts_path = alphamelts_path_apollo
elif location == 'starlite':
    source = opp_starlite
    alphamelts_path = alphamelts_path_starlite
for ce in core_eff:
    for Xf in X_ferric:
        output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(Xf * 100)) + 'ferric_ext/'
        output_parent_path = source + output_sub

        # generate planet compositions from hypatia and calcualte mantle fo2
        mf.fo2_from_hypatia(pressures_of_interest, n_sample=-1, core_efficiency=ce, X_ferric=Xf,
                            T_final=T_iso, verbose=True,
                            oxide_list=oxide_list, oxides=oxide_list,  # fucked this up somewhere just give both names lol
                            planet_kwargs={},
                            # compare_buffer='qfm',
                            output_parent_path=output_parent_path,
                            alphamelts_path=alphamelts_path_starlite,
                            # perplex_path='/raid1/cmg76/perple_x/',  # for qfm data table
                            names_file='/home/claire/Works/rocky-water/py/host_names.txt',
                            # use_local_composition=False,
                            use_local_composition=True, existing_dir='hypatia_88coreeff_3ferric_ext/',  # try local first
                            existing_output_parent=opp_starlite, # '/raid1/cmg76/perple_x/output/rocky-fo2/',
                            suffix=str(Xf*100).replace('.', ',') + 'fer',
                            skip_names=skip_stars,
                            # restart='HIP 22627'
                            )
        mf.common_Tmin(output_parent_path)



