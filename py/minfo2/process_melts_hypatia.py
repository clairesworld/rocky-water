import numpy as np
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(os.path.dirname(PARENT_DIR))

import py.minfo2.meltsfugacitydata as mfug

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
T_of_interest = 1373.15  # 1673.15
# X_ferric = [0.03, 0.02, 0.01, 0.05, 0.07, 0.09]  #[0.01, 0.03, 0.05, 0.07, 0.09]
# core_eff = [0.88, 0.85, 0.8, 0.7, 0.95, 0.99, 0.75, 0.9,  0.65]
core_eff = [0.70, 0.95, 0.99, 0.65]
X_ferric = [0.03]
location = 'apollo'  # 'starlite'

# run
if location == 'apollo':
    source = opp_apollo
    alphamelts_path = alphamelts_path_apollo
    perplex_path = '/raid1/cmg76/perple_x/'
elif location == 'starlite':
    source = opp_starlite
    alphamelts_path = alphamelts_path_starlite
    perplex_path = '/home/claire/Works/perple_x/'

for ce in core_eff:
    for Xf in X_ferric:
        output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(Xf * 100)) + 'ferric_ext/'
        output_parent_path = source + output_sub
        # output_parent_path = '/raid1/cmg76/alphamelts/output/rocky-fo2/hypatia_88coreeff_3ferric_ext_Cr/'

        mfug.common_Tmin(output_parent_path)

        # calculate mantle fo2 only
        mfug.fo2_from_local(output_parent_path, core_efficiency=ce, X_ferric=Xf, alphamelts_path=alphamelts_path,
                          compare_buffer='qfm', perplex_path=perplex_path,
                          T_of_interest=T_of_interest,  # reload_TP=True,
                          verbose=False)



""" custom """
# output_parent_path = '/raid1/cmg76/alphamelts/output/rocky-fo2/mgsi_from_earth/'
# mfug.common_Tmin(output_parent_path)
#
# # calculate mantle fo2 only
# ce, X_ferric = 0.88, 0.03
# mfug.fo2_from_local(output_parent_path, core_efficiency=ce, X_ferric=X_ferric, alphamelts_path=alphamelts_path,
#                     compare_buffer='qfm', perplex_path=perplex_path,
#                     T_of_interest=T_of_interest,  # reload_TP=True,
#                     verbose=False)




