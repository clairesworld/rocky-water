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
git pull
python3 py/minfo2/process_melts_hypatia.py
"""

# set paths
alphamelts_path_apollo = '/raid1/cmg76/alphamelts/'
opp_apollo = '/raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/'
alphamelts_path_starlite = '/home/claire/Works/alphamelts/'
opp_starlite = '/home/claire/Works/min-fo2/alphamelts_output/hypatia_local2/'
opp_galson = '/home/claire/Works/min-fo2/alphamelts_output/earth-tea23/'

""" vvv UNCOMMENT TO RUN vvv """

# set these
T_of_interest = 1373.15  # 1673.15
X_ferric = [float(x) for x in input('Enter X_ferric, separated by spaces (e.g. 0.03): ').split()] #[0.07]  #, 0.03, 0.05, 0.07, 0.09]  # , 0.01, 0.05, 0.07, 0.09]  #[0.01, 0.03, 0.05, 0.07, 0.09]
core_eff = [float(x) for x in input('Enter core_eff, separated by spaces (e.g. 0.88): ').split()] #[0.88]
# core_eff = [0.88]
# X_ferric = [0.03]
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

        print('\n\n\nfinding common T min\n--------------------\n')
        mfug.common_Tmin(output_parent_path, include_p=[1e4, 4e4])

        print('\n\n\nstarting fO2 calculation\n------------------------\n')

        # calculate mantle fo2 only
        mfug.fo2_from_local(output_parent_path, core_efficiency=ce, X_ferric=Xf, alphamelts_path=alphamelts_path,
                            compare_buffer='qfm', perplex_path=perplex_path, T_of_interest=T_of_interest, save=True,
                            verbose=False,
                            # restart='2MASS19155319+4437283'
                            )

""" ^^^ UNCOMMENT TO RUN ^^^ """

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

#
""" for a min T look at all the cases that got there """
# source = opp_galson
# ce = 0.88
# Xf = 0.03
# T = 1373.15
# P = 1
# output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(Xf * 100)) + 'ferric_ext/'
# output_parent_path = source + output_sub
#
# mfug.create_isothermal_csv(output_parent_path, T, P, Xf, ce)
#
