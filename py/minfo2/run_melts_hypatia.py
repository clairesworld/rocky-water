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
python3 py/minfo2/run_melts_hypatia.py

cd /raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/hypatia_95coreeff_3ferric_ext/
sh runstars.sh

"""

# set paths
alphamelts_path_apollo = '/raid1/cmg76/alphamelts/'
opp_apollo = '/raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/'
alphamelts_path_starlite = '/home/claire/Works/alphamelts/'
opp_starlite = '/home/claire/Works/min-fo2/alphamelts_output/hypatia_local2/'
opp_galson = '/home/claire/Works/min-fo2/alphamelts_output/earth-tea23/'

# set these - unlikely to change
T_final = 1373.15
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
pressures_of_interest = np.linspace(p_min, p_max, 15)  # bar, for alphaMELTS
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O']  #, 'Cr2O3']
skip_stars = ['2MASS 19141179+3833548']  #['HIP 522', 'HIP 801', 'HIP 102409']
location = 'apollo'

# user input
X_ferric = [float(x) for x in input('Enter X_ferric, separated by spaces (default 0.03): ').split()] #[0.07]  #, 0.03, 0.05, 0.07, 0.09]  # , 0.01, 0.05, 0.07, 0.09]  #[0.01, 0.03, 0.05, 0.07, 0.09]
core_eff = [float(x) for x in input('Enter core_eff, separated by spaces (default 0.88): ').split()] #[0.88]
if not X_ferric:
    X_ferric = [0.03]
if not core_eff:
    core_eff = [0.88]

# run
if location == 'apollo':
    source = opp_apollo
    alphamelts_path = alphamelts_path_apollo
    perplex_path = '/raid1/cmg76/perple_x/'
    names_file = '/home/cmg76/Works/rocky-water/py/host_names.txt'
elif location == 'starlite':
    source = opp_starlite
    alphamelts_path = alphamelts_path_starlite
    perplex_path = '/home/claire/Works/perple_x/'
    names_file = '/home/claire/Works/rocky-water/py/host_names.txt'
for ce in core_eff:
    for Xf in X_ferric:
        output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(Xf * 100)) + 'ferric_ext/'
        output_parent_path = source + output_sub

        # generate planet compositions from hypatia and calcualte mantle fo2
        mf.fo2_from_hypatia(pressures_of_interest, n_sample=-1, core_efficiency=ce, X_ferric=Xf,
                            T_final=T_final, verbose=True,
                            oxide_list=oxide_list, oxides=oxide_list,  # fucked this up somewhere just give both names lol
                            planet_kwargs={},
                            compare_buffer='qfm',
                            output_parent_path=output_parent_path,
                            alphamelts_path=alphamelts_path,
                            perplex_path=perplex_path,  # for qfm data table
                            names_file=names_file,
                            # use_local_composition=False,
                            use_local_composition=True, existing_dir='hypatia_88coreeff_1ferric_ext/',  # try local first
                            # use_local_composition=True, existing_dir=output_sub,
                            existing_output_parent=source,  # '/raid1/cmg76/perple_x/output/rocky-fo2/',
                            suffix=str(Xf*100).replace('.', ',') + 'fer',
                            skip_names=skip_stars,  # []
                            # restart='2MASS 19155319+4437283',
                            dry_setup=True,  # dry_setup always True for batch melts calculations
                            )
        restart = None  # reset






""" other """
# wt_oxides_DMM_ext = {'SiO2': 44.71, 'MgO': 38.73, 'CaO': 3.17, 'Al2O3': 3.98, 'FeO': 8.17,
#                      'Na2O': 0.28, 'TiO2': 0.13}  # Workman & Hart depleted mantle, extended
# from py.main import update_MgSi
# ce, Xf = 0.88, 0.03
# output_parent_path = '/raid1/cmg76/alphamelts/output/rocky-fo2/mgsi_from_earth/'
# for mgsi in np.linspace(0.75, 1.5, num=25):
#     test_oxides = update_MgSi(mgsi, wt_oxides_DMM_ext)
#     name = 'mgsi' + str(mgsi).replace('.', ',')
#     # generate planet compositions from mgsi
#     mf.fo2_from_oxides(name, pressures_of_interest, core_efficiency=ce, X_ferric=Xf,
#                         T_final=T_iso, verbose=True,
#                         oxide_list=oxide_list, oxides=oxide_list,  # fucked this up somewhere just give both names lol
#                         test_oxides=test_oxides,
#                         planet_kwargs={},
#                         compare_buffer='qfm',
#                         output_parent_path=output_parent_path,
#                         alphamelts_path=alphamelts_path,
#                         perplex_path=perplex_path,  # for qfm data table
#                         suffix=str(Xf * 100).replace('.', ',') + 'fer',
#                         skip_names=skip_stars,
#                         # restart='HIP 22627'
#                         dry_setup=True,
#                         )
#     # mf.common_Tmin(output_parent_path)
#
# #
# import numpy as np
# import sys
# import os
#
# SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# PARENT_DIR = os.path.dirname(SCRIPT_DIR)
# sys.path.append(os.path.dirname(PARENT_DIR))
#
# import py.minfo2.meltsfugacitydata as mf
#
# """
# source /raid1/cmg76/venv/bin/activate
# cd ~/Works/rocky-water/
# """
#
# # set paths
# alphamelts_path_apollo = '/raid1/cmg76/alphamelts/'
# opp_apollo = '/raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/'
# alphamelts_path_starlite = '/home/claire/Works/alphamelts/'
# opp_starlite = '/home/claire/Works/min-fo2/alphamelts_output/hypatia_local2/'
# opp_galson = '/home/claire/Works/min-fo2/alphamelts_output/earth-tea23/'
#
# # set these
# T_iso = 1373
# p_min, p_max = 1e4, 4e4
# T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
# pressures_of_interest = np.linspace(p_min, p_max, 15)  # bar, for alphaMELTS
# oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O']  #, 'Cr2O3']
# X_ferric = [0.03]  # [0.01, 0.03, 0.05, 0.07, 0.09]
# # core_eff = [0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99]
# core_eff = [0.99]
# skip_stars = []  #['HIP 522', 'HIP 801', 'HIP 102409']
# location = 'apollo'  # 'starlite'
#
# # run
# if location == 'apollo':
#     source = opp_apollo
#     alphamelts_path = alphamelts_path_apollo
#     perplex_path = '/raid1/cmg76/perple_x/'
#     names_file = '/home/cmg76/Works/rocky-water/py/host_names.txt'
# elif location == 'starlite':
#     source = opp_starlite
#     alphamelts_path = alphamelts_path_starlite
#     perplex_path = '/home/claire/Works/perple_x/'
#     names_file = '/home/claire/Works/rocky-water/py/host_names.txt'
# for ce in core_eff:
#     for Xf in X_ferric:
#         output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(Xf * 100)) + 'ferric_ext/'
#         output_parent_path = source + output_sub
#
#         # generate planet compositions from hypatia and calcualte mantle fo2
#         mf.fo2_from_hypatia(pressures_of_interest, n_sample=-1, core_efficiency=ce, X_ferric=Xf,
#                             T_final=T_iso, verbose=True,
#                             oxide_list=oxide_list, oxides=oxide_list,  # fucked this up somewhere just give both names lol
#                             planet_kwargs={},
#                             compare_buffer='qfm',
#                             output_parent_path=output_parent_path,
#                             alphamelts_path=alphamelts_path,
#                             perplex_path=perplex_path,  # for qfm data table
#                             names_file=names_file,
#                             # use_local_composition=False,
#                             use_local_composition=True, existing_dir='hypatia_88coreeff_9ferric_ext/',  # try local first
#                             existing_output_parent=source,  # '/raid1/cmg76/perple_x/output/rocky-fo2/',
#                             suffix=str(Xf*100).replace('.', ',') + 'fer',
#                             skip_names=skip_stars,
#                             # restart='HIP 22627'
#                             dry_setup=True,
#                             )
#         # mf.common_Tmin(output_parent_path)
# #
# # wt_oxides_DMM_ext = {'SiO2': 44.71, 'MgO': 38.73, 'CaO': 3.17, 'Al2O3': 3.98, 'FeO': 8.17,
# #                      'Na2O': 0.28, 'TiO2': 0.13}  # Workman & Hart depleted mantle, extended
# # from py.main import update_MgSi
# # ce, Xf = 0.88, 0.03
# # output_parent_path = '/raid1/cmg76/alphamelts/output/rocky-fo2/mgsi_from_earth/'
# # for mgsi in np.linspace(0.75, 1.5, num=25):
# #     test_oxides = update_MgSi(mgsi, wt_oxides_DMM_ext)
# #     name = 'mgsi' + str(mgsi).replace('.', ',')
# #     # generate planet compositions from mgsi
# #     mf.fo2_from_oxides(name, pressures_of_interest, core_efficiency=ce, X_ferric=Xf,
# #                         T_final=T_iso, verbose=True,
# #                         oxide_list=oxide_list, oxides=oxide_list,  # fucked this up somewhere just give both names lol
# #                         test_oxides=test_oxides,
# #                         planet_kwargs={},
# #                         compare_buffer='qfm',
# #                         output_parent_path=output_parent_path,
# #                         alphamelts_path=alphamelts_path,
# #                         perplex_path=perplex_path,  # for qfm data table
# #                         suffix=str(Xf * 100).replace('.', ',') + 'fer',
# #                         skip_names=skip_stars,
# #                         # restart='HIP 22627'
# #                         dry_setup=True,
# #                         )
# #     # mf.common_Tmin(output_parent_path)
# #
# #
