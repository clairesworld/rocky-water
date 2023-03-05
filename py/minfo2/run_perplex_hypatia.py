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
python3 py/minfo2/run_perplex_hypatia.py
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
# TiO2 needs to be in oxides to find Na, it won't get inxluded in final px calculation
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O']  #, 'Cr2O3']
px_melt_phases = ['ctjL', 'dijL', 'enL', 'geik']

# star compositions which had weird phases in melts, going to ignore here too
skip_stars = ['2MASS 19141179+3833548']

for T_of_interest in T_of_interests:
    for ce in core_eff:
        for Xf in X_ferric:
            output_sub = 'hypatia_' + str(int(ce * 100)) + 'coreeff_' + str(int(Xf * 100)) + 'ferric_ext/'
            output_parent_path = pf.output_parent_apollo + output_sub

            # only need to run this once to get build files (i.e., bulk composition) and vertex output files
            # note this will first download the stellar composisions for the whole sample and then run perple_x
            pf.fo2_from_hypatia(p_min, p_max, n_sample=-1, T_min=T_min, T_max=T_max, T_iso=T_of_interest,
                                X_ferric=Xf, core_efficiency=ce,
                                planet_kwargs={'Tp': 999, 'oxides': oxide_list, 'excluded_phases': px_melt_phases},
                                #solve_interior=False, --> already a parameter
                                do_system_comp=True,  # uses werami_command_comp after fo2 calc
                                suppress_output=False, run=True, verbose=True,
                                output_parent_path=output_parent_path, perplex_path=perplex_path,
                                mu0_file='data_tables/mu_o2_standard.tab', compare_buffer='qfm',
                                names_file='/home/cmg76/Works/rocky-water/py/host_names.txt',
                                # use_local_compositon=False,
                                use_local_composition=True, existing_dir='earth-tea23/hypatia_95coreeff_3ferric_ext/',
                                # try local first, 1 ferric
                                existing_output_parent='/raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/',
                                # restart='2MASS 19243554+4040098',
                                run_vertex='auto',  # overwrite existing vertex files if exist
                                skip_existing=True,  # do not do anything if directory exists with *_results.csv
                                suffix=str(Xf*100).replace('.', ',') + 'fer',
                                skip_stars=skip_stars,
                                # names=['1M_95Ceff_2MASS19421779+4248231_999K_3,0fer']
                                )

            # pf.fo2_from_hypatia_1D(p_min, p_max, n_sample=-1, T_min=T_min, T_max=T_max, T_iso=T_iso,
            #                     X_ferric=X_ferric, core_efficiency=ce,
            #                     planet_kwargs={'Tp': 999, 'oxides': oxide_list, 'excluded_phases': px_melt_phases},
            #
            #                     #solve_interior=False, --> already a parameter
            #                     check_comp=True,  # uses werami_command_comp after fo2 calc
            #                     suppress_output=False, run=True, verbose=True,
            #                     output_parent_path=output_parent_path, perplex_path=perplex_path,
            #                     mu0_file='data_tables/mu_o2_standard.tab', compare_buffer='qfm',
            #                     names_file='/home/cmg76/Works/rocky-water/py/host_names.txt',
            #                     # use_local_compositon=False,
            #                     use_local_composition=True, existing_dir='hypatia_88coreeff_3ferric_ext_Cr/',  # try local first
            #                     existing_output_parent='/raid1/cmg76/alphamelts/output/rocky-fo2/',
            #                     # existing_output_parent='/raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/',
            #                     # existing_output_parent=pfug.output_parent_apollo,  # <== existing kwarg
            #                     restart='HIP 3479',
            #                     run_vertex='auto',  # overwrite existing vertex files if exist
            #                     skip_existing=True,  # do not do anything if directory exists with *_results.csv
            #                     suffix=str(X_ferric*100).replace('.', ',') + 'fer',
            #                     )
