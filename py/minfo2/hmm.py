import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
import numpy as np
import meltsfugacitydata as mfug
import perplexfugacitydata as pfug
import py.bulk_composition as bulk

""" check on certain cases """

core_eff = 88
Xf = 3.0
star = '2MASS19153994+3935409'
# name = '1M_' + str(core_eff) + 'Ceff_' + star + '_999K_' + str(Xf).replace('.', ',') + 'fer'
name = 'dmm_ext'  #'Stolper'

opp_mlt = mfug.output_parent_default  # fo2plt.output_parent_mlt + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'
opp_px = pfug.output_parent_default  #fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'

""" look @ melts stuff """

# print('\n\n-----------------------------------------------------------\nmelts data')
# mdat = mfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp_mlt, verbose=False)
# mdat.print_comp()


""" look @ px stuff """

print('\n\n-----------------------------------------------------------\nplerplex data')
pdat = pfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp_px, load_results_csv=True, verbose=False)
pdat.print_comp()
# pdat.ferric_composition_calc(T_iso=1373,  verbose=True)
# df = pdat.read_ferric_phase_comp(phase='Opx', T_of_interest=1373, p_of_interest=1)
# print(df.head())

d = pdat.read_phase_comp(p_of_interest=1, T_of_interest=1373, component='Fe2O3', phases=pfug.solution_phases_default,
                        absolute_abundance=True, verbose=False)
print(d)

bulk.get_element_ratio('Mg/Si', pdat.wt_oxides)