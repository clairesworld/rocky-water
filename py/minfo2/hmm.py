import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
import numpy as np
import meltsfugacitydata as mfug
import perplexfugacitydata as pfug

""" check on certain cases """

core_eff = 88
Xf = 3.0
star = '2MASS19153994+3935409'

# name = '1M_' + str(core_eff) + 'Ceff_' + star + '_999K_' + str(Xf).replace('.', ',') + 'fer'
name = 'Stolper'

opp_mlt = mfug.output_parent_default  # fo2plt.output_parent_mlt + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'
opp_px = pfug.output_parent_default  #fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'

""" look @ melts stuff """

print('\n\n-----------------------------------------------------------\nmelts data')
mdat = mfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp_mlt, verbose=False)
mdat.print_comp()



""" look @ px stuff """

print('\n\n-----------------------------------------------------------\nplerplex data')
pdat = pfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp_px, verbose=False)
pdat.print_comp()
