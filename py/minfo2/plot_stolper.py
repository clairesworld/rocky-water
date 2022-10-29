import oxygen_fugacity_plots as pltfo2
import numpy as np
import perplexfugacitydata as perplex
import meltsfugacitydata as melts
import matplotlib.pyplot as plt

" try to recreate Stolper fig 3 "

T_iso = 1373
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
pressures_of_interest = np.linspace(p_min, p_max, 12)  # for alphaMELTS
X_ferric = 0.03
# core_eff = 0.88
star = None
name = 'Stolper'

wt_oxides_DMM_ext = {'SiO2': 44.71, 'MgO': 38.73, 'CaO': 3.17, 'Al2O3': 3.98, 'FeO': 8.17,
                     'Na2O': 0.28, 'TiO2': 0.13, 'Cr2O3': 0.57}  # Workman & Hart depleted mantle, extended

melts.fo2_from_oxides(name=name, pressures_of_interest=pressures_of_interest, T_final=T_iso,
                      test_oxides=wt_oxides_DMM_ext, X_ferric=X_ferric, star=star,
                      # core_efficiency=core_eff,
                      output_parent_path=melts.output_parent_default, oxide_list=list(wt_oxides_DMM_ext.keys()),
                      compare_buffer='qfm', )  # need oxide_list for MELTS?

plot_kwargs = {}
pltfo2.stolper_subplot(name=name, output_parent_path=melts.output_parent_default, show_buffer=True, save=True,
                       # p_min=1.1,  # in GPa
                       lw=2, model='melts', **plot_kwargs)

plt.show()