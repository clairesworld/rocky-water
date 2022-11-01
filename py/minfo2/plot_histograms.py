import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt

""" compare fo2 hist for multiple runs """

core_eff = 88
X_ferric_list = [1, 3, 5, 7, 9]
dirs = [fo2plt.output_parent_default + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
        X_ferric_list]
fo2plt.compare_pop_hist(dirs, x_var='logfo2_1GPa', z_var='X_ferric', x_scale=1, z_scale=1,
                        fname='compare_pop_hist_Xferric',
                        legtitle=r'Fe$^{3+}$ percentage', save=True,
                        ls='-', cmap='spring', vmin=min(X_ferric_list), vmax=max(X_ferric_list), xlabel='log(fO$_2$) at 1 GPa', bins=20)

# coreeff_list = [70, 88, 99]
# dirs = [output_parent_default + 'hypatia_' + str(ce) + 'coreeff_' + str(X_ferric) + 'ferric/' for ce in coreeff_list]
# compare_pop_hist(dirs, x_var='logfo2_1GPa', z_var='core_eff', x_scale=1, z_scale=1, fname='compare_pop_hist_coreeff', figpath=figpath,
#                  save=True, exclude_names=[], legtitle=r'Fe percentage into core',
#                  ls='-', cmap='cool', vmin=0.65, vmax=1, xlabel='log(fO$_2$) at 1 GPa', bins=40)


plt.show()
