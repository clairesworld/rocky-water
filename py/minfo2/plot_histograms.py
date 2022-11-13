import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt

""" compare fo2 hist for multiple runs """

# set these
core_eff = 88
X_ferric_list = [1, 3, 5, 7, 9]
p_of_interest = 4

dirs = [fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
        X_ferric_list]
if p_of_interest == 1:
    x_var = 'logfo2_1GPa'
    xlabel = 'log(fO$_2$) at 1 GPa'
    fname = 'compare_pop_hist_Xferric_1GPa'
elif p_of_interest == 4:
    x_var = 'logfo2_4GPa'
    xlabel = 'log(fO$_2$) at 4 GPa'
    fname = 'compare_pop_hist_Xferric_4GPa'
fo2plt.compare_pop_hist(dirs, x_var=x_var, z_var='X_ferric', x_scale=1, z_scale=1, fname=fname,  xlabel=xlabel,
                        legtitle=r'Fe$^{3+}/\Sigma$Fe (%)', bins=10,
                        exclude_silica=True,
                        ls='-', cmap='spring', vmin=min(X_ferric_list), vmax=max(X_ferric_list), save=True)

""" across core partitioning """

# # set these
# X_ferric = 3
# coreeff_list = [70, 85, 88, 99]
# p_of_interest = 4
#
# dirs = [fo2plt.output_parent_px + 'hypatia_' + str(ce) + 'coreeff_' + str(X_ferric) + 'ferric_ext/' for ce in coreeff_list]
# if p_of_interest == 1:
#     x_var = 'logfo2_1GPa'
#     xlabel = 'log(fO$_2$) at 1 GPa'
#     fname = 'compare_pop_hist_coreeff_1GPa'
# elif p_of_interest == 4:
#     x_var = 'logfo2_4GPa'
#     xlabel = 'log(fO$_2$) at 4 GPa'
#     fname = 'compare_pop_hist_coreeff_4GPa'
# fo2plt.compare_pop_hist(dirs, x_var=x_var, z_var='core_eff', x_scale=1, z_scale=1, fname=fname,
#                         save=True, exclude_names=[], legtitle=r'Fe percentage into core',
#                         ls='-', cmap='cool', vmin=0.65, vmax=1, xlabel='log(fO$_2$) at 1 GPa', bins=40)

plt.show()
