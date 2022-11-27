import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
from py.useful_and_bespoke import dark_background

""" compare fo2 hist for multiple runs """

# set these
core_eff = 88
X_ferric_list = [1, 3, 5, 7, 9]
p_of_interest = 1
model = 'melts'

if model == 'melts':
    dirs = [fo2plt.output_parent_mlt_earth + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
            X_ferric_list]
    z_scale=100
elif model == 'perplex':
    dirs = [fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
            X_ferric_list]
    z_scale=1
if p_of_interest == 1:
    x_var = 'delta_qfm_1GPa'  # 'logfo2_1GPa'
    xlabel = 'log(fO$_2$) at 1 GPa'
    fname = 'pop_hist_Xferric_1GPa_' + model
elif p_of_interest == 4:
    x_var = 'delta_qfm_4GPa'
    xlabel = 'log(fO$_2$) at 4 GPa'
    fname = 'pop_hist_Xferric_4GPa_' + model

for ii, Xf in enumerate(X_ferric_list):
    fig, ax = fo2plt.compare_pop_hist(dirs, x_var=x_var, z_var='X_ferric', x_scale=1, fname=fname,  xlabel=xlabel,
                                      hilite=ii, figsize=(11, 3),
                            legtitle=r'Fe$^{3+}/\Sigma$Fe', bins=25, labelsize=18, ticksize=14,
                            exclude_silica=True, model=model, labelpad=12,z_scale=z_scale, legloc='upper left',
                            ls='-', cmap='autumn', vmin=min(X_ferric_list), vmax=max(X_ferric_list), save=False, verbose=False,
                                  # x_shift=-14.431326430547
                                      )

    ax.set_yticks([])
    ax.set_xlabel('log$f$O$_2$ ($\Delta$QFM)', fontsize=18, labelpad=12)
    ax.set_title('1 GPa', fontsize=18)
    ax.set_xlim(-4, 3)
    ax.set_ylim(0, 1.8)

    fig, ax = dark_background(fig, ax)
    fig.savefig(fo2plt.figpath + 'hist_fer_dark' + str(ii) + '.png', bbox_inches='tight', facecolor=fig.get_facecolor())


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