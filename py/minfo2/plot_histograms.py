import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
from py.useful_and_bespoke import dark_background

""" compare fo2 hist for multiple runs """

# set these
core_eff = 88
X_ferric_list = [1, 3, 5, 7, 9]
p_of_interest = 1
model = 'perplex'
legsize = 10
labelsize = 14
ticksize = 10
labelpad = 10

if model == 'melts':
    dirs = [fo2plt.output_parent_mlt_earth + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
            X_ferric_list]
    z_scale = 100
    model_str = 'MELTS'
elif model == 'perplex':
    dirs = [fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
            X_ferric_list]
    z_scale = 1
    model_str = 'Perple_x'
if p_of_interest == 1:
    x_var = 'delta_qfm_1GPa'  # 'logfo2_1GPa'
    xlabel = 'log(fO$_2$) at 1 GPa'
    fname = 'pop_hist_Xferric_1GPa_' + model
elif p_of_interest == 4:
    x_var = 'delta_qfm_4GPa'
    xlabel = 'log(fO$_2$) at 4 GPa'
    fname = 'pop_hist_Xferric_4GPa_' + model

fig, axes = plt.subplots(len(dirs), 1, figsize=(6, 4))

for ii, ax in enumerate(axes):
    if ii == 0:
        legtitle = r'Fe$^{3+}/\Sigma$Fe'
    else:
        legtitle = ''
    fig, ax = fo2plt.compare_pop_hist([dirs[ii]], x_var=x_var, z_var='X_ferric', x_scale=1, fname=fname,
                                      hilite=0,  # figsize=(11, 3),
                                      title=None, xlabel='', legtitle=legtitle, bins=25, labelsize=labelsize,
                                      ticksize=ticksize, legsize=legsize, lw=2,
                                      exclude_silica=True, model=model, labelpad=labelpad, z_scale=z_scale,
                                      legloc='upper left',
                                      ls='-', cmap='autumn', vmin=min(X_ferric_list), vmax=max(X_ferric_list) + 3,
                                      save=False, verbose=False,
                                      fig=fig, ax=ax,
                                      # x_shift=-14.431326430547
                                      )

    ax.set_yticks([])
    if ii < len(dirs) - 1:
        ax.set_xticks([])
    ax.set_xlim(-6, 2)  # (-4, 3) for melts, (-6, 2) for perplex
    ax.set_ylim(0, 1.3)  # 2 for melts 1 GPa; 1.8 for melts 4; 1.3 for px 4

axes[-1].set_xlabel('log$f$O$_2$ ($\Delta$QFM)', fontsize=labelsize, labelpad=labelpad)
axes[0].set_title(model_str + ' at ' + str(p_of_interest) + ' GPa', fontsize=labelsize)
# ax.set_title(str(p_of_interest) + ' GPa', fontsize=18)

plt.tight_layout()
plt.subplots_adjust(hspace=0)
fig.savefig(fo2plt.figpath + 'hist_fer_' + model + str(p_of_interest) + '.png', bbox_inches='tight', facecolor=fig.get_facecolor())


""" across core partitioning """

# set these
# X_ferric = 3
# coreeff_list = [80, 85, 88, 99]
#
# dirs = [fo2plt.output_parent_px + 'hypatia_' + str(ce) + 'coreeff_' + str(X_ferric) + 'ferric_ext/' for ce in coreeff_list]
# if p_of_interest == 1:
#     x_var = 'delta_qfm_1GPa'
#     xlabel = 'log(fO$_2$) at 1 GPa'
#     fname = 'compare_pop_hist_coreeff_1GPa'
# elif p_of_interest == 4:
#     x_var = 'delta_qfm_4GPa'
#     xlabel = 'log(fO$_2$) at 4 GPa'
#     fname = 'compare_pop_hist_coreeff_4GPa'
# fig, ax = fo2plt.compare_pop_hist(dirs, x_var=x_var, z_var='core_eff', x_scale=1, z_scale=1, fname=fname,
#                         save=False, exclude_names=[], legtitle=r'Fe percentage into core',
#                         ls='-', cmap='cool', vmin=0.65, vmax=1, xlabel='log(fO$_2$) at 1 GPa', bins=40)
# ax.set_yticks([])
# ax.set_xlabel('log$f$O$_2$ ($\Delta$QFM)', fontsize=18, labelpad=12)
# ax.set_title(str(p_of_interest) + ' GPa', fontsize=18)
# ax.set_xlim(-4, 3)
# ax.set_ylim(0, 1.8)
#
# fig.savefig(fo2plt.figpath + 'hist_core.png', bbox_inches='tight', facecolor=fig.get_facecolor())

plt.show()
