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
    xlim = (-4, 3)
elif model == 'perplex':
    dirs = [fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
            X_ferric_list]
    z_scale = 1
    model_str = 'Perple_X'
    xilm = (-6, 2)
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
                                      legloc='lower left',
                                      ls='-', cmap='Greys', alpha=0.7,  #'autumn',
                                      vmin=min(X_ferric_list) - 9, vmax=max(X_ferric_list),
                                      # vmin=min(X_ferric_list), vmax=max(X_ferric_list) + 3,
                                      save=False, verbose=False,
                                      fig=fig, ax=ax,
                                      x_shift=-3.5,
                                      # x_shift=-14.431326430547
                                      )

    ax.set_yticks([])
    if ii < len(dirs) - 1:
        ax.set_xticks([])
    ax.set_xlim(-3, 4)  # (-4, 3) for melts, (-6, 2) for perplex
    ax.set_ylim(0, 1.3)  # 2 for melts 1 GPa; 1.8 for melts 4; 1.3 for px 4

axes[-1].set_xlabel('log$f$O$_2$ ($\Delta$FMQ) at ' + str(p_of_interest) + ' GPa', fontsize=labelsize, labelpad=labelpad)
# axes[0].set_title(model_str + ' at ' + str(p_of_interest) + ' GPa', fontsize=labelsize)
# ax.set_title(str(p_of_interest) + ' GPa', fontsize=18)

plt.tight_layout()
plt.subplots_adjust(hspace=0)
fig.savefig(fo2plt.figpath + 'hist_fer_' + model + str(p_of_interest) + '.png', bbox_inches='tight', facecolor=fig.get_facecolor(),
            dpi=400)





""" across core partitioning """
#
# # set these
# Xf = 3
# p_of_interest = 1
# model = 'perplex'
# legsize = 10
# labelsize = 14
# ticksize = 10
# labelpad = 10
#
# if model == 'melts':
#     core_eff = [99, 95, 88, 85, 80, 70]
#     dirs = [fo2plt.output_parent_mlt_earth + 'hypatia_' + str(ce) + 'coreeff_' + str(Xf) + 'ferric_ext/' for ce in
#             core_eff]
#     model_str = 'MELTS'
#     xlim = -4, 3
# elif model == 'perplex':
#     core_eff = [99, 88, 85, 80, 70]
#     dirs = [fo2plt.output_parent_px + 'hypatia_' + str(ce) + 'coreeff_' + str(Xf) + 'ferric_ext/' for ce in
#             core_eff]
#     model_str = 'Perple_X'
#     xlim = -6, 2
# if p_of_interest == 1:
#     x_var = 'delta_qfm_1GPa'  # 'logfo2_1GPa'
#     xlabel = 'log(fO$_2$) at 1 GPa'
#     # fname = 'pop_hist_Xferric_1GPa_' + model
# elif p_of_interest == 4:
#     x_var = 'delta_qfm_4GPa'
#     xlabel = 'log(fO$_2$) at 4 GPa'
#     # fname = 'pop_hist_Xferric_4GPa_' + model
#
# fig, axes = plt.subplots(len(dirs), 1, figsize=(6, 4))
#
# for ii, ax in enumerate(axes):
#     if ii == 0:
#         legtitle = '' #r'$\chi^{\rm Fe}_{\rm core}$'
#     else:
#         legtitle = ''
#     fig, ax = fo2plt.compare_pop_hist([dirs[ii]], x_var=x_var, z_var='core_eff', z_scale=100, x_scale=1, #fname=fname,
#                                       hilite=0,  # figsize=(11, 3),
#                                       title=None, xlabel='', legtitle=legtitle, bins=25, labelsize=labelsize,
#                                       ticksize=ticksize, legsize=legsize, lw=2,
#                                       exclude_silica=True, model=model, labelpad=labelpad,
#                                       legloc='upper left', #alpha=0.9,
#                                       ls='-', cmap='autumn',
#                                       # vmin=min(core_eff) - 15, vmax=max(core_eff) - 1,
#                                       vmin=min(core_eff), vmax=max(core_eff) + 3,
#                                       save=False, verbose=False,
#                                       fig=fig, ax=ax,
#                                       # x_shift=-3.5,
#                                       # x_shift=-14.431326430547
#                                       )
#
#     ax.set_yticks([])
#     if ii < len(dirs) - 1:
#         ax.set_xticks([])
#     ax.set_xlim(xlim)  # (-4, 3) for melts, (-6, 2) for perplex
#     ax.set_ylim(0, 1.3)  # 2 for melts 1 GPa; 1.8 for melts 4; 1.3 for px 4
#
# axes[-1].set_xlabel('log$f$O$_2$ ($\Delta$FMQ) at ' + str(p_of_interest) + ' GPa', fontsize=labelsize, labelpad=labelpad)
# # axes[0].set_title(model_str + ' at ' + str(p_of_interest) + ' GPa', fontsize=labelsize)
# # ax.set_title(str(p_of_interest) + ' GPa', fontsize=18)
#
# plt.tight_layout()
# plt.subplots_adjust(hspace=0)
# fig.savefig(fo2plt.figpath + 'hist_core_' + model + str(p_of_interest) + '.png', bbox_inches='tight', facecolor=fig.get_facecolor(),
#             dpi=300)

plt.show()
