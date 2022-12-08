import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
from py.useful_and_bespoke import cornertext

""" make cross plot of oxides / elements"""

labelsize = 22
core_eff = 88
Xf = 3
exclude_silica = True
model = 'melts'
ncols = 4

""" plot fo2 vs. element ratios (cols) for two pressures (rows) """
""" oolour by Fe3+ content in opx?"""

fig = plt.figure(figsize=(ncols * 4, 8))
gs = fig.add_gridspec(2, ncols + 1, width_ratios=[10] * ncols + [1], height_ratios=[10, 10],
                      left=0.1, right=0.9,
                      wspace=0.05, hspace=0.05)

for jj, p_of_interest in enumerate((1, 4)):
    if model == 'melts':
        source = fo2plt.output_parent_mlt_earth
        ylims = [(-11, -8), (-8, -5)]  # 1 GPa, 4 GPa
        c = 'xkcd:midnight blue'
        title = 'pMELTS'
    elif model == 'perplex':
        source = fo2plt.output_parent_px
        ylims = [(-13.5, -8.7), (-10, -6.5)]
        c = 'xkcd:brick red'
        title = 'Perple_X'
    output_parent_path = source + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/'

    axes = []  # create new ax list just for this row
    [axes.append(fig.add_subplot(gs[jj, ii])) for ii in range(ncols)]

    fig, axs = fo2plt.element_xplot(p_of_interest=p_of_interest, components=['Mg/Si', 'Fe/Si', 'Al/Si', 'Ca/Si'],
                                    output_parent_path=output_parent_path,
                                    xlim=[(0.8, 1.65), (0.05, 0.2), (0, 0.18), (0.025, 0.13)],
                                    ylim=ylims[jj], c=c, s=30, alpha=0.4, labelsize=labelsize, save=False, fig=fig, axes=axes,
                                    model=model, verbose=False,
                                    exclude_silica=exclude_silica)

    axs[0] = cornertext(axs[0], str(p_of_interest) + ' GPa', size=labelsize, pos='top left')
    axs[1].set_xticks([0.075, 0.125, 0.175])  # Fe/Si
    if jj == 0:
        for ax in axs:
            ax.set_xticks([])
            ax.set_xlabel(None)

fig.suptitle(title, fontsize=labelsize, y=0.92)
fig.savefig(fo2plt.figpath + 'crossplot_elements_' + model + '.png', bbox_inches='tight')
# plt.show()

""" other phase abundances """

# if p_of_interest == 1:
#     ylim = (-14, -9)
#     phcomps = ['Ol', 'Opx', 'Cpx', 'Sp']
# elif p_of_interest == 4:
#     ylim = (-11.5, -7)
#     phcomps = ['Ol', 'Opx', 'Cpx', 'Gt']
# if not exclude_silica:
#     phcomps.append('q')
# fo2plt.element_xplot(p_of_interest=p_of_interest, components=phcomps,
#                      output_parent_path=output_parent_path,
#                      ylim=ylim, linec='k', labelsize=16, save=True, fname='crossplot_phases_' + str(p_of_interest),
#                      model=model, verbose=True,
#                      exclude_silica=exclude_silica)


# fo2plt.element_xplot(p_of_interest=p_of_interest, components=['MgO', 'SiO2', 'Al2O3', 'FeO', 'CaO'],
#                      xlim=[(30, 45), (43, 57), (1, 5), (3,10), (1,5)],
#                      y_name='X_Gt', ylabel='Gt abundance (wt%)',
#                      output_parent_path=output_parent_path,
#                      ylim=(0, 16),
#                      linec='k', labelsize=16, save=True, fname='crossplot_gt_' + str(p_of_interest),
#                      model=model, verbose=False,
#                      exclude_silica=exclude_silica)
#
# fo2plt.element_xplot(p_of_interest=p_of_interest, components=['MgO', 'SiO2', 'Al2O3', 'FeO', 'CaO'],
# xlim=[(30, 45), (43, 57), (1, 5), (3,10), (1,5)],
#                      y_name='X_Opx', ylabel='Opx abundance (wt%)',
#                      output_parent_path=output_parent_path,
#                      ylim=(0,80),
#                      linec='k', labelsize=16, save=True, fname='crossplot_opx_' + str(p_of_interest),
#                      model=model, verbose=False,
#                      exclude_silica=exclude_silica)
