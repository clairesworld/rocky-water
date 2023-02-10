import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
from py.useful_and_bespoke import cornertext

""" make cross plot of oxides / elements"""

labelsize = 22
core_eff = 88
Xf = 3
exclude_silica = True
ncols = 2

""" plot fo2 vs. element ratios (cols) for two pressures (rows) """
""" oolour by Fe3+ content in opx?"""

fig = plt.figure(figsize=(ncols * 4, 8))
gs = fig.add_gridspec(2, ncols + 1, width_ratios=[10] * ncols + [1], height_ratios=[10, 10],
                      left=0.1, right=0.9,
                      wspace=0.05, hspace=0.05)

for jj, p_of_interest in enumerate((1, 4)):
    axes = []  # create new ax list just for this row
    [axes.append(fig.add_subplot(gs[jj, iii])) for iii in range(ncols)]
    for ii, model in enumerate(['melts', 'perplex']):
        if model == 'melts':
            source = fo2plt.output_parent_mlt_earth
            ylims = [(-11, -8.5), (-8, -5.3)]  # 1 GPa, 4 GPa
            c = 'xkcd:midnight blue'
            title = 'pMELTS'
        elif model == 'perplex':
            source = fo2plt.output_parent_px
            ylims = [(-13.5, -8.7), (-10, -6.5)]
            c = 'xkcd:brick red'
            title = 'Perple_X'
        output_parent_path = source + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/'

        fig, axs = fo2plt.element_xplot(p_of_interest=p_of_interest, components=['Fe/H'],
                                        output_parent_path=output_parent_path,
                                        output_parent_path_melts=fo2plt.output_parent_mlt_earth + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/',
                                        # xlim=[(0.8, 1.65), (0.05, 0.2), (0, 0.18), (0.025, 0.13)],
                                        ylim=ylims[jj], c=c, s=30, alpha=0.4, labelsize=labelsize, save=False, fig=fig,
                                        axes=[axes[ii]],
                                        model=model, verbose=False,
                                        # z_name = 'X_Fe3_Opx', vmin=0, vmax=0.6, cmap='viridis',
                                        exclude_silica=exclude_silica)

        if ii == 0:
            axs[0] = cornertext(axs[0], str(p_of_interest) + ' GPa', size=labelsize, pos='top left')

        elif ii == 1:
            axs[0].set_yticks([])
            axs[0].set_ylabel('')
        # axs[1].set_xticks([0.075, 0.125, 0.175])  # Fe/Si
        if jj == 0:
            for ax in axs:
                ax.set_xticks([])
                ax.set_xlabel(None)
                ax.set_title(title, fontsize=labelsize)
        elif jj == 1:
            for ax in axs:
                ax.set_xlabel('[Fe/H]', fontsize=labelsize)


fig.suptitle('', fontsize=labelsize, y=0.92)
fig.savefig(fo2plt.figpath + 'crossplot_metallicity_' + model + '.png', bbox_inches='tight')
# plt.show()
