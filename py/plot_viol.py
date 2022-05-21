import numpy as np
import parameters as p
import perplexdata as px
import plot_perplex as plotpx
import main as rw
import parameters as p
import matplotlib.pyplot as plt
from useful_and_bespoke import cornertext, colorize, colourbar, colourised_legend
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec

""" get earth benchmark """
Tp = 1600
earth = rw.build_planet(M_p=1 * p.M_E, test_CMF=0.325, test_oxides=px.wt_oxides_Earth,
                        maxIter=30, tol=1e-4, n=800, Tp=Tp,  # core_efficiency=0.8,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Earth_' + str(Tp) + 'K',
                        )


""" upper mantle water mass - violin plot """


def viol_saturation_mass(which='um', masses=None, yscale=p.TO ** -1, xlabel=None, ylabel=None,
                         labelsize=14, legsize=12, ticksize=12, yticks=None,
                         cmap='rainbow', save=True, xlim='default', ylim='default', fig=None, ax=None,
                         fig_path=plotpx.fig_path, extension='.png', show_legend=True, earth=None):
    if masses is None:
        masses = [0.1, 0.5, 1, 2, 3, 4, 5]
    if xlim == 'default':
        xlim = (-0.25, max(masses) + 0.5)
    if xlabel is None:
        xlabel = 'Planet mass ($M_\oplus$)'
    if ylabel is None:
        ylabel = 'Water capacity (OM)'
    if which == 'um':
        key = 'mass_h2o_um'
        textstr = 'Upper mantle'
        fname = 'viol-Mp_w_um'
        fc = 'xkcd:seafoam'
        if ylim == 'default':
            ylim = (0, 4)
    elif which == 'total':
        key = 'mass_h2o_total'
        textstr = 'Whole mantle'
        fname = 'viol-Mp_w_tot'
        fc = 'xkcd:blood orange'
        if ylim == 'default':
            ylim = (0, 15)

    if ax is None:
        fig, ax = plt.subplots(1, 1)
    c = colorize(masses, cmap)[0]
    for ii, Mp in enumerate(masses):
        if isinstance(Mp, float):
            mass_str = str(Mp).replace('.', ',')
        elif isinstance(Mp, int):
            mass_str = str(Mp)
        else:
            print('MP is', type(Mp))
        planets = rw.read_dir(px.perplex_path_default + 'output/hypatia' + mass_str + 'M/')

        y = []
        for pl in planets:
            try:
                y.append(eval('pl.' + key) * yscale)
            except AttributeError:
                pass  # blank data, didn't work because star not measured maybe

        # plot violin at x=mass
        parts = ax.violinplot(y, positions=[Mp], showmedians=False, showextrema=True, widths=0.25)

        # edit colours
        for pc in parts['bodies']:
            pc.set_facecolor(fc)
            pc.set_edgecolor('k')
            pc.set_alpha(1)
        for prop in ['cmins', 'cmaxes', 'cbars']:
            parts[prop].set_edgecolor('k')

        # show median
        quartile1, medians, quartile3 = np.percentile(y, [25, 50, 75])
        print(Mp, 'M_E', 'med', medians, '| range', min(y), max(y))
        ax.scatter(Mp, medians, marker='o', color='k', s=30, zorder=3)

    # add 1:1 line
    if earth is not None:
        y_earth = eval('earth.' + key) * yscale
        ax.scatter(1, y_earth, label='Earth', c='xkcd:dark blue', marker=r'$\oplus$', s=200, zorder=100)
        ax.plot([min(masses), max(masses)], [min(masses) * y_earth, max(masses) * y_earth], c='k', ls=':', lw=1,
                label='1:1 scaling')

    # limits and legends
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.set_xticks(masses)
    if yticks is not None:
        ax.set_yticks(yticks)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    if show_legend:
        ax.legend(frameon=False, fontsize=legsize, loc='upper left')
    ax.text(0.97, 0.03, textstr, transform=ax.transAxes, fontsize=legsize,
            va='bottom', ha='right', c='0.35', bbox=dict(boxstyle='round', fc='w', ec='0.35', alpha=1))
    if save:
        plt.tight_layout()
        plt.savefig(fig_path + fname + extension, bbox_inches='tight')
    return fig, ax


fig, axes = plt.subplots(2, 1, figsize=(6, 8), sharex=True)
fig, ax = viol_saturation_mass(which='um', earth=earth, fig=fig, ax=axes[0], save=False, xlabel='', yticks=(0, 1, 2, 3, 4))
fig, ax = viol_saturation_mass(which='total', earth=earth, fig=fig, ax=axes[1], show_legend=False)

plt.show()

""" examine extrema """
# for name in mins:
#     dat = rw.read_name(output_path=px.output_parent_default + 'hypatia0,1M/', name=name, M_p=0.1, core_efficiency=0.8, Tp=1600)
#     dat.mass_h2o_um = sat.total_water_mass(dat.df_all, i_min=0, i_max=dat.find_lower_mantle() - 1)
#     i_um_base = dat.find_lower_mantle() - 1
#     print('i_um_base', i_um_base, 'len', len(dat.df_all))
#
#     plotpx.composition_subfig(dat, 'c_h2o', var_scale=1e6, var_log=True, var_label='Water capacity\n(ppm)',
#                               vertical_pressure=False, title=dat.star, show_UM=True, p_max=None,
#                               annotation='total = {0:3.1f} earth oceans'.format(dat.mass_h2o_total / p.TO),
#                               save=False, show=True, phase_order=None, labelsize=12, legsize=10, ticksize=12, xpad=10,
#                               ypad=10, linec='xkcd:navy', cmap_phases='tab20', linew=2, ymin=1e0, ymax=1e5)

