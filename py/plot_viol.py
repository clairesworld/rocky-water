import numpy as np
import parameters as p
import perplexdata as px
import plot_perplex as plotpx
import main as rw
import parameters as p
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
from useful_and_bespoke import cornertext, colorize, colourbar, colourised_legend
plt.rc('text', usetex=True)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})  # this only works if outside function for some reason


""" get earth benchmark """
Tp = 1600
earth = rw.build_planet(M_p=1 * p.M_E, test_CMF=0.325, test_oxides=px.wt_oxides_Earth,
                        maxIter=30, tol=1e-4, n=800, Tp=Tp,  # core_efficiency=0.8,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Earth_' + str(Tp) + 'K',
                        )
earth_5m = rw.build_planet(M_p=5 * p.M_E, test_oxides=px.wt_oxides_Earth,
                        maxIter=30, tol=1e-4, n=800, Tp=Tp,  core_efficiency=0.88,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Earth_5M_' + str(Tp) + 'K',
                        )

""" upper mantle water mass - violin plot """
def viol_saturation_mass(which='um', masses=None, yscale=p.TO ** -1, xlabel=None, ylabel=None,
                         labelsize=18, legsize=14, ticksize=14, yticks=None, save=True, xlim='default', ylim='default', fig=None, ax=None,
                         fig_path=plotpx.fig_path, extension='.png', show_legend=True, earth=None):
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
    if masses is None:
        masses = [0.1, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5]
    if xlim == 'default':
        xlim = (min(masses) - 0.2, max(masses) + 0.2)
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
        points = 200
    elif which == 'total':
        key = 'mass_h2o_total'
        textstr = 'Whole mantle'
        fname = 'viol-Mp_w_tot'
        fc = 'xkcd:blood orange'
        if ylim == 'default':
            ylim = (0, 15)
        points = 500  # because extremes look weird

    if ax is None:
        fig, ax = plt.subplots(1, 1)
    mearth_med = None
    for ii, Mp in enumerate(masses):
        if isinstance(Mp, float):
            mass_str = str(Mp).replace('.', ',')
        elif isinstance(Mp, int):
            mass_str = str(Mp)
        else:
            print('MP is', type(Mp))
        planets = rw.read_dir(px.perplex_path_default + 'output/hypatia' + mass_str + 'M_1600K_88Fe/')

        y = []
        for pl in planets:
            try:
                y.append(eval('pl.' + key) * yscale)
            except AttributeError:
                pass  # blank data, didn't work because star not measured maybe

        # plot violin at x=mass
        parts = ax.violinplot(y, positions=[Mp], showmedians=False, showextrema=True, widths=0.25, points=points)

        # edit colours
        if Mp <= 3:
            for pc in parts['bodies']:
                pc.set_alpha(None)  # needed to fix mpl bug
                pc.set_facecolor(mcolors.to_rgba(fc, alpha=1))
                pc.set_edgecolor((0, 0, 0, 0.1))  # 4th component for alpha
                pc.set_linewidth(0.1)
            for prop in ['cmins', 'cmaxes']:
                parts[prop].set_alpha(None)
                parts[prop].set_edgecolor((0, 0, 0, 1))
                parts[prop].set_linewidth(0.7)
            for prop in ['cbars']:
                parts[prop].set_alpha(None)
                parts[prop].set_edgecolor((0, 0, 0, 0.5))
                parts[prop].set_linewidth(0.7)
        else:  # above 3 ME get ppv dissociation and it gets weird
            for pc in parts['bodies']:
                pc.set_alpha(None)
                pc.set_facecolor(mcolors.to_rgba(fc, alpha=0.2))
                pc.set_edgecolor((0, 0, 0, 0.5))  # 4th component for alpha
                pc.set_linewidth(0.5)
                pc.set_linestyle((0, (5, 10)))
            for prop in ['cmins', 'cmaxes', 'cbars']:
                parts[prop].set_alpha(None)
                parts[prop].set_edgecolor((0, 0, 0, 1))
                parts[prop].set_linewidth(0.7)
            for prop in ['cbars']:
                parts[prop].set_alpha(None)
                parts[prop].set_edgecolor((0, 0, 0, 0.3))
        print(Mp, 'fc:', [pc.get_fc() for pc in parts['bodies']])

        # show median
        quartile1, median, quartile3 = np.percentile(y, [25, 50, 75])  # [2.27, 50, 97.73]
        print(Mp, 'M_E', 'med', median, '| range', min(y), max(y), '| sd', np.std(y), ' IQR', quartile3 - quartile1)
        ax.scatter(Mp, median, marker='o', color='0.5', s=30, zorder=3)
        if Mp == 1:
            mearth_med = median

    # add 1:1 line
    if earth is not None:
        y_earth = eval('earth.' + key) * yscale
        ax.scatter(1, y_earth, label='Earth', c='xkcd:dark blue', marker=r'$\oplus$', s=150,
                   zorder=200)  # 200 for old font

        if which == 'total':
            ax.plot([min(masses), max(masses)], [min(masses) * mearth_med, max(masses) * mearth_med], c='k', ls=':', lw=1.5,
                    alpha=0.5, label='Mass scaling')

            # try scaling with area (m_w \propto g**-1) as in Cowan & Abbot
            mg = [0.1, 0.3, 0.5, 0.7, 1, 1.5, 2, 2.5, 3, 4, 5]
            gravs = [3.9167161840595828, 6.083229307426359, 7.511612991847443, 8.637912586721864, 9.8,
                     12.041850104921483, 13.704242549268521, 15.180141341244289, 16.525783688468746, 18.952818930786943,
                     21.13518499447704]
            ax.plot(mg, [mearth_med / (g / 9.8) * m for g, m in zip(gravs, mg)], c='k', ls='-.', lw=1, alpha=0.7,
                    label='Area scaling')


            # add 5M
            # y_earth5 = eval('earth_5m.' + key) * yscale
            # ax.scatter(5, y_earth5, label='Earth5', c='xkcd:cyan', marker=r'$\oplus$', s=150,
            #            zorder=200)  # 200 for old font
            # ax.plot([min(masses), max(masses)], [min(masses)/5 * y_earth5, max(masses)/5 * y_earth5], c='xkcd:cyan', ls=':',
            #         lw=1.5,alpha=0.5, label='5M Mass scaling')

    # limits and legends
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.set_xticks(masses)
    if yticks is not None:
        ax.set_yticks(yticks)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.yaxis.set_ticks_position('both')
    minor_locator = AutoMinorLocator(2)
    ax.yaxis.set_minor_locator(minor_locator)
    if show_legend:
        ax.legend(frameon=False, fontsize=legsize, loc='lower right')
    ax.text(0.03, 0.95, textstr, transform=ax.transAxes, fontsize=labelsize,
            va='top', ha='left', c='0.35', bbox=dict(boxstyle='round, pad=0.3', fc='0.97', ec='0.35', alpha=1, lw=0.5))
    if save:
        plt.tight_layout()
        plt.savefig(fig_path + fname + extension, bbox_inches='tight')
    return fig, ax


def make_subplot(save=True, ylabel='Water capacity (OM)', labelsize=18, fig_path=plotpx.fig_path, extension='.png', fname='viol-Mp', **kwargs):
    fig, axes = plt.subplots(2, 1, figsize=(6, 8), sharex=True)
    fig, axes[0] = viol_saturation_mass(which='um', earth=earth, fig=fig, ax=axes[0], show_legend=False, save=False, xlabel='',
                                   ylabel='', yticks=(0, 1, 2, 3, 4), labelsize=labelsize, **kwargs)
    fig, axes[1] = viol_saturation_mass(which='total', earth=earth, fig=fig, ax=axes[1], show_legend=True, save=False,
                                   ylabel='', labelsize=labelsize, **kwargs)

    fig.supylabel(ylabel, fontsize=labelsize, x=0.05)
    plt.subplots_adjust(hspace=0.1)
    if save:
        plt.tight_layout()
        plt.savefig(fig_path + fname + extension, bbox_inches='tight', dpi=300)
    return fig, axes


make_subplot(save=True, fname='viol-Mp')

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
