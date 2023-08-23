import numpy as np
import parameters as p
import perplexdata as px
import plot_perplex as plotpx
import main as rw
import parameters as p
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as mcolors
from useful_and_bespoke import dark_background, cornertext
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
from useful_and_bespoke import cornertext, colorize, colourbar, colourised_legend
plt.rc('text', usetex=True)
# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})  # this only works if outside function for some reason


""" get earth benchmark """
Tp = 1600
earth = rw.build_planet(M_p=1 * p.M_E, test_CMF=0.325, test_oxides=px.wt_oxides_MD95,
                        maxIter=30, tol=1e-4, n=1200,
                        Tp=Tp,  # core_efficiency=0.8,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Earth_' + str(Tp) + 'K',
                        )
sun = rw.build_planet(M_p=1 * p.M_E,star='sun',
                        maxIter=30, tol=1e-4, n=1200,
                        Tp=Tp,  core_efficiency=0.88,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Sun_' + str(Tp) + 'K',
                        )
# earth_5m = rw.build_planet(M_p=5 * p.M_E, test_oxides=px.wt_oxides_Earth,
#                         maxIter=30, tol=1e-4, n=800, Tp=Tp,  core_efficiency=0.88,
#                         plot_all=False, get_saturation=True, verbose=True, clean=True,
#                         vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                         name='Earth_5M_' + str(Tp) + 'K',
#                         )

""" upper mantle water mass - violin plot """
def viol_saturation_mass(which='um', masses=None, yscale=p.TO ** -1, xlabel=None, ylabel=None,
                         labelsize=18, legsize=14, ticksize=14, yticks=None, save=True, xlim='default', ylim='default', fig=None, ax=None,
                         fig_path=plotpx.fig_path, extension='.png', show_legend=True, earth=None, sun=None, exclude=True):
    plt.rc('text', usetex=True)
    # plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
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
        fc = 'xkcd:sage'
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
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    mearth_med = None
    for ii, Mp in enumerate(masses):
        if isinstance(Mp, float):
            mass_str = str(Mp).replace('.', ',')
        elif isinstance(Mp, int):
            mass_str = str(Mp)
        else:
            print('MP is', type(Mp))
        planets = rw.read_dir(px.perplex_path_default + 'output/apollo/hypatia' + mass_str + 'M_1600K_88Fe_hires/')

        y = []
        for pl in planets:
            try:
                y.append(eval('pl.' + key) * yscale)
            except AttributeError:
                pass  # blank data, didn't work because star not measured maybe

        # exclude extrema?
        if exclude:
            y_min, y_max = np.percentile(y, [2.27, 97.93], method='nearest')
            y = [a for a in y if a <= y_max]
            y = [a for a in y if a >= y_min]

        # plot violin at x=mass
        parts = ax.violinplot(y, positions=[Mp], showmedians=False, showextrema=True, widths=0.2, points=points)

        # edit colours
        if Mp <= 3 or which=='um':
            for pc in parts['bodies']:
                pc.set_alpha(None)  # needed to fix mpl bug
                pc.set_facecolor(mcolors.to_rgba(fc, alpha=0.5))
                pc.set_edgecolor((1, 1, 1, 0.1))  # 4th component for alpha
                pc.set_linewidth(2)
            for prop in ['cmins', 'cmaxes']:
                parts[prop].set_alpha(None)
                parts[prop].set_edgecolor((1, 1, 1, 1))
                parts[prop].set_linewidth(2)
            for prop in ['cbars']:
                parts[prop].set_alpha(None)
                parts[prop].set_edgecolor((1, 1, 1, 0.5))
                parts[prop].set_linewidth(2)
                parts[prop].set_zorder(2)
        else:  # above 3 ME get ppv dissociation and it gets weird
            for pc in parts['bodies']:
                pc.set_alpha(None)
                pc.set_facecolor(mcolors.to_rgba(fc, alpha=0.2))
                pc.set_edgecolor((1, 1, 1, 0.5))  # 4th component for alpha
                pc.set_linewidth(2)
                pc.set_linestyle((0, (5, 10)))
            for prop in ['cmins', 'cmaxes', 'cbars']:
                parts[prop].set_alpha(None)
                parts[prop].set_edgecolor((1, 1, 1, 1))
                parts[prop].set_linewidth(2)
            for prop in ['cbars']:
                parts[prop].set_alpha(None)
                parts[prop].set_edgecolor((1, 1, 1, 0.3))

        # show median
        quartile1, median, quartile3 = np.percentile(y, [25, 50, 75])  # [2.27, 50, 97.73]
        print(Mp, 'M_E', 'med', median, '| range', min(y), max(y), '| sd', np.std(y), ' IQR', quartile3 - quartile1,
              '| quartiles', quartile1, quartile3)
        ax.scatter(Mp, median, marker='o', color='0.5', s=30, zorder=3)
        if Mp == 1:
            mearth_med = median

    # add 1:1 line
    if sun:
        y_sun = eval('sun.' + key) * yscale
        ax.scatter(1, y_sun, label='Solar', c='#b13c02', #alpha=0.4,
                   marker='*', s=200, zorder=100)
    if earth is not None:
        y_earth = eval('earth.' + key) * yscale
        ax.scatter(1, y_earth, label='Earth', c='#b13c02', #'xkcd:dark blue',
                   marker=r'$\oplus$', s=380,
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
    ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=15)
    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=15)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.yaxis.set_ticks_position('both')
    # minor_locator = AutoMinorLocator(2)
    # ax.yaxis.set_minor_locator(minor_locator)
    if show_legend:
        ax.legend(frameon=False, fontsize=legsize, loc='lower left', bbox_to_anchor=(1.01, 0.001))
    ax.text(0.03, 0.95, textstr, transform=ax.transAxes, fontsize=labelsize, zorder=500,
            va='top', ha='left', c='xkcd:off white', #bbox=dict(boxstyle='round, pad=0.3', fc='0.97', ec='0.35', alpha=1, lw=0.5)
            )
    if save:
        plt.tight_layout()
        plt.savefig(fig_path + fname + extension, bbox_inches='tight')
    return fig, ax




fig, ax = viol_saturation_mass(masses = [0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5],
                               which='um', earth=earth, sun=sun,  show_legend=True, save=False,
                                   ylabel='Water capacity (Earth oceans)', labelsize=28, legsize=28, ticksize=18)
fig, ax = dark_background(fig, ax)
plt.savefig('/home/claire/Desktop/viol.png', bbox_inches='tight', dpi=300)
plt.show()
