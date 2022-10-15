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
from matplotlib import rc

# activate latex text rendering
rc('text', usetex=True)
plt.rc('font', **{'family': 'serif',
                  'serif': ['Computer Modern Roman']})  # this only works if outside function for some reason


""" get earth benchmark """
Tp = 1600
earth = rw.build_planet(M_p=1 * p.M_E, test_CMF=0.325, test_oxides=px.wt_oxides_MD95,
                        maxIter=30, tol=1e-4, n=1200, Tp=Tp,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Earth_' + str(Tp) + 'K',
                        )
sun = rw.build_planet(M_p=1 * p.M_E, star='sun',
                      maxIter=30, tol=1e-4, n=1200,
                      Tp=Tp, core_efficiency=0.88,
                      plot_all=False, get_saturation=True, verbose=True, clean=True,
                      vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                      name='Sun_' + str(Tp) + 'K',
                      )


def load_exoplanets(directory, param, scale=1):
    dats = rw.read_dir(px.perplex_path_default + 'output/' + directory)
    names = []
    masses = []
    ys = []
    for dat in dats:
        # print('mgsi',
        #       10 ** dat.nH_star[dat.oxide_list.index('MgO')] / 10 ** dat.nH_star[dat.oxide_list.index('SiO2')]
        #       # molar ratio
        #       )
        names.append(dat.name)
        masses.append(dat.M_p / p.M_E)
        ys.append(eval('dat.' + param) * scale)
    return names, masses, ys


""" scatter hist of mgsi effect on upper mantle water but over all masses """


def scatterhist_mgsi_allmasses(masses=[0.1, 0.3, 0.5, 1, 2, 3, 4, 5], cmap='rainbow', earth=None, alpha=0.4, ylim=None,
                               labelsize=14, ticksize=12, legsize=12, xlim=(0, 3), vmin=None, vmax=None,
                               show_histy=True,
                               which='total', ms=200, output_path=px.perplex_path_default + 'output/',
                               end='_1600K_88Fe',
                               fname=None):
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family': 'serif',
                      'serif': ['Computer Modern Roman']})  # this only works if outside function for some reason
    ii_1M = masses.index(1)

    if which == 'total':
        y_var = "mass_h2o_total"
        ylabel = 'Water capacity (ppm wt)'
        if fname is None:
            fname = 'mgsi_c_h2o_mantle_scatter_all'
        if ylim is None:
            ylim = (100, 2000)  # (0, 0.3),
    elif which == 'um':
        y_var = "mass_h2o_um"
        ylabel = 'UM water storage capacity (ppm)'
        if fname is None:
            fname = 'mgsi_c_h2o_um_scatter_all'
        if ylim is None:
            ylim = (0, 2000)
    mass_str = []
    for Mp in masses:
        if isinstance(Mp, float):  # Tp = 1600
            mass_str.append(str(Mp).replace('.', ','))
        elif isinstance(Mp, int):
            mass_str.append(str(Mp))
        else:
            print('MP is', type(Mp))
    dirs = ['hypatia' + ms + 'M' + end for ms in mass_str]
    colours = colorize(np.arange(len(dirs)), cmap=cmap, vmin=vmin, vmax=vmax)[0]

    # get colour for 1ME
    c_1M_alpha = colours[ii_1M]
    c_1M_alpha[3] = 0.8  # 0.25  # set alpha

    # set up scatter & histogram gridspec
    dats = rw.read_dir(output_path + dirs[ii_1M] + '/')
    fig, ax, ax_histx, ax_histy = plotpx.pop_scatterhist(dats, "mgsi", y_var, y_scale=1e6 / (masses[ii_1M] * p.M_E),
                                                         xlabel='Mg/Si', ylabel='', show_histx=True,
                                                         show_histy=show_histy,
                                                         histx_kwargs={'color': '0.9', 'edgecolor': 'k',
                                                                       'linewidth': 0.5, 'bins': 55},
                                                         histy_kwargs={'color': c_1M_alpha, 'edgecolor': colours[ii_1M],
                                                                       'linewidth': 0.9, 'bins': 40},
                                                         earth=False, xlim=xlim, ylim=ylim, ms=ms,
                                                         save=False, show=False, lim_histx=(0, 6), lim_histy=(0, 0.01),
                                                         annotate_n=False, labelsize=labelsize, legsize=legsize,
                                                         c=colours[ii_1M], alpha=alpha)
    # show Earth vline
    for axx in [ax, ax_histx]:
        axx.axvline(earth.mgsi, c='k', alpha=0.2, lw=0.5, ls=(10, (15, 10)))
    for axy in [ax, ax_histy]:
        axy.axhline(eval('earth.' + y_var) * 1e6 / (masses[ii_1M] * p.M_E),
                    c='k', alpha=0.2, lw=0.5, ls=(10, (15, 10)))

    # overplot other masses
    for ii, d in enumerate(dirs):
        if ii != ii_1M:
            dats = rw.read_dir(output_path + d + '/', verbose=True)
            fig, ax = plotpx.pop_scatter(dats, "mgsi", y_var, y_scale=1e6 / (masses[ii] * p.M_E), xlabel='Mg/Si',
                                         ylabel=ylabel, annotate_n=False,
                                         earth=False, fig=fig, ax=ax, xlim=xlim, ms=ms,
                                         ylim=ylim, labelsize=labelsize, ticksize=ticksize, legsize=legsize, save=False,
                                         show=False, c=colours[ii], alpha=alpha)

    # show Earth scaled Mg/Si
    dats_e = rw.read_dir(px.perplex_path_default + 'output/MgSi_from_earth/')
    x_e, y_e = [], []
    for dat in dats_e:
        x_e.append(dat.mgsi)
        y_e.append(eval('dat.' + y_var) / p.M_E * 1e6)
    x_e, y_e = zip(*sorted(zip(x_e, y_e)))
    ax.plot(x_e, y_e, c='k', ls=':', label='Compositionally-scaled Earth')
    ax.scatter(sun.mgsi, eval('sun.' + y_var) / sun.M_p * 1e6, label='Solar', c='k',  # alpha=0.4,
               marker='*', s=100, zorder=100)
    ax.scatter(earth.mgsi, eval('earth.' + y_var) / p.M_E * 1e6, label='Earth', c='xkcd:dark blue',
               marker='$\oplus$', s=200, zorder=200)

    # make legends
    if show_histy:
        leg_bbox_to_anchor = (1.01, 1)
    else:
        leg_bbox_to_anchor = (1.01, 1)
    cax = colourised_legend(ax, colours, [str(m) + ' $M_\oplus$' for m in masses],  # alpha=alpha,
                            title=None,  # r'$\bf{Planet}$ $\bf{mass}$',
                            legsize=legsize, titlesize=legsize, markersize=5,
                            bbox_to_anchor=leg_bbox_to_anchor, handletextpad=0.1)
    ax_histy.set_zorder(-1)  # to go below leg
    ax.legend(handles=[mlines.Line2D([], [], color='k', marker='*',
                                     markersize=8, lw=0, label='Solar'),
                       mlines.Line2D([], [], color='xkcd:dark blue', marker='$\oplus$',
                                     markersize=10, lw=0, label='Earth'),
                       mlines.Line2D([], [], color='k', marker=None,
                                     lw=1, ls=':', label='Compositionally-scaled Earth'),
                       ], frameon=False, fontsize=legsize - 1)
    fig.savefig(plotpx.fig_path + fname + '.png', bbox_inches='tight', dpi=300)



def add_exo_axis(ax, xlim=None, cmap=None, vmin=None, vmax=None, y_var=None, labelsize=None, ticksize=None,
                 legsize=None):
    exo_names, exo_masses, exo_mgsi = load_exoplanets('earthsize_planets_1600K_88Fe/', 'mgsi', scale=1)
    _, _, exo_w = load_exoplanets('earthsize_planets_1600K_88Fe/', y_var, scale=1)

    exo_colours = colorize(exo_masses, cmap=cmap, vmin=vmin, vmax=vmax)[0]

    ax.set_ylim([250, 700])
    ax.set_xlim(xlim)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.set_xlabel('Mg/Si', fontsize=labelsize)

    # errorbars
    exo_names_min, _, exo_ys_min = load_exoplanets('earthsize_planets_1600K_88Fe_mgsimin/', y_var, scale=1)
    exo_names_max, _, exo_ys_max = load_exoplanets('earthsize_planets_1600K_88Fe_mgsimax/', y_var, scale=1)
    _, _, exo_xs_min = load_exoplanets('earthsize_planets_1600K_88Fe_mgsimin/', 'mgsi', scale=1)
    _, _, exo_xs_max = load_exoplanets('earthsize_planets_1600K_88Fe_mgsimax/', 'mgsi', scale=1)

    medy, minw, maxw = [], [], []
    for ii, pl in enumerate(exo_names):
        print('name', pl, '=', exo_names_min[ii], '=', exo_names_max[ii])
        medy.append(exo_w[ii] / (exo_masses[ii] * p.M_E) * 1e6)
        minw.append(exo_ys_min[ii] / (exo_masses[ii] * p.M_E) * 1e6)
        maxw.append(exo_ys_max[ii] / (exo_masses[ii] * p.M_E) * 1e6)
        print(pl, 'wmf =', medy[ii], '_', minw[ii], '^', maxw[ii], 'ppm')
    yerr_max = (np.array(maxw) - np.array(medy))
    yerr_min = (np.array(medy) - np.array(minw))
    yerr = np.vstack((yerr_min, yerr_max))
    # print('yerr', np.shape(yerr))

    xerr_max = (np.array(exo_xs_max) - np.array(exo_mgsi))
    xerr_min = (np.array(exo_mgsi) - np.array(exo_xs_min))
    xerr = np.vstack((xerr_min, xerr_max))

    for jj, xx in enumerate(exo_mgsi):
        # print('yerr[: , jj]', np.shape(yerr[:, jj].T), yerr[:, jj])
        print(['xerr', [xerr[0, jj]], [xerr[1, jj]]])
        print(['yerr', [yerr[0, jj]], [yerr[1, jj]]])
        ax.errorbar(xx, exo_w[jj] / (exo_masses[jj] * p.M_E) * 1e6,
                    xerr=[[xerr[0, jj]], [xerr[1, jj]]],
                    yerr=[[yerr[0, jj]], [yerr[1, jj]]],
                    elinewidth=1, ecolor=exo_colours[jj], alpha=0.6,
                    marker='.', ms=0, zorder=99,  # cmap=cmap,
                    # vmin=0.1, vmax=5
                    )
        ax.scatter(exo_mgsi[jj], exo_w[jj] / (exo_masses[jj] * p.M_E) * 1e6,
                   c=exo_colours[jj],  # alpha=0.4,
                   marker='.', s=100, zorder=100,  # cmap=cmap, vmin=vmin, vmax=vmax
                   )

    # annotate exoplanets
    for en, em, emg, ew in zip(exo_names, exo_masses, exo_mgsi, exo_w):
        print(en, em, 'ME, mg/si =', emg, 'wmf', ew / (em * p.M_E) * 1e6, 'ppm')

        if en in ['Kepler-128c', 'Kepler-11c']:  # label subset  , 'Kepler-51b'
            ax.text(emg + 0.02, ew / (em * p.M_E) * 1e6 - 10, en, fontsize=legsize, va='top', ha='left', zorder=110)
        # elif en in ['Kepler-11b']:  # label subset
        #     ax.text(emg + 0.02, ew / (em * p.M_E) * 1e6 + 2, en, fontsize=legsize, va='bottom', ha='left')
        elif en in ['Kepler-167d']:  # label subset  'Kepler-68c'
            ax.text(emg - 0.02, ew / (em * p.M_E) * 1e6 - 2, en, fontsize=legsize, va='top', ha='right', zorder=110)
        elif en in ['Kepler-128b', 'Kepler-78b']:  # label subset
            ax.text(emg - 0.01, ew / (em * p.M_E) * 1e6 + 2, en, fontsize=legsize, va='bottom', ha='right', zorder=110)
        return ax


def scatterhist_mgsi_allmasses_kepler(masses=[0.1, 0.3, 0.5, 1, 2, 3, 4, 5], cmap='rainbow', earth=None, sun=None,
                                      alpha=0.4, ylim=None, ylims=None, fformat='.png',
                                      labelsize=14, ticksize=12, legsize=12, xlim=(0, 3), vmin=None, vmax=None,
                                      ms=200, ratios=[7, 1], output_path=px.perplex_path_default + 'output/', **kwargs):
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family': 'serif',
                      'serif': ['Computer Modern Roman']})  # this only works if outside function for some reason
    ii_1M = masses.index(1)
    nrows = len(masses) + 1  # each mass + kepler

    y_var = "mass_h2o_total"
    ylabel = 'Water capacity (ppm wt)'
    fname = 'mgsi_c_h2o_exo_2'
    if ylim is None:
        ylim = (100, 2000)
    mass_str = []
    for Mp in masses:
        if isinstance(Mp, float):  # Tp = 1600
            mass_str.append(str(Mp).replace('.', ','))
        elif isinstance(Mp, int):
            mass_str.append(str(Mp))
        else:
            print('MP is', type(Mp))

    fig = plt.figure(figsize=(8, 8))

    gs = fig.add_gridspec(1 + nrows, 2, width_ratios=[ratios[0] + 2, ratios[1]],
                          height_ratios=[ratios[1] + 1] + [ratios[0]] * len(masses) + [1.5],
                          # change last one for kepler panel
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    ax_histx = fig.add_subplot(gs[0, 0])
    axes = []
    axes_histy = []
    for ii in range(1, nrows):
        axes.append(fig.add_subplot(gs[ii, 0], sharex=ax_histx))  # model
        axes_histy.append(fig.add_subplot(gs[ii, 1], sharey=axes[ii - 1]))

    dirs = ['hypatia' + ms + 'M_1600K_88Fe_hires' for ms in mass_str]
    cmfs = np.linspace(vmin, vmax, num=7)

    colours = colorize(cmfs, cmap=cmap, vmin=vmin, vmax=vmax)[0]

    # get colour for 1ME
    c_1M_alpha = colours[ii_1M]
    c_1M_alpha[3] = 0.8  # 0.25  # set alpha

    # set up scatter & histogram gridspec
    for ii, d in enumerate(dirs):
        if ylims is not None:
            ylim = ylims[ii]
        dats = rw.read_dir(output_path + d + '/', verbose=True)
        fig, axes[ii], ax_histx, axes_histy[ii] = plotpx.pop_scatterhist(dats, "mgsi", y_var,
                                                                         y_scale=1e6 / (masses[ii] * p.M_E),
                                                                         xlabel='', ylabel='',
                                                                         histx_kwargs={'color': '0.9', 'edgecolor': 'k',
                                                                                       'linewidth': 0.5, 'bins': 40},
                                                                         histy_kwargs={'color': c_1M_alpha,
                                                                                       'edgecolor': 'k',  # colours[ii],
                                                                                       'linewidth': 0.9, 'bins': 40},
                                                                         earth=False, xlim=xlim, ylim=ylim, ms=ms,
                                                                         save=False, show=False, lim_histx=(0, 4),
                                                                         lim_histy=(0, 0.02),
                                                                         annotate_n=False, labelsize=labelsize,
                                                                         legsize=legsize,
                                                                         # c=colours[ii],
                                                                         alpha=alpha,
                                                                         fig=fig, ax=axes[ii], ax_histy=axes_histy[ii],
                                                                         ax_histx=ax_histx,
                                                                         show_histx=True, show_histy=True,
                                                                         c='CMF', cmap=cmap, vmin=vmin, vmax=vmax)
        axes[ii] = cornertext(axes[ii], text=str(masses[ii]) + r' $M_E$', pos='top left', )

    # # show real planets
    # axes.append(fig.add_subplot(gs[nrows, 0]))  # exo
    # axes[-1] = add_exo_axis(axes[-1], xlim=xlim, cmap=cmap, vmin=vmin, vmax=vmax, y_var=y_var, labelsize=labelsize, ticksize=ticksize,
    #            legsize=legsize)
    # for axx in axes[:-1]:
    #     axx.set_xticks([])

    # show Earth scaled Mg/Si
    dats_e = rw.read_dir(px.perplex_path_default + 'output/MgSi_from_earth/')
    x_e, y_e = [], []
    for dat in dats_e:
        x_e.append(dat.mgsi)
        y_e.append(eval('dat.' + y_var) / p.M_E * 1e6)
    x_e, y_e = zip(*sorted(zip(x_e, y_e)))
    for ax in axes:
        ax.plot(x_e, y_e, c='k', ls=':', label='Compositionally-scaled Earth', alpha=0.5)
        ax.scatter(sun.mgsi, eval('sun.' + y_var) / sun.M_p * 1e6, label='Solar', c='k',  # alpha=0.4,
                   marker='*', s=100, zorder=100)
        ax.scatter(earth.mgsi, eval('earth.' + y_var) / p.M_E * 1e6, label='Earth', c='xkcd:dark blue',
                   marker='$\oplus$', s=200, zorder=200)

    # show Earth vline
    for axx in axes + [ax_histx]:
        axx.axvline(earth.mgsi, c='k', alpha=0.2, lw=0.5, ls=(10, (15, 10)))
    for axy in axes + axes_histy:
        axy.axhline(eval('earth.' + y_var) * 1e6 / (masses[ii_1M] * p.M_E),
                    c='k', alpha=0.2, lw=0.5, ls=(10, (15, 10)))

    # make legends
    leg_bbox_to_anchor = (1.05, 1)
    cax = colourised_legend(axes[0], colours, ['{:.2f}'.format(cmf) for cmf in cmfs],  # alpha=alpha,
                            title=r'$\bf{CMF}$',
                            legsize=legsize, titlesize=legsize, markersize=5,
                            bbox_to_anchor=leg_bbox_to_anchor, handletextpad=0.1)

    # cbar = colourbar(mappable=None, vector=np.linspace(vmin, vmax), ax=axes[0], vmin=vmin, vmax=vmax, label='CMF',
    #                  labelsize=14, ticksize=14,
    #           ticks=None, ticklabels=None, labelpad=17, loc='right', cax=None,
    #           rot=None, discrete=False, cmap=cmap, tickformatter=None, pad=0.05,  **kwargs)

    for ax_histy in axes_histy:
        ax_histy.set_zorder(-1)  # to go below leg
    axes[0].legend(handles=[mlines.Line2D([], [], color='k', marker='*',
                                          markersize=8, lw=0, label='Solar'),
                            mlines.Line2D([], [], color='xkcd:dark blue', marker='$\oplus$',
                                          markersize=10, lw=0, label='Earth'),
                            mlines.Line2D([], [], color='k', marker=None,
                                          lw=1, ls=':', label='Compositionally-scaled Earth'),
                            ], frameon=False, fontsize=legsize, bbox_to_anchor=(1.05, -1), loc='upper left')

    fig.supylabel(ylabel, fontsize=labelsize, y=0.45)
    fig.savefig(plotpx.fig_path + fname + fformat, bbox_inches='tight', dpi=200)


masses = [0.1, 0.5, 1, 2, 3]


scatterhist_mgsi_allmasses_kepler(masses=masses, cmap='viridis', earth=earth, sun=sun, alpha=0.2, vmin=0.2, vmax=0.35,
                                  xlim=(0.7, 1.6), ms=20,
                                  ylims=[(210, 1700), (210, 1300), (210, 1000), (210, 800), (210, 550)],
                                  output_path=px.perplex_path_default + 'output/apollo/', fformat='.png')


# isolate 3 M_E - what is causing scatter here?
def spread_test(mass=3, cmap='viridis', alpha=0.4, ylim=(200, 600), xlims=None, fformat='.png',
                labelsize=14, ticksize=12, legsize=12, xlim=None, vmin=None, vmax=None, suptitle=r'3 $M_\oplus$',
                ms=50, output_path=px.perplex_path_default + 'output/apollo/', fig=None, axes=None, save=True,
                xvars=['mgsi', 'CMF', 'wt_oxides["FeO"]', 'nH_star[-1]'],
                xlabels=['Mg/Si', 'CMF', r'Mantle FeO (wt\%)', '[Fe/H]$_*$'],
                y_var="mass_h2o_total", ylabel='Water capacity (ppm wt)',
                **kwargs):

    if isinstance(mass, float):  # Tp = 1600
        mass_str = str(mass).replace('.', ',')
    elif isinstance(mass, int):
        mass_str = str(mass)
    else:
        print('MP is', type(mass))
    dats = rw.read_dir(output_path + 'hypatia' + mass_str + 'M_1600K_88Fe_hires/', verbose=True)

    if axes is None:
        fig, axes = plt.subplots(1, len(xvars), figsize=(10, 3))

    for ii, var in enumerate(xvars):
        if xlims is not None:
            xlim = xlims[ii]
        fig, axes[ii] = plotpx.pop_scatter(dats, var, y_var,
                                           y_scale=1e6 / (mass * p.M_E),
                                           xlabel=xlabels[ii], ylabel=ylabel,
                                           earth=False, xlim=xlim, ylim=ylim, ms=ms,
                                           save=False, show=False,
                                           annotate_n=False, labelsize=labelsize,
                                           legsize=legsize,
                                           # c=colours[ii],
                                           alpha=alpha, fig=fig, ax=axes[ii],
                                           # c='CMF', cmap=cmap, vmin=vmin, vmax=vmax,
                                           **kwargs)

    for ax in axes[1:]:
        ax.set_ylabel('')
        ax.axes.yaxis.set_ticklabels([])

    if suptitle:
        fig.suptitle(suptitle, fontsize=20)

    if save:
        plt.tight_layout()
        fig.savefig(plotpx.fig_path + '3ME_spread.png', bbox_inches='tight')
    return fig, axes


def spread_test_subplot(masses, xvars=None, ylims=[(300, 1500), (200, 1000), (200, 700)], labelsize=14, fformat='.pdf',
                        **kwargs):
    if xvars is None:
        xvars = ['mgsi', 'CMF', 'wt_oxides["FeO"]', 'nH_star[-1]']

    fig, axes = plt.subplots(len(masses), len(xvars), figsize=(10, 7))

    for ii, mass in enumerate(masses):
        fig, axes[ii] = spread_test(mass=masses[ii], fig=fig, axes=axes[ii], suptitle=None, labelsize=labelsize,
                                    ylim=ylims[ii], ylabel='', **kwargs)

        axes[ii][0] = cornertext(axes[ii][0], str(mass) + ' $M_\oplus$', pos='top left', size=labelsize)

    # label bottom row
    for axe in axes[:-1][:]:
        for ax in axe:
            ax.set_xlabel('')
            ax.axes.xaxis.set_ticklabels([])

    fig.supylabel('Water capacity mass fraction (ppm)', fontsize=labelsize)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1, wspace=0.11)
    fig.savefig(plotpx.fig_path + 'masses_spread' + fformat, bbox_inches='tight')


# spread_test_subplot(masses=[0.1, 1, 3], xlims=[(0.8, 1.75), (0.1, 0.45), (2, 12), None], ms=10, alpha=0.2, labelsize=16,
#                     xlabelpad=9, show_corr=True, min_mgsi=0.85, c='xkcd:grass', transparent_edge=True)

plt.show()
