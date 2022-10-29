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


# def scatterhist_mgsi_allmasses_multipanel(masses=[0.1, 0.3, 0.5, 1, 2, 3, 4, 5], cmap='rainbow', earth=None, alpha=0.4,
#                                           ylims=None,
#                                           labelsize=14, ticksize=12, legsize=12, xlim=(0, 3), vmin=None, vmax=None,
#                                           show_histy=True,
#                                           ms=200):
#     plt.rc('text', usetex=True)
#     plt.rc('font', **{'family': 'serif',
#                       'serif': ['Computer Modern Roman']})  # this only works if outside function for some reason
#     ii_1M = masses.index(1)
#
#     y_vars = ["mass_h2o_total", "c_h2o_obm"]
#     ylabels = ['Total H$_2$O storage capacity (ppm)', 'Olivine-bearing mantle H$_2$O content (ppm)']
#     fname = 'mgsi_c_h2o_obm'
#     if ylims is None:
#         ylims = [(100, 2000), (0, 450)]
#     mass_str = []
#     for Mp in masses:
#         if isinstance(Mp, float):  # Tp = 1600
#             mass_str.append(str(Mp).replace('.', ','))
#         elif isinstance(Mp, int):
#             mass_str.append(str(Mp))
#         else:
#             print('MP is', type(Mp))
#     dirs = ['hypatia' + ms + 'M_1600K_88Fe' for ms in mass_str]
#     colours = colorize(np.arange(len(dirs)), cmap=cmap, vmin=vmin, vmax=vmax)[0]
#
#     # get colour for 1ME
#     c_1M_alpha = colours[ii_1M]
#     c_1M_alpha[3] = 0.25  # set alpha
#
#     # set up scatter & histogram gridspec
#     fig = plt.figure(figsize=(8, 12))
#     dats = rw.read_dir(px.perplex_path_default + 'output/' + dirs[ii_1M] + '/')
#     fig, axes, ax_histx, axes_histy = plotpx.pop_scatterhist_subplotwrapper(dats, "mgsi", y_vars,
#                                                                             y_scales=[1e6 / (masses[ii_1M] * p.M_E),
#                                                                                       1e6],
#                                                                             xlabel='Mg/Si', ylabels=ylabels,
#                                                                             histx_kwargs={'color': '0.9',
#                                                                                           'edgecolor': 'k',
#                                                                                           'linewidth': 0.5, 'bins': 55},
#                                                                             histy_kwargs={'color': c_1M_alpha,
#                                                                                           'edgecolor': colours[ii_1M],
#                                                                                           'linewidth': 0.9, 'bins': 40},
#                                                                             earth=False, xlim=xlim, ylims=ylims, ms=ms,
#                                                                             save=False, show=False, lim_histx=(0, 6),
#                                                                             lims_histy=[(0, 0.01)] * 2,
#                                                                             annotate_n=False, labelsize=labelsize,
#                                                                             legsize=legsize,
#                                                                             c=colours[ii_1M], alpha=alpha, fig=fig)
#     # show Earth vline
#     for axx in [axes[0], axes[1], ax_histx]:
#         axx.axvline(earth.mgsi, c='k', alpha=0.2, lw=0.5, ls=(10, (15, 10)))
#
#     # overplot other masses
#     for ii, d in enumerate(dirs):
#         for iax, ax in enumerate(axes):
#             if ii != ii_1M:
#                 dats = rw.read_dir(px.perplex_path_default + 'output/' + d + '/')
#                 fig, ax = plotpx.pop_scatter(dats, "mgsi", y_vars[iax], y_scale=[1e6 / (masses[ii] * p.M_E), 1e6][iax],
#                                              xlabel='Mg/Si',
#                                              ylabel='', annotate_n=False,
#                                              earth=False, fig=fig, ax=ax, xlim=xlim, ms=ms,
#                                              ylim=ylims[iax], labelsize=labelsize, ticksize=ticksize, legsize=legsize,
#                                              save=False,
#                                              show=False, c=colours[ii], alpha=alpha)
#
#     # show Earth scaled Mg/Si
#     dats_e = rw.read_dir(px.perplex_path_default + 'output/MgSi_from_earth/')
#     for ii, ax in enumerate(axes):
#         x_e, y_e = [], []
#         for dat in dats_e:
#             x_e.append(dat.mgsi)
#             y_e.append(eval('dat.' + y_vars[ii]) / p.M_E * 1e6)
#         x_e, y_e = zip(*sorted(zip(x_e, y_e)))
#         ax.plot(x_e, y_e, c='k', ls=':', label='Compositionally-scaled Earth')
#         ax.scatter(earth.mgsi, eval('earth.' + y_vars[ii]) * [1e6 / (masses[ii_1M] * p.M_E), 1e6][ii], label='Earth',
#                    c='k',
#                    marker='$\oplus$', s=200, zorder=200)
#
#     # make legends
#     leg_bbox_to_anchor = (1.01, 1)
#     cax = colourised_legend(axes[0], colours, [str(m) + ' $M_\oplus$' for m in masses], title=None,
#                             # r'$\bf{Planet}$ $\bf{mass}$',
#                             legsize=legsize, titlesize=legsize, markersize=5, bbox_to_anchor=leg_bbox_to_anchor,
#                             handletextpad=0.1)
#     for ax_histy in axes_histy:
#         ax_histy.set_zorder(-1)  # to go below leg
#     # axes[0].set_xticks([])
#     axes[0].legend(handles=[mlines.Line2D([], [], color='k', marker='$\oplus$',
#                                           markersize=13, lw=0, label='Earth'),
#                             mlines.Line2D([], [], color='k', marker=None,
#                                           lw=1, ls=':', label='Compositionally-scaled Earth'),
#                             ], frameon=False, fontsize=legsize)
#     # fig.supylabel(ylabel, fontsize=labelsize, y=0.45)
#     fig.savefig(plotpx.fig_path + fname + '.png', bbox_inches='tight', dpi=300)


def scatterhist_mgsi_allmasses_kepler(masses=[0.1, 0.3, 0.5, 1, 2, 3, 4, 5], cmap='rainbow', earth=None, sun=None,
                                      alpha=0.4, ylim=None, fformat='.png', xlabelpad=9,
                                      labelsize=14, ticksize=12, legsize=12, xlim=(0, 3), vmin=None, vmax=None, thin=None,
                                      ms=200, ratios=[7, 1], output_path=px.perplex_path_default + 'output/apollo/', **kwargs):
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family': 'serif',
                      'serif': ['Computer Modern Roman']})  # this only works if outside function for some reason
    ii_1M = masses.index(1)
    nrows = 2

    y_var = "mass_h2o_total"
    ylabel = 'Water capacity (ppm wt)'
    fname = 'mgsi_c_h2o_exo_r2'
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
                          height_ratios=[ratios[1], ratios[0], 1.5],  # change last one for kepler panel
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    ax_histx = fig.add_subplot(gs[0, 0])
    axes = []
    axes.append(fig.add_subplot(gs[1, 0], sharex=ax_histx))  # model
    axes.append(fig.add_subplot(gs[2, 0]))  # exo
    ax_histy = fig.add_subplot(gs[1, 1], sharey=axes[0])

    dirs = ['hypatia' + ms + 'M_1600K_88Fe_hires' for ms in mass_str]

    colours = colorize(masses, cmap=cmap, vmin=vmin, vmax=vmax)[0]

    # get colour for 1ME
    c_1M_alpha = colours[ii_1M]
    c_1M_alpha[3] = 0.8  # 0.25  # set alpha

    # set up scatter & histogram gridspec
    dats = rw.read_dir(output_path + dirs[ii_1M] + '/')
    fig, axes[0], ax_histx, ax_histy = plotpx.pop_scatterhist(dats, "mgsi", y_var,
                                                              y_scale=1e6 / (masses[ii_1M] * p.M_E),
                                                              xlabel='', ylabel='',
                                                              histx_kwargs={'color': '0.9', 'edgecolor': 'k',
                                                                            'linewidth': 0.5, 'bins': 55},
                                                              histy_kwargs={'color': c_1M_alpha,
                                                                            'edgecolor': colours[ii_1M],
                                                                            'linewidth': 0.9, 'bins': 40},
                                                              earth=False, xlim=xlim, ylim=ylim, ms=ms,
                                                              save=False, show=False, lim_histx=(0, 6),
                                                              lim_histy=(0, 0.01),
                                                              annotate_n=False, labelsize=labelsize, legsize=legsize,
                                                              c=colours[ii_1M], alpha=alpha,
                                                              fig=fig, ax=axes[0], ax_histy=ax_histy, ax_histx=ax_histx,
                                                              show_histx=True, show_histy=True, thin=thin)
    # show Earth vline
    for axx in [axes[0], axes[1], ax_histx]:
        axx.axvline(earth.mgsi, c='k', alpha=0.2, lw=0.5, ls=(10, (15, 10)))
    for axy in [axes[0], ax_histy]:
        axy.axhline(eval('earth.' + y_var) * 1e6 / (masses[ii_1M] * p.M_E),
                    c='k', alpha=0.2, lw=0.5, ls=(10, (15, 10)))

    # overplot other masses
    for ii, d in enumerate(dirs):
        if ii != ii_1M:
            dats = rw.read_dir(output_path + d + '/', verbose=True)
            fig, axes[0] = plotpx.pop_scatter(dats, "mgsi", y_var, y_scale=1e6 / (masses[ii] * p.M_E), xlabel='',
                                              ylabel='', annotate_n=False,
                                              earth=False, fig=fig, ax=axes[0], xlim=xlim, ms=ms,
                                              ylim=ylim, labelsize=labelsize, ticksize=ticksize, legsize=legsize,
                                              save=False,
                                              show=False, c=colours[ii], alpha=alpha, thin=thin)

            # # show Earth-scaled Mg/Si
            # dats_e = rw.read_dir(px.perplex_path_default + 'output/MgSi_from_earth_M' + str(masses[ii]) + '/')
            # x_e, y_e = [], []
            # for dat in dats_e:
            #     x_e.append(dat.mgsi)
            #     y_e.append(eval('dat.' + y_var) / (p.M_E * masses[ii]) * 1e6)
            # x_e, y_e = zip(*sorted(zip(x_e, y_e)))
            # axes[0].plot(x_e, y_e, c=colours[ii], ls=':', label=None, alpha=0.5)

    # only 1 ME comp scaling for lower plot & show markers
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

    # show real planets
    exo_names, exo_masses, exo_mgsi = load_exoplanets('earthsize_planets_1600K_88Fe/', 'mgsi', scale=1)
    _, _, exo_w = load_exoplanets('earthsize_planets_1600K_88Fe/', y_var, scale=1)

    exo_colours = colorize(exo_masses, cmap=cmap, vmin=vmin, vmax=vmax)[0]

    axes[0].set_xticks([])
    ax = axes[1]
    ax.set_ylim([250, 700])
    ax.set_xlim(xlim)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.set_xlabel('Mg/Si', fontsize=labelsize, labelpad=xlabelpad)

    # errorbars
    exo_names_min, _, exo_ys_min = load_exoplanets('earthsize_planets_1600K_88Fe_mgsimin/', y_var, scale=1)
    exo_names_max, _, exo_ys_max = load_exoplanets('earthsize_planets_1600K_88Fe_mgsimax/', y_var, scale=1)
    _, _, exo_xs_min = load_exoplanets('earthsize_planets_1600K_88Fe_mgsimin/', 'mgsi', scale=1)
    _, _, exo_xs_max = load_exoplanets('earthsize_planets_1600K_88Fe_mgsimax/', 'mgsi', scale=1)

    medy, minw, maxw = [], [], []
    for ii, pl in enumerate(exo_names):
        # print('name', pl, '=', exo_names_min[ii], '=', exo_names_max[ii])
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
        # print(['xerr', [xerr[0, jj]], [xerr[1, jj]]])
        # print(['yerr', [yerr[0, jj]], [yerr[1, jj]]])
        ax.errorbar(xx, exo_w[jj] / (exo_masses[jj] * p.M_E) * 1e6,
                    xerr=[[xerr[0, jj]], [xerr[1, jj]]],
                    yerr=[[yerr[0, jj]], [yerr[1, jj]]],
                    elinewidth=1, ecolor=exo_colours[jj], alpha=0.6,
                    marker='.', ms=0, zorder=99,  # cmap=cmap,
                    # vmin=0.1, vmax=5
                    )
        ax.scatter(exo_mgsi[jj], exo_w[jj] / (exo_masses[jj] * p.M_E) * 1e6,
                   c=exo_colours[jj],  # alpha=0.4,
                   marker='.', s=100, zorder=100, #cmap=cmap, vmin=vmin, vmax=vmax
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

    # make legends
    leg_bbox_to_anchor = (1.01, 1)
    cax = colourised_legend(axes[0], colours, [str(m) + ' $M_\oplus$' for m in masses],  # alpha=alpha,
                            title=None,  # r'$\bf{Planet}$ $\bf{mass}$',
                            legsize=legsize, titlesize=legsize, markersize=5,
                            bbox_to_anchor=leg_bbox_to_anchor, handletextpad=0.1)
    ax_histy.set_zorder(-1)  # to go below leg
    axes[0].legend(handles=[mlines.Line2D([], [], color='k', marker='*',
                                          markersize=8, lw=0, label='Solar'),
                            mlines.Line2D([], [], color='xkcd:dark blue', marker='$\oplus$',
                                          markersize=10, lw=0, label='Earth'),
                            mlines.Line2D([], [], color='k', marker=None,
                                          lw=1, ls=':', label='Compositionally-scaled Earth'),
                            ], frameon=False, fontsize=legsize)

    fig.supylabel(ylabel, fontsize=labelsize, y=0.46, x=0.005)
    fig.savefig(plotpx.fig_path + fname + fformat, bbox_inches='tight', dpi=200)


masses = [0.1, 0.5, 1, 2, 3]

# scatterhist_mgsi_allmasses(masses=masses, cmap='YlOrBr_r', earth=earth, alpha=0.4, vmax=len(masses), xlim=(0.25, 2), ms=100,
#                            fname='test', end='_1600K_88Fe_hires', output_path=px.perplex_path_default + 'output/apollo/')

# scatterhist_mgsi_allmasses_multipanel(masses=masses, cmap='YlOrBr_r', earth=earth, alpha=0.4, vmax=len(masses), xlim=(0.25, 2), ms=100)

scatterhist_mgsi_allmasses_kepler(masses=masses, cmap='YlOrBr_r', earth=earth, sun=sun, alpha=0.4, vmin=0.1, vmax=4,
                                  xlim=(0.25, 2), ms=30, thin=3,
                                  output_path=px.perplex_path_default + 'output/apollo/', fformat='.pdf')

# some of these 0.1 M_E planets with low Mg/Si have no ol/wad, some ring at bottom but lots of qtz/coes/stv
# so have a kind of gradual incline of increased c_h2o from 10-20 GPa rather than sharp mtz at 1 GPa, same as 1ME planets in figure 1 just cut at 20 GPa.
# probably won't have pool of hydrous melting because no sharp boundary. no clue about melt buoyancy of si-rich melt at like 15 GPa
# still a clear shift at breakdown of opx - so use this to define "OBM"
# even with this there's some super weird cases e.g. 2MASS 19394601-2544539 has ol breaks down to fp and gt, no wad, short opx layer

plt.show()
