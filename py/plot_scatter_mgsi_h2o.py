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


""" scatter hist of mgsi effect on upper mantle water but over all masses """
def scatterhist_mgsi_allmasses(masses=[0.1, 0.3, 0.5, 1, 2, 3, 4, 5], cmap='rainbow', earth=None, alpha=0.4, ylim=None,
                               labelsize=14, ticksize=12, legsize=12, vmin=None, vmax=None, show_histy=True,
                               which='total'):
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family': 'serif',
                      'serif': ['Computer Modern Roman']})  # this only works if outside function for some reason
    ii_1M = masses.index(1)

    if which == 'total':
        y_var = "mass_h2o_total"
        ylabel = 'Water storage capacity (ppm)'
        fname = 'mgsi_c_h2o_mantle_scatter_all'
        if ylim is None:
            ylim = (100, 2000)  # (0, 0.3),
    elif which == 'um':
        y_var = "mass_h2o_um"
        ylabel = 'UM water storage capacity (ppm)'
        fname = 'mgsi_c_h2o_um_scatter_all'
        if ylim is None:
            ylim = (0, 2000)
    mass_str = []
    for Mp in masses:
        if isinstance(Mp, float):
            mass_str.append(str(Mp).replace('.', ','))
        elif isinstance(Mp, int):
            mass_str.append(str(Mp))
        else:
            print('MP is', type(Mp))
    dirs = ['hypatia' + ms + 'M_1600K_88Fe' for ms in mass_str]
    colours = colorize(np.arange(len(dirs)), cmap=cmap, vmin=vmin, vmax=vmax)[0]

    # get colour for 1ME
    print('colour', colours[ii_1M])
    c_1M_alpha = colours[ii_1M]
    c_1M_alpha[3] = 0.25  # set alpha
    print('updated colour', colours[ii_1M])

    # set up scatter & histogram gridspec
    dats = rw.read_dir(px.perplex_path_default + 'output/' + dirs[ii_1M] + '/')
    fig, ax, ax_histx, ax_histy = plotpx.pop_scatterhist(dats, "mgsi", y_var, y_scale=1e6 / (masses[ii_1M] * p.M_E),
                                                         xlabel='Mg/Si', ylabel='', show_histx=True, show_histy=show_histy,
                                                         histx_kwargs={'color': '0.9', 'edgecolor': 'k',
                                                                       'linewidth': 0.5, 'bins': 50},
                                                         histy_kwargs={'color': c_1M_alpha, 'edgecolor': colours[ii_1M],
                                                                       'linewidth': 0.9, 'bins': 230},
                                                         earth=False, xlim=(0, 3), ylim=ylim,
                                                         bins=50, save=False, show=False, lim_histx=(0, 250), lim_histy=(0, 450),
                                                         annotate_n=False, labelsize=labelsize, legsize=legsize,
                                                         c=colours[ii_1M], alpha=alpha)
    # show Earth vline
    for axx in (ax, ax_histx):
        axx.axvline(earth.mgsi, c='k', alpha=0.2, lw=0.5, ls=(10, (15, 10)))

    # overplot other masses
    for ii, d in enumerate(dirs):
        if ii != ii_1M:
            dats = rw.read_dir(px.perplex_path_default + 'output/' + d + '/')
            fig, ax = plotpx.pop_scatter(dats, "mgsi", y_var, y_scale=1e6/(masses[ii]*p.M_E), xlabel='Mg/Si',
                                         ylabel=ylabel, annotate_n=False,
                                         earth=False, fig=fig, ax=ax, xlim=(0, 3),
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
    ax.scatter(earth.mgsi, eval('earth.' + y_var) / p.M_E * 1e6, label='Earth', c='k',
               marker='$\oplus$', s=200, zorder=200)

    # make legends
    if show_histy:
        leg_bbox_to_anchor = (1.01, 1)
    else:
        leg_bbox_to_anchor = (1.01, 1)
    cax = colourised_legend(ax, colours, [str(m) + ' $M_\oplus$' for m in masses], title=None,  # r'$\bf{Planet}$ $\bf{mass}$',
                            legsize=legsize, titlesize=legsize, markersize=5, bbox_to_anchor=leg_bbox_to_anchor, handletextpad=0.1)
    ax_histy.set_zorder(-1)  # to go below leg
    ax.legend(handles=[mlines.Line2D([], [], color='k', marker='$\oplus$',
                                     markersize=13, lw=0, label='Earth'),
                       mlines.Line2D([], [], color='k', marker=None,
                                     lw=1, ls=':', label='Compositionally-scaled Earth'),
                       ], frameon=False, fontsize=legsize)
    fig.savefig(plotpx.fig_path + fname + '.png', bbox_inches='tight', dpi=300)


masses = [0.1, 0.5, 1, 2, 3]
scatterhist_mgsi_allmasses(masses=masses, cmap='YlOrBr_r', earth=earth, alpha=0.4, vmax=len(masses))

plt.show()
