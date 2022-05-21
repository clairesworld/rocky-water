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
                               hist_kwargs={'color': '0.5', 'edgecolor': 'k'}, labelsize=14, vmin=None, vmax=None):
    mass_str = []
    for Mp in masses:
        if isinstance(Mp, float):
            mass_str.append(str(Mp).replace('.', ','))
        elif isinstance(Mp, int):
            mass_str.append(str(Mp))
        else:
            print('MP is', type(Mp))
    dirs = ['hypatia' + ms + 'M' for ms in mass_str]
    colours = colorize(np.arange(len(dirs)), cmap=cmap, vmin=vmin, vmax=vmax)[0]

    # set up scatter & histogram gridspec
    dats = rw.read_dir(px.perplex_path_default + 'output/' + dirs[0] + '/')
    fig, ax, _, _ = plotpx.pop_scatterhist(dats, "mgsi", "mass_h2o_total", y_scale=1e6/(masses[0]*p.M_E), xlabel='Mg/Si',
                                           ylabel='', show_histx=True, bins=50,
                                           hist_kwargs=hist_kwargs, earth=False, xlim=(0, 3),
                                           ylim=ylim, lim_histx=(0, 300), labelsize=labelsize, save=False, show=False, c=colours[0],
                                           alpha=alpha)

    # overplot other masses
    for ii, d in enumerate(dirs):
        if ii > 0:
            dats = rw.read_dir(px.perplex_path_default + 'output/' + d + '/')
            fig, ax = plotpx.pop_scatter(dats, "mgsi", "mass_h2o_total", y_scale=1e6/(masses[ii]*p.M_E), xlabel='Mg/Si',
                                         ylabel='Water storage capacity (ppm)', annotate_n=False,
                                         earth=False, fig=fig, ax=ax, xlim=(0, 3),
                                         ylim=ylim, labelsize=labelsize, save=False, show=False, c=colours[ii], alpha=alpha)

    # show Earth scaled Mg/Si
    dats_e = rw.read_dir(px.perplex_path_default + 'output/MgSi_from_earth/')
    x_e, y_e = [], []
    for dat in dats_e:
        x_e.append(dat.mgsi)
        y_e.append(dat.mass_h2o_total / p.M_E * 1e6)
    x_e, y_e = zip(*sorted(zip(x_e, y_e)))
    ax.plot(x_e, y_e, c='k', ls=':', label='Compositionally-scaled Earth')
    ax.scatter(earth.mgsi, earth.mass_h2o_total / p.M_E * 1e6, label='Earth', c='k',
               marker='$\oplus$', s=200, zorder=200)

    # make legends
    cax = colourised_legend(ax, colours, [str(m) + ' $M_\oplus$' for m in masses], title=r'$\bf{Planet}$ $\bf{mass}$', legsize=12,
                            titlesize=12, markersize=5)
    ax.legend(handles=[mlines.Line2D([], [], color='k', marker='$\oplus$',
                                     markersize=13, lw=0, label='Earth'),
                       mlines.Line2D([], [], color='k', marker=None,
                                     lw=1, ls=':', label='Compositionally-scaled Earth'),
                       ], frameon=False, fontsize=12)
    fig.savefig(plotpx.fig_path + 'mgsi_c_h2o_mantle_scatter_all.png', bbox_inches='tight')


scatterhist_mgsi_allmasses(masses=[0.1, 0.3, 0.5, 1, 2, 3, 4, 5], earth=earth, cmap='YlOrBr_r',
                           alpha=0.4, vmax=9, ylim=(0, 3000),  # (0, 0.3),
                           hist_kwargs={'color': '0.9', 'edgecolor': 'k', 'linewidth': 0.5 })

plt.show()
