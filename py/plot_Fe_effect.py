import numpy as np
import parameters as p
import ask_hypatia as hyp
import perplexdata as px
import plot_perplex as plotpx
import saturation as sat
import main as rw
import parameters as p
import matplotlib.pyplot as plt
from useful_and_bespoke import cornertext, colorize, colourbar, colourised_legend
from UM_mass_scaling import um_mass_scaling
import saturation as sat
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
print(earth.core_eff_from_cmf())

""" plot Fe dependence """
def plot_XFe_dependence(M_p=1, Tp=1600, x_Fe=None, sigma=1, labelsize=16, figsize=(4, 4), xlabel='Mantle Fe fraction',
                        legsize=12, ticksize=12, xpad=20, ypad=20, title='', earth=None, colours=['xkcd:gold', 'k'],
                        xlims=None, save=True, fig=None, ax=None, **kwargs):
    """ make plot for several masses showing mantle structure and relative water in each zone"""

    plt.rc('text', usetex=True)
    plt.rc('font', **{'family': 'serif',
                      'serif': ['Computer Modern Roman']})  # this only works if outside function for some reason


    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    if x_Fe is None:
        x_Fe = [0.5, 0.6, 0.7, 0.8, 0.88, 0.9, 0.9999999]
    x_Fe_mtl = [x * -1 + 1 for x in x_Fe]
    dirs = [px.perplex_path_default + 'output/hypatia' + str(M_p) + 'M_' + str(Tp) + 'K_' + str(int(x * 100)) + 'Fe/'
            for x in x_Fe]

    datalabels = ['Upper mantle', 'Whole mantle']
    y_names = ['mass_h2o_um', 'mass_h2o_total']
    lws = [2, 2]
    earthlabel = 'Earth'
    for ii in range(len(y_names)):
        if earth:
            earth.core_eff = earth.core_eff_from_cmf()
            print('earth core eff', earth.core_eff)
            ax.scatter(earth.core_eff * -1 + 1, eval('earth.' + y_names[ii]) * p.TO ** -1, label=earthlabel, c=colours[ii],
                       marker='$\oplus$', s=100, zorder=100)
            earthlabel = None
        fig, ax = plotpx.compare_pop_fillbetween(dirs, 'core_eff*-1 + 1', y_names[ii], x_scale=1, y_scale=p.TO ** -1,
                                                 xlabel=xlabel, ylabel='Water capacity (OM)', title=title,
                                                 save=False, show=False, labelsize=labelsize, show_med=True,
                                                 legsize=legsize, ticksize=ticksize, earth=None, head=-1,
                                                 patch_kwargs={'color': colours[ii], 'alpha': 0.3}, c=colours[ii],
                                                 line_kwargs={'color': colours[ii], 'lw': lws[ii]}, fig=fig, ax=ax,
                                                 figsize=(6, 4),
                                                 sigma=sigma, show_n=False, datalabel=datalabels[ii], **kwargs)
    ax.set_xlim(0, np.max(x_Fe_mtl))
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    h1, l1 = ax.get_legend_handles_labels()
    ax.legend(handles=[h1[1], h1[2], h1[0]], labels=[l1[1], l1[2], l1[0]], frameon=False, fontsize=legsize, loc='lower right',
              bbox_to_anchor=(0.9, 1), ncol=2, mode='expand'
              )
    plt.tight_layout()
    # fig, *ax = dark_background(fig, ax, )
    if save:
        fig.savefig(plotpx.fig_path + 'x_Fe_dependence.png', bbox_inches='tight',
                    facecolor=fig.get_facecolor())


plot_XFe_dependence(M_p=1, Tp=1600, sigma=1, x_Fe = [0.65, 0.7, 0.75, 0.8, 0.85, 0.88, 0.90, 0.93, 0.95, 0.97, 0.9999999],
                    title='', #r'1 $M_\oplus$, 1600 K',
                    xlabel='Mantle/bulk Fe',  # 'Mantle Fe fraction'
                    save=True, earth=earth, colours=['xkcd:seafoam', 'xkcd:blood orange'])

plt.show()
