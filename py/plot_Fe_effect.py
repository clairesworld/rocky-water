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
                        maxIter=30, tol=1e-4,  n=1200,
                        Tp=Tp,  # core_efficiency=0.8,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Earth_' + str(Tp) + 'K',
                        )
sun = rw.build_planet(M_p=1 * p.M_E, star='sun',
                        maxIter=30, tol=1e-4, n=1200,
                        Tp=Tp,  core_efficiency=0.8838,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Sun_' + str(Tp) + 'K',
                        )
print(earth.core_eff_from_cmf())
xe, ye = earth.core_eff*-1 + 1, earth.mass_h2o_total / p.TO
print('earth out of fn, x, y', xe, ye)

""" plot Fe dependence """

def plot_XFe_dependence(M_p=1, Tp=1600, x_Fe=None, sigma=1, labelsize=16, figsize=(4, 4), xlabel='Mantle Fe fraction',
                        legsize=12, ticksize=12, xpad=20, ypad=20, title='', colours=['xkcd:gold', 'k'], fformat='.png',
                        xlims=None, save=True, fig=None, ax=None, earth=None, sun=None, **kwargs):
    """ make plot for several masses showing mantle structure and relative water in each zone"""

    plt.rc('text', usetex=True)
    plt.rc('font', **{'family': 'serif',
                      'serif': ['Computer Modern Roman']})  # this only works if outside function for some reason

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    if x_Fe is None:
        x_Fe = [0.5, 0.6, 0.7, 0.8, 0.88, 0.9, 0.9999999]
    x_Fe_mtl = [x * -1 + 1 for x in x_Fe]
    dirs = [px.perplex_path_default + 'output/apollo/hypatia' + str(M_p) + 'M_' + str(Tp) + 'K_' + str(int(x * 100)) + 'Fe_hires/'
            for x in x_Fe]

    datalabels = ['Upper mantle', 'Whole mantle']
    y_names = ['mass_h2o_um', 'mass_h2o_total']
    lws = [1, 1]
    earthlabel, sunlabel = 'Earth', 'Solar'
    for ii in range(len(y_names)):
        if earth:
            earth.core_eff = earth.core_eff_from_cmf()
            print('earth core eff', earth.core_eff)
            ax.scatter(sun.core_eff * -1 + 1, eval('sun.' + y_names[ii]) * p.TO ** -1, label=sunlabel,
                       c=colours[ii],
                       marker='*', s=80, zorder=99)
            ax.scatter(earth.core_eff * -1 + 1, eval('earth.' + y_names[ii]) * p.TO ** -1, label=earthlabel,
                       c=colours[ii],
                       marker='$\oplus$', s=150, zorder=100)
            earthlabel, sunlabel = None, None


        fig, ax = plotpx.compare_pop_fillbetween(dirs, 'core_eff*-1 + 1', y_names[ii], x_scale=1, y_scale=p.TO ** -1,
                                                 xlabel=xlabel, ylabel='Water capacity (OM)', title=title,
                                                 save=False, show=False, labelsize=labelsize, show_med=True,
                                                 legsize=legsize, ticksize=ticksize, earth=None, sun=None, head=-1,
                                                 patch_kwargs={'color': colours[ii], 'alpha': 0.3}, c=colours[ii],
                                                 line_kwargs={'color': colours[ii], 'lw': lws[ii]}, fig=fig, ax=ax,
                                                 figsize=(6, 4),
                                                 sigma=sigma, show_n=False, datalabel=None, legend=False, **kwargs)
    ax.set_xlim(0, np.max(x_Fe_mtl))
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    # ax.scatter(x_Fe_mtl, [1]*len(x_Fe_mtl), c='k')

    # label UM and whole mantle and create legend
    dx, dy = 0.35 - 0.25, 3.1819100995387464 - 2.8729511088776563  # run plotting to print these values out
    angle = np.degrees(np.arctan2(dy, dx))
    ax.text(0.34, 2.73, 'Whole mantle', fontsize=legsize, ha='right', rotation=angle, rotation_mode='anchor',transform_rotates_text=True,
            c=colours[1])

    dx, dy = 0.35 - 0.25, 2.116008344883379 - 1.7853485383846979  # run plotting to print these values out
    angle = np.degrees(np.arctan2(dy, dx))
    ax.text(0.34, 1.67, 'Upper mantle', fontsize=legsize, ha='right', rotation=angle, rotation_mode='anchor',transform_rotates_text=True,
            c=colours[0])

    # label earth and sun
    xe, ye = earth.core_eff * -1 + 1, earth.mass_h2o_total / p.TO
    xs, ys = sun.core_eff * -1 + 1, sun.mass_h2o_total / p.TO
    print('earth, x, y', xe, ye)
    ax.text(xe - 0.01, ye, 'Earth', fontsize=legsize, c=colours[1], ha='right', va='bottom')  # whole mantle
    ax.text(xs + 0.01, ys - 0.02, 'Solar', fontsize=legsize, c=colours[1], ha='left', va='top')  # whole mantle

    # h1, l1 = ax.get_legend_handles_labels()
    # ax.legend(#handles=[h1[0], h1[1]], labels=[l1[0], l1[1]],
    #     handles=[mlines.Line2D([], [], color='k', marker='*',
    #                            markersize=8, lw=0, label='Solar'),
    #              mlines.Line2D([], [], color='k', marker='$\oplus$',
    #                            markersize=10, lw=0, label='Earth'),
    #              ],
    #           frameon=False, fontsize=legsize, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #    ncol=2, mode="expand", borderaxespad=0.
    #           )

    # # test LM isolated
    # for folder in dirs:
    #     dats = rw.update_dir(folder, px.PerplexData.get_lm_water, store=True)
    # fig, ax = plotpx.compare_pop_fillbetween(dirs, 'core_eff*-1 + 1', 'mass_h2o_lm', x_scale=1, y_scale=p.TO ** -1,
    #                                          xlabel=xlabel, ylabel='Water capacity (OM)', title=title,
    #                                          save=False, show=False, labelsize=labelsize, show_med=True,
    #                                          legsize=legsize, ticksize=ticksize, earth=None, sun=None, head=-1,
    #                                          patch_kwargs={'color': 'k', 'alpha': 0.01, 'hatch': '//'}, c=colours[ii],
    #                                          line_kwargs={'color': 'k', 'lw': 1, 'ls':'--'}, fig=fig, ax=ax,
    #                                          figsize=(6, 4),
    #                                          sigma=sigma, show_n=False, datalabel=None, legend=False, **kwargs)

    plt.tight_layout()
    # fig, *ax = dark_background(fig, ax, )
    if save:
<<<<<<< HEAD
        fig.savefig(plotpx.fig_path + 'x_Fe_dependence' + fformat, bbox_inches='tight', dpi=400,
=======
        fig.savefig(plotpx.fig_path + 'x_Fe_dependence.png', bbox_inches='tight', dpi=400,
>>>>>>> 780fc79a30ed92265a34ee2e0dc76d862124a3a0
                    facecolor=fig.get_facecolor())


plot_XFe_dependence(M_p=1, Tp=1600, sigma=1, labelsize=14, fformat='.pdf',
                    x_Fe=[0.65, 0.7, 0.75, 0.8, 0.85, 0.88, 0.90, 0.93, 0.95, 0.97, 0.98, 0.9999999],
                    title='',  # r'1 $M_\oplus$, 1600 K',
                    xlabel='Mantle/bulk Fe (mol/mol)',  # 'Mantle Fe fraction'
                    save=True, earth=earth, sun=sun, colours=['xkcd:sea', 'xkcd:blood orange'])  # ocean

plt.show()
