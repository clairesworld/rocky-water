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

""" plot water profiles and compositoin gridspec for median and extrema """

def demo_composition_gridspec(stars=None, mgsi=None, colours=None, labelsize=12, legsize=10, M_p=1, Tp=1600, save=True):
    if stars is not None:
        dats = [rw.read_name(output_path=px.output_parent_default + 'hypatia1M/', star=st, M_p=p.M_E,
                             core_efficiency=0.8, Tp=1600) for st in stars]
    elif mgsi is not None:  # update Earth value
        dats = []
        for z in mgsi:
            oxides = rw.update_MgSi(MgSi=z, oxides=px.wt_oxides_Earth)
            dats.append(rw.build_planet(test_oxides=oxides, test_CMF=0.325, M_p=M_p * p.M_E, Tp=Tp, get_saturation=True,
                                        plot_all=False))
    if colours is None:
        colours = ['xkcd:silver', 'xkcd:olive', 'xkcd:black']
    bbox_colours = ['k', 'w', 'w']

    fig = plt.figure()
    gs = fig.add_gridspec(2, 1,
                          # width_ratios=(2, 1, 1, 1),
                          height_ratios=(2, 3),
                          # left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.15)

    gs0 = gridspec.GridSpecFromSubplotSpec(3, 1, height_ratios=(1, 1, 1), hspace=0.07, subplot_spec=gs[1])
    ax_profs = fig.add_subplot(gs[0])
    ax_comp0 = fig.add_subplot(gs0[0])
    ax_comp1 = fig.add_subplot(gs0[1])
    ax_comp2 = fig.add_subplot(gs0[2])

    comp_ax_legend = True
    # compositions
    for ii, (dat, ax) in enumerate(zip(dats, [ax_comp0, ax_comp1, ax_comp2])):
        if ii == 1:
            comp_label = 'Phase modality'#\n(wt.%)'
        else:
            comp_label = ''
        if ii == 2:
            xlabel = 'Pressure (GPa)'
        else:
            xlabel = ''
        fig, ax = plotpx.single_composition(dat, which='pressure', modality_type='phase', comp_stacked=True, save=False,
                                            ylabel=comp_label, xlabel=xlabel, cmap='tab20', override_ax_arrow=True,
                                            labelsize=labelsize, legsize=legsize,
                                            plot_phases_order=['gt', 'cpx', 'ol', 'opx', 'qtz', 'coes',
                                                               'hpcpx', 'wad', 'capv', 'ring', 'pv',
                                                               'st', 'wus',  #'ppv', 'fapv'
                                                               ], p_max=30,
                                            fig=fig, ax=ax, title='', make_legend=comp_ax_legend, ylabelpad=29.4,
                                            leg_bbox_to_anchor=(1.02, 2.1), legtitle=r'$\bf{Phase}$')
        ax.set_yticks([])
        ax.text(0.97, 0.85, 'Mg/Si = {:.1f}'.format(dat.mgsi), transform=ax.transAxes, fontsize=legsize, #weight='bold',
                va='top', ha='right', c=['k', 'k', 'w'][ii], bbox=dict(boxstyle='round', fc=colours[ii], ec='k', alpha=0.5))
        if ii < 2:
            # ax.tick_params(axis="x", labelbottom=False)
            ax.set_xticks([])

        fig, ax_profs = plotpx.profile(dat, 'c_h2o', independent_ax='pressure', reverse_y=False,
                                       ax_label='Saturation (ppm H$_2$O)', scale=1e6,
                                       c=colours[ii], lw=2, alpha=0.8,
                                       xmin=1000e-4, xmax=30, ymin=10, ymax=2e4, ax_ticks=None, label_x=False,
                                       label_y=True, labelsize=labelsize, legsize=legsize, legtitle=r'$\bf{Mg/Si}$',
                                       fig=fig, ax=ax_profs, log=True, leg_label='{:.1f}'.format(dat.mgsi),
                                       orientation='horizontal', leg_bbox_to_anchor=(1.02, 1.1), save=False)
        ax_profs.tick_params(axis="x", labelbottom=False)
        comp_ax_legend = False

    if save:
        plt.tight_layout()
        fig.savefig(plotpx.fig_path + 'phase_demo' + '.png', bbox_inches='tight')


demo_composition_gridspec(  # stars=['HIP 29295', '2MASS 19375133+4945541', '2MASS 19472449+4905040']
    mgsi=[1.5, 1, 0.5], labelsize=11, legsize=10)