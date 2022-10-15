import numpy as np
import parameters as p
import perplexdata as px
import plot_perplex as plotpx
import main as rw
import parameters as p
import matplotlib.pyplot as plt
import matplotlib.colors
from useful_and_bespoke import cornertext, colorize, colourbar, colourised_legend
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes

""" plot water profiles and compositoin gridspec for median and extrema """


def demo_composition_gridspec(stars=None, mgsi=None, dats=None, colours_comp=None, ls_comp=None, labelsize=14,
                              legsize=10, M_p=1, Tp=1600, prof_lw=1, labelpad=6, ticklabelsize=8,
                              p_max=35, phases_order='default', log_profile=False, orientation='horizontal',
                              comp_var='mgsi', cmap='tab20', figtitle='phase_demo', save=True, fformat='.png'):
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
    if phases_order == 'default':
        phases_order = ['gt', 'cpx', 'opx', 'hpcpx', 'ol', 'wad', 'ring', 'pv', 'qtz', 'coes',
                        'st', 'wus', 'dvm',  # 'ppv', 'fapv'
                        ]  # may want an (plag) at top
    if dats is not None:
        pass
    elif stars is not None:
        dats = [rw.read_name(output_path=px.output_parent_default + 'hypatia1M_1600K_88Fe/', star=st, M_p=p.M_E,
                             core_efficiency=0.8831461545794602, Tp=1600) for st in stars]
    elif mgsi is not None:  # update Earth value
        dats = []
        for z in mgsi:
            oxides = rw.update_MgSi(MgSi=z, oxides=px.wt_oxides_MD95)
            dats.append(rw.build_planet(test_oxides=oxides, test_CMF=0.325, M_p=M_p * p.M_E, Tp=Tp, get_saturation=True,
                                        output_parent_path=px.output_parent_default + '/MgSi_from_earth_fig1/', n=1200,
                                        plot_all=False))
    if colours_comp is None:  # default is 3 mgsi ratios
        cmap_comp = matplotlib.cm.get_cmap('bone_r')
        colours_comp = [cmap_comp(0.3), cmap_comp(0.5), cmap_comp(0.999)]
        # colours_mgsi = ['xkcd:silver', 'xkcd:olive', 'xkcd:black']
    if ls_comp is None:
        ls_comp = ['--', '-', '-.']

    if cmap is None:
        # create custom cmap - 13 phases
        c_phases = [plt.cm.tab20b(14),  # gt
                    plt.cm.tab20b(10),  # cpx
                    plt.cm.tab20b(11),  # opx
                    plt.cm.tab20b(9),  # hpcpx
                    plt.cm.tab20b(7),  # ol
                    plt.cm.tab20b(5),  # wad
                    plt.cm.tab20b(4),  # ring
                    plt.cm.tab20c(19),  # pv
                    plt.cm.tab20b(19),  # qtz
                    plt.cm.tab20b(18),  # coes
                    plt.cm.tab20b(17),  # st
                    plt.cm.tab20b(3),  # wus
                    plt.cm.tab20b(0),  # capv
                    ]
        cmap = matplotlib.colors.ListedColormap(c_phases)

    fig = plt.figure()
    if orientation == 'horizontal':
        gs = fig.add_gridspec(2, 1,
                              # width_ratios=(2, 1, 1, 1),
                              height_ratios=(2, 3),
                              # left=0.1, right=0.9, bottom=0.1, top=0.9,
                              wspace=0.05, hspace=0.15)
        gs0 = gridspec.GridSpecFromSubplotSpec(len(dats), 1, height_ratios=[1] * len(dats), hspace=0.07,
                                               subplot_spec=gs[1])
    elif orientation == 'vertical':
        gs = fig.add_gridspec(1, 2,
                              width_ratios=(7, 1.5),
                              hspace=0.05, wspace=0.1)
        gs0 = gridspec.GridSpecFromSubplotSpec(1, len(dats), width_ratios=[1] * len(dats), wspace=0.05,
                                               subplot_spec=gs[0])

    ax_profs = fig.add_subplot(gs[1])
    axes_comp = [fig.add_subplot(gs0[ii]) for ii in range(len(dats))]

    comp_ax_legend = True
    # compositions
    for ii, (dat, ax) in enumerate(zip(dats, axes_comp)):
        if ii == 1:
            comp_ax_ylabel = r'Mineral mode (wt\%)'
        else:
            comp_ax_ylabel = ''
        if ii == 0:
            xlabel = 'Pressure (GPa)'
        else:
            xlabel = ''
        if comp_var == 'mgsi':
            comp_leg_title = r'$\bf{Mg/Si}$'
            prof_leg_label = '{:.1f}'.format(dat.mgsi)
            comp_ax_title = r'$\bf{Mg/Si} =$ ' + prof_leg_label
        elif comp_var == 'mgfe':
            comp_leg_title = r'$\bf{Mg\#}$'
            prof_leg_label = '{:.1f}'.format(dat.wt_oxides['MgO'] / (dat.wt_oxides['FeO'] + dat.wt_oxides['MgO']))
            comp_ax_title = r'Mg\# = ' + prof_leg_label
        elif comp_var == 'T':
            comp_leg_title = r'$\bf{T_p}$'
            prof_leg_label = '{:.1f}'.format(dat.temperature[-1])
            comp_ax_title = r'$T_p$ = ' + prof_leg_label

        if orientation == 'horizontal':
            leg_bbox_to_anchor = (1.02, 2.1)
            comptitle = ''
        elif orientation == 'vertical':
            leg_bbox_to_anchor = (1.02, 0.5)
            comptitle = ''  # comp_ax_title
        fig, ax = plotpx.single_composition(dat, which='pressure', modality_type='phase', comp_stacked=True, save=False,
                                            ylabel=comp_ax_ylabel, xlabel=xlabel, cmap=cmap, override_ax_arrow=True,
                                            labelsize=labelsize, legsize=legsize,
                                            plot_phases_order=phases_order, p_max=p_max, orientation=orientation,
                                            fig=fig, ax=ax, title=comptitle, make_legend=False, ylabelpad=labelpad,
                                            xlabelpad=labelpad,
                                            leg_bbox_to_anchor=leg_bbox_to_anchor, legtitle=r'$\bf{Phase}$')
        # ax.set_yticks([])
        if orientation == 'horizontal':
            # only label Mg/Si in bbox for horizontal, otherwise it's title
            ax.text(0.98, 0.4, comp_ax_title, transform=ax.transAxes, fontsize=legsize,  # weight='bold',
                    va='center', ha='right', c=['k', 'k', 'w'][ii],
                    bbox=dict(boxstyle='round', fc=colours_comp[ii], ec='k', alpha=0.4, ls=ls_comp[ii]))
        elif orientation == 'vertical':
            ax.text(0.5, 0.97, comp_ax_title, transform=ax.transAxes, fontsize=legsize,
                    va='top', ha='center', c='k',
                    # bbox=dict(boxstyle='round', fc=(0, 0, 0, 0), ec=colours_comp[ii], ls=ls_comp[ii], lw=0.5)
                    )
            for spine in ['bottom', 'left', 'top', 'right']:
                ax.spines[spine].set_linestyle(ls_comp[ii])
                ax.spines[spine].set_color(colours_comp[ii])
                ax.spines[spine].set_linewidth(prof_lw)
            ax.tick_params(axis='both', colors=colours_comp[ii])
            plt.setp(ax.get_xticklabels(), color='k')
            plt.setp(ax.get_yticklabels(), color='k')

        if orientation == 'horizontal':
            profile_kwargs = {'reverse_y': False, 'xmin': 1000e-4, 'xmax': p_max, 'ymin': 0, 'ymax': 1,
                              'label_x': False,
                              'label_y': True, 'leg_bbox_to_anchor': (1.02, 1.1), 'ax_ticks': None}
        elif orientation == 'vertical':
            profile_kwargs = {'reverse_y': True, 'ymin': 1000e-4, 'ymax': p_max, 'xmin': 0, 'xmax': 1,
                              'label_x': True,
                              'label_y': False,  'leg_bbox_to_anchor': (0.5, 0.999),
                              'xlabelpad': labelpad,
                              'legloc': 'upper center'
                              }
        fig, ax_profs = plotpx.profile(dat, 'c_h2o', independent_ax='pressure',
                                       ax_label=r'Saturation (wt\% H$_2$O)', scale=1e2,
                                       c=colours_comp[ii], ls=ls_comp[ii], lw=prof_lw, alpha=0.8,
                                       labelsize=labelsize, legsize=legsize,
                                       fig=fig, ax=ax_profs, log=log_profile,
                                       orientation=orientation, save=False,
                                       leg_label=prof_leg_label, #leg_bbox_to_anchor=leg_bbox_to_anchor,
                                       legtitle=comp_leg_title,
                                       leg_kwargs={'handlelength': 1},
                                       **profile_kwargs)

        ax.tick_params(axis='both', which='major', labelsize=ticklabelsize)
        ax_profs.tick_params(axis='both', which='major', labelsize=ticklabelsize)

        if orientation == 'horizontal':
            if ii < len(dats) - 1:
                ax.tick_params(axis="x", labelbottom=False)
                ax.set_xticks([])
            ax_profs.tick_params(axis="x", labelbottom=False)
        if orientation == 'vertical':
            if ii > 0:
                ax.tick_params(axis="y", labelleft=False)
                ax.set_yticks([])
                ax.set_xticks([50, 100])
            ax_profs.yaxis.tick_right()
            ax_profs.yaxis.set_ticklabels([])
            ax_profs.set_xticks([0.2, 0.8]) # [0.1, 0.5, 0.9]
            # ax_profs.tick_params(axis="y", labelleft=False)
        comp_ax_legend = False  # only make legend once

        print('\n', dat.name, 'um water mass', dat.mass_h2o_um / p.TO, 'total water mass', dat.mass_h2o_total / p.TO,
              'TO')

    if save:
        plt.tight_layout()
        fig.savefig(plotpx.fig_path + figtitle + fformat, bbox_inches='tight', dpi=300)


demo_composition_gridspec(mgsi=[0.7244359600749921, 1.0715193052376069, 1.4125375446227555],
                          comp_var='mgsi', labelsize=10, legsize=10, cmap=None, fformat='.pdf',
                          figtitle='phase_demo_ver_blank', orientation='vertical')

# plt.show()
