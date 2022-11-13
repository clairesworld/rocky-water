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
                              legsize=10, M_p=1, Tp=1600, cumulative=False, profile_var='c_h2o', ylim=(0, 1),
                              p_max=30, phases_order='default', log_profile=False, orientation='horizontal',
                              comp_var='mgsi', cmap='tab20', figtitle='phase_demo', save=True, fformat='.png', comp_labels=None,
                              prof_ax_label=r'Saturation (wt.\% H$_2$O)', profile_scale=1e2, show_cum_mass=False):
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
        if show_cum_mass:
            gs = fig.add_gridspec(3, 1,
                                  # width_ratios=(2, 1, 1, 1),
                                  height_ratios=(2, 2, 3),
                                  # left=0.1, right=0.9, bottom=0.1, top=0.9,
                                  wspace=0.05, hspace=0.15)
        else:
            gs = fig.add_gridspec(2, 1,
                                  # width_ratios=(2, 1, 1, 1),
                                  height_ratios=(2, 3),
                                  # left=0.1, right=0.9, bottom=0.1, top=0.9,
                                  wspace=0.05, hspace=0.15)
        gs0 = gridspec.GridSpecFromSubplotSpec(len(dats), 1, height_ratios=[1] * len(dats), hspace=0.07,
                                               subplot_spec=gs[-1])
        ax_profs = fig.add_subplot(gs[0])
        ax_mass = fig.add_subplot(gs[1])
    elif orientation == 'vertical':
        gs = fig.add_gridspec(1, 2,
                              width_ratios=(5, 2),
                              hspace=0.05, wspace=0.15)
        gs0 = gridspec.GridSpecFromSubplotSpec(1, len(dats), width_ratios=[1] * len(dats), wspace=0.07,
                                               subplot_spec=gs[0])

        ax_profs = fig.add_subplot(gs[1])
    axes_comp = [fig.add_subplot(gs0[ii]) for ii in range(len(dats))]

    comp_ax_legend = True
    # compositions
    for ii, (dat, ax) in enumerate(zip(dats, axes_comp)):
        if ii == 1:
            comp_ax_ylabel = r'Mineral mode (wt.\%)'
        else:
            comp_ax_ylabel = ''
        if ii == 2:
            xlabel = 'Pressure (GPa)'
        else:
            xlabel = ''
        if comp_var == 'mgsi':
            comp_leg_title = r'$\bf{Mg/Si}$'
            prof_leg_label = '{:.1f}'.format(dat.mgsi)
            comp_ax_title = 'Mg/Si = ' + prof_leg_label
        elif comp_var == 'mgfe':
            comp_leg_title = r'$\bf{Mg\#}$'
            prof_leg_label = '{:.1f}'.format(dat.wt_oxides['MgO'] / (dat.wt_oxides['FeO'] + dat.wt_oxides['MgO']))
            comp_ax_title = r'Mg\# = ' + prof_leg_label
        elif comp_var == 'T':
            comp_leg_title = r'$\bf{T_p}$'
            prof_leg_label = '{:.1f}'.format(dat.temperature[-1])
            comp_ax_title = r'$T_p$ = ' + prof_leg_label
        else:
            comp_ax_title = ''
            comp_leg_title = ''
            prof_leg_label = comp_labels[ii]


        if orientation == 'horizontal':
            leg_bbox_to_anchor = (1.02, 2.1)
            comptitle = ''
        elif orientation == 'vertical':
            leg_bbox_to_anchor = (1.02, 0.5)
            comptitle = comp_ax_title
        fig, ax = plotpx.single_composition(dat, which='pressure', modality_type='phase', comp_stacked=True, save=False,
                                            ylabel=comp_ax_ylabel, xlabel=xlabel, cmap=cmap, override_ax_arrow=True,
                                            labelsize=labelsize, legsize=legsize,
                                            plot_phases_order=phases_order, p_max=p_max, orientation=orientation,
                                            fig=fig, ax=ax, title=comptitle, make_legend=comp_ax_legend, ylabelpad=10,
                                            xlabelpad=10,
                                            leg_bbox_to_anchor=leg_bbox_to_anchor, legtitle=r'$\bf{Phase}$')
        # ax.set_yticks([])
        if orientation == 'horizontal':
            # only label Mg/Si in bbox for horizontal, otherwise it's title
            ax.text(0.98, 0.4, comp_ax_title, transform=ax.transAxes, fontsize=legsize,  # weight='bold',
                    va='center', ha='right', c=['k', 'k', 'w'][ii],
                    bbox=dict(boxstyle='round', fc=colours_comp[ii], ec='k', alpha=0.4, ls=ls_comp[ii]))
        # ax.text(0.98, 0.4, comp_ax_title, transform=ax.transAxes, fontsize=legsize,
        #         va='center', ha='right', c='k',
        #         bbox=dict(boxstyle='round', fc=(0, 0, 0, 0), ec=colours_comp[ii], ls=ls_comp[ii], lw=2))

        if orientation == 'horizontal':
            profile_kwargs = {'reverse_y': False, 'xmin':1000e-4, 'xmax':p_max, 'ymin':ylim[0], 'ymax':ylim[1], 'label_x':False,
                              'label_y': True, 'leg_bbox_to_anchor':(1.02, 1.1), 'ax_ticks':None}
            if p_max is None:
                print('WARNING: possible problem setting axes limits with p_max = None')
        elif orientation == 'vertical':
            profile_kwargs = {'reverse_y': True, 'ymin': 1000e-4, 'ymax': p_max, 'xmin': ylim[0], 'xmax': ylim[1],
                              'label_x': True,
                              'label_y': False, #'leg_bbox_to_anchor': (1.02, 1.1)
                              'xlabelpad':10,
                              }
        fig, ax_profs = plotpx.profile(dat, profile_var, independent_ax='pressure', cumulative=cumulative,
                                       ax_label=prof_ax_label,
                                       scale=profile_scale,
                                       c=colours_comp[ii], ls=ls_comp[ii], lw=1.5, alpha=0.8,
                                       labelsize=labelsize, legsize=legsize, legtitle=comp_leg_title,
                                       fig=fig, ax=ax_profs, log=log_profile, leg_label=prof_leg_label,
                                       orientation=orientation, save=False, **profile_kwargs)

        if show_cum_mass:
            profile_kwargs2 = profile_kwargs.copy()
            profile_kwargs2['ymax'] = dat.M_p
            fig, ax_mass = plotpx.profile(dat, 'mass(kg)', independent_ax='pressure', cumulative=True,
                                           ax_label='Cumulative\nmass (kg)',
                                           scale=1, show_legend=False,
                                           c=colours_comp[ii], ls=ls_comp[ii], lw=1.5, alpha=0.8,
                                           labelsize=labelsize, legsize=legsize, legtitle=comp_leg_title,
                                           fig=fig, ax=ax_mass, log=log_profile, leg_label=prof_leg_label,
                                           orientation=orientation, save=False, **profile_kwargs2)
            ax_mass.set_xticks([])

        if orientation == 'horizontal':
            if ii < len(dats) - 1:
                ax.tick_params(axis="x", labelbottom=False)
                ax.set_xticks([])
            ax_profs.tick_params(axis="x", labelbottom=False)
        if orientation == 'vertical':
            if ii > 0:
                ax.tick_params(axis="y", labelleft=False)
                ax.set_yticks([])
            ax_profs.tick_params(axis="y", labelleft=False)
        comp_ax_legend = False  # only make legend once

        print('\n', dat.name, 'um water mass', dat.mass_h2o_um / p.TO, 'total water mass', dat.mass_h2o_total / p.TO,
              'TO', '\n')
        # print('c_ppv', dat.data['sat_corr_ppv'].iloc[-1] * 1e6, 'ppm', 'D_ppv_pv',
        #       dat.data['sat_corr_ppv'].iloc[-1] / dat.data['sat_corr_pv'].iloc[-1])

    if save:
        plt.tight_layout()
        fig.savefig(plotpx.fig_path + figtitle + fformat, bbox_inches='tight', dpi=300)


# demo_composition_gridspec(mgsi=[1.4125375446227555, 1.0715193052376069, 0.7244359600749921],
#                           comp_var='mgsi', labelsize=12, legsize=10, cmap=None, fformat='.pdf',
#                           figtitle='phase_demo_ver', orientation='vertical')

""" test with LM included """
# demo_composition_gridspec(mgsi=[1.4125375446227555, 1.0715193052376069, 0.7244359600749921],
#                           comp_var='mgsi', labelsize=12, legsize=10,  # cmap=None,
#                           save=True, p_max=300, figtitle='LM_mgsi_comp',
#                           phases_order=['gt', 'cpx', 'opx', 'hpcpx', 'ol', 'wad', 'ring', 'pv', 'qtz', 'coes',
#                                         'st', 'wus', 'capv', 'ppv', 'fapv'],  # for testing
#                           )

""" temperature """
# sample = [rw.read_name(output_path=px.output_parent_px,
#                        name='Earth_' + t + 'K') for t in ['1900', '1600']]
# star = 'HIP29295'
# sample = [rw.read_name(output_path=px.output_parent_px + 'hypatia1M_' + t + 'K_80Fe/',
#                        name='1M_80Ceff_' + star + '_' + t + 'K') for t in ['1900', '1600']]
# demo_composition_gridspec(dats=sample,
#                           comp_var='T', labelsize=12, legsize=10, log_profile=True, #cmap=None,
#                           save=True, p_max=150, figtitle='T_comp_test', phases_order=['gt', 'cpx', 'opx','hpcpx',  'ol', 'wad', 'ring','pv', 'qtz', 'coes',
#                                                                'st', 'wus', 'capv', 'ppv', 'fapv'] ,  # for testing
#                           )


""" Mg/Fe ? """
# # sample = [rw.read_name(output_path=px.output_parent_px + 'MgFe_from_earth/',
# #                        name='1M_32CMF_' + m + 'Mg_1600K') for m in ['47', '31', '29']]  # vary mantle FeO, constant CMF
# sample = [rw.read_name(output_path=px.output_parent_px + 'MgFe_from_sun/',
#                        name='1M_' + x + 'Ceff_sun_1600K') for x in ['99', '88', '70']]  # vary core_eff, constant bulk Fe
# sample = [rw.read_name(output_path=px.output_parent_px + 'bulk_Fe_test/',
#                        name='Fe_' + s) for s in ['min', 'max']]  # vary only total bulk Fe
# for dat in sample:
#     print(dat.core_eff_from_cmf())
#     plotpx.single_structure(dat, fig_path=dat.output_path)
# demo_composition_gridspec(dats=sample, comp_var='mgfe', labelsize=12, legsize=10, p_max=None, #cmap=None,
#                           save=True, figtitle='phases_vs_bulkFe',
#                           phases_order=['gt', 'cpx', 'ol', 'opx', 'qtz', 'coes',
#                                                                'hpcpx', 'wad', 'capv', 'ring', 'pv',
#                                                                'st', 'wus',  'ppv', 'fapv']
#                           )
#


plt.show()
