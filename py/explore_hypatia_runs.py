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
earth = rw.build_planet(M_p=1 * p.M_E, test_oxides=px.wt_oxides_Earth,
                        maxIter=30, tol=1e-4, n=800, Tp=Tp, core_efficiency=0.88,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Earth_Si_test', x_Si_core=9,  # 'Earth_' + str(Tp) + 'K',
                        )
# earth.get_garnet_composition()
# plotpx.single_composition(earth, which='pressure', modality_type='water', comp_stacked=True, save=True,
#                           show=True, cmap='tab20', labelsize=16, plot_phases_order=None, p_max=None, make_legend=True)
# plotpx.single_phase_subfig(earth, 'c_h2o', var_scale=1e6, var_log=True, vertical_pressure=False,
#                                            title=earth.name,
#                                            var_label='Water capacity\n(ppm)', save=False, show=True, phase_order=None,
#                                            labelsize=12, legsize=10, ticksize=12, xpad=10, ypad=10,
#                                            # ax_ticks=[1e1, 1e3, 1e5],
#                                            annotation='total = {0:3.1f} earth oceans'.format(earth.mass_h2o_total / p.TO),
#                                            cmap_phases='tab20', linec='xkcd:navy', linew=2, ymin=1e0, ymax=1e5)

""" check compositon and sat profile of individual cases """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia2/', subsample=10)
# stars = ['2MASS 19141179+3833548', 'HIP 86087', 'HIP 79126']
# dats_1600 = [rw.read_name(output_path=px.output_parent_default + 'hypatia1M/', star=st, M_p=p.M_E,
#                      core_efficiency=0.8, Tp=1600) for st in stars]
# dats_1900 = [rw.read_name(output_path=px.output_parent_default + 'hypatia1M_1900K/', star=st, M_p=p.M_E,
#                      core_efficiency=0.8, Tp=1900) for st in stars]

# dat0 = rw.read_name(output_path=px.output_parent_default + 'hypatia1M_1900K_60Fe/', star='2MASS 19375133+4945541', M_p=1*p.M_E, core_efficiency=0.6, Tp=1900)
# dat1 = rw.read_name(output_path=px.output_parent_default + 'hypatia1M_1600K_70Fe/', star='2MASS 19375133+4945541', M_p=1*p.M_E, core_efficiency=0.7, Tp=1600)
# #
# # # m_mtl = dat0.cum_mass[-1] - dat0.cum_mass[dat0.i_cmb]
# # # print('mass ratio UM/mantle', dat0.mass_um / dat0.M_p)
# for dat in (dat0, dat1):
#     print('pmax', dat.pressure[dat.i_cmb + 1] * 1e-9, 'GPa', 'Tmax', dat.temperature[dat.i_cmb + 1], 'K')
#     m_w_tot = dat.mass_h2o_total / p.TO
#     plotpx.single_composition(dat, which='pressure', modality_type='water',
#                        xlabel=None, ylabel=None, cmap='tab20', labelsize=16, plot_phases_order=['gt', 'cpx', 'opx','hpcpx',  'ol', 'wad', 'ring','pv', 'qtz', 'coes',
#                                                                'st', 'wus', 'capv', 'ppv', 'fapv'] , p_max=30,)
# fig, axes = plotpx.composition_subfig(dat, 'c_h2o', var_scale=1e6, var_log=False, vertical_pressure=False,
#                                        title='',
#                                       phase_order=['gt', 'cpx', 'ol', 'opx', 'hpcpx', 'wad', 'capv', 'ring', 'st', 'wus', 'pv', 'ppv', 'fapv'],
#                                        var_label='Water capacity\n(ppm)', save=True, show=False,
#                                        labelsize=12, legsize=10, ticksize=12, xpad=10, ypad=10,
#                                        # ax_ticks=[1e1, 1e3, 1e5],
#                                        annotation='total = {0:3.1f} earth oceans'.format(m_w_tot),
#                                        cmap_phases='tab20', linec='xkcd:navy', linew=2, ymin=1e0, ymax=6e2)


""" try some crossplots - compositional effects """
# earth.get_phase_masses()

# UM water vs. mg/si
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M/')
# plotpx.pop_scatter(dats, 'mgsi', 'mass_h2o_um', x_scale=1, y_scale=p.TO**-1, c='k', alpha=0.4,
#                         xlabel='Mg/Si', ylabel='Upper mantle water capacity (Earth oceans)', filename=None,
#                         earth=earth, save=True)

# # mass of phase vs. water
# dats = rw.update_dir(px.perplex_path_default + 'output/hypatia2/', px.PerplexData.get_phase_masses, store=False)
# plotpx.pop_scatter(dats[:500], "phase_mass['O']", "mass_h2o_um", x_scale=1, y_scale=p.TO**-1, c='k', alpha=0.4,
#                         ylabel='Upper mantle water capacity (Earth oceans)', xlabel='mass Ol (kg)', filename='mass_h2o_ol',
#                         earth=earth, save=True)

# mgsi effect on upper mantle mass
# dats = rw.update_dir(px.perplex_path_default + 'output/hypatia2M/', px.PerplexData.get_phase_masses, store=False)
# plotpx.pop_scatter(dats, "mgsi", "mass_um", c='xkcd:bordeaux', alpha=0.4,
#                         ylabel='Upper mantle mass (kg)', xlabel='Mg/Si',
#                         earth=earth, save=True)


""" upper mantle water mass - histograms """


def hist_saturation_subfig(which='um', masses=None, xlim='default', bins=None, labelsize=14, cmap='rainbow', save=True,
                           show=True, earth=None, **kwargs):
    if masses is None:
        masses = [0.1, 0.5, 1, 2, 3, 4]
    if which == 'um':
        if xlim == 'default':
            xlim = (0.1, 4.0)
        key = 'mass_h2o_um'
        xlabel = 'Mantle water capacity (Earth oceans)'
        data_label = None
        ls = '-'
        fname = 'histsubplot-Mp_w_um'
    elif which == 'total':
        if xlim == 'default':
            xlim = (0.1, 8)
        key = 'mass_h2o_total'
        xlabel = 'UM water capacity (Earth oceans)'
        data_label = None
        ls = '-'
        fname = 'histsubplot-Mp_w_tot'
    elif which == 'both':
        if xlim == 'default':
            xlim = (0.1, 8)
        key = ['mass_h2o_um', 'mass_h2o_total']
        xlabel = 'Water capacity (Earth oceans)'
        data_label = ['Upper mantle', 'Whole mantle']
        ls = ['-', '--']
        fname = 'histsubplot-Mp_w_all'

    fig, axes = plt.subplots(len(masses), 1)
    c = colorize(masses, cmap)[0]
    for ii, Mp in enumerate(masses):
        if isinstance(Mp, float):
            mass_str = str(Mp).replace('.', ',')
        elif isinstance(Mp, int):
            mass_str = str(Mp)
        else:
            print('MP is', type(Mp))
        planets = rw.read_dir(px.perplex_path_default + 'output/hypatia' + mass_str + 'M/')
        try:
            ax = axes[ii]
        except:
            ax = axes
        fig, ax = plotpx.pop_hist1D(planets, key, scale=p.TO ** -1, earth=None, xlabel='', c_hist=c[ii], ls_stats=ls,
                                    save=False, show=False, data_label=data_label, fig=fig, ax=ax, xlim=xlim,
                                    labelsize=labelsize, legsize=10, bins=bins, alpha=0.7, histtype='step', **kwargs)
        # add naive Earth scaling
        if earth is not None:
            w_um_earth = earth.mass_h2o_um / p.TO  # in earth oceans
            print('w UM earth', w_um_earth, 'scaled', w_um_earth * Mp)
        ax.axvline(w_um_earth * Mp, c='k', ls='--', lw=0.5, label=str(Mp) + r' $\times$ Earth value')
        ax.legend(fontsize=10, frameon=False)
        if ii == len(masses) - 1:
            ax.set_xlabel(xlabel)
        else:
            ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(str(Mp) + ' $M_\oplus$', fontsize=labelsize)
    # fig.suptitle('Mantle wantle water capacities, $T_p$ = 1600')
    if save:
        fig.savefig(plotpx.fig_path + fname + '.png')
    if show:
        plt.show()


# hist_saturation_subfig(which='um', masses=[0.1], xlim=(0, 1),
#                        bins=100, showmedian=True, show=False)


def find_percentile(q, dir_name, prop):
    if np.size(q) > 1:
        raise NotImplementedError('find_percentile() - multiple percentiles at once not implemented. use scalar q.')
    dats = rw.read_dir(px.perplex_path_default + 'output/' + dir_name + '/')
    x = []
    names = []
    test_vals = []
    for dat in dats:
        try:
            x.append(eval('dat.' + prop))
            names.append(dat.star)
        except AttributeError:
            pass
    pcen = np.percentile(x, q, method='nearest')
    i_near = abs(x - pcen).argmin()
    print(q, 'pcen', pcen, 'val')
    return names[i_near]


folder = 'hypatia1M_1600K_80Fe'
# print(find_percentile(2.27, folder, 'mgsi'))  # HIP 5643, 0.7244
# print(find_percentile(10, folder, 'mgsi'))  # HIP 29295, 0.8912509381337477
# print(find_percentile(48, folder, 'mgsi'))  # HIP 48455
# print(find_percentile(50, folder, 'mgsi'))  # 2MASS 19375133+4945541, 1.0715
# print(find_percentile(52, folder, 'mgsi'))  # 2MASS 19330262+4452080
# print(find_percentile(90, folder', 'mgsi'))
# print(find_percentile(97.73, folder, 'mgsi'))  # HIP 114933, 1.4125

""" c_h2o profiles with T dependence, for an average star"""

# star = '2MASS 19375133+4945541'
# dat = rw.read_name(px.perplex_path_default + 'output/hypatia1M/', star=star, M_p=1*p.M_E, core_efficiency=0.8, Tp=1600)
# fig, ax = plotpx.profile(dat, 'mass_h2o(kg)', independent_ax='pressure', reverse_y=True, ax_label='Water capacity (Earth oceans)', scale=p.TO**-1,
#                          c='xkcd:deep blue', lw=2, alpha=0.9,
#             xmin=None, xmax=None, ymin=0, ymax='max', ax_ticks=None, label_x=True, label_y=True, labelsize=14, legsize=12,
#             fig=None, ax=None, figsize=(4, 6), log=True, leg_label='1600 K', orientation='vertical',
#             y2var=None, y2label=None, y2scale=None, save=False)
#
# dat = rw.read_name(px.perplex_path_default + 'output/hypatia1M_1900K/', star=star, M_p=1*p.M_E, core_efficiency=0.8, Tp=1900)
# fig, ax = plotpx.profile(dat, 'mass_h2o(kg)', independent_ax='pressure', reverse_y=True, ax_label='Water capacity (Earth oceans)', scale=p.TO**-1,
#                          c='xkcd:peach', lw=2, alpha=0.9,
#             xmin=None, xmax=None, ymin=0, ymax='max', ax_ticks=None, label_x=True, label_y=True, labelsize=14, legsize=12,
#             fig=fig, ax=ax, log=True, leg_label='1900 K', orientation='vertical', save=True, fname='massh2o_profiles_hotcold')


""" compare Mg/Si vs. Mg/Fe - across hypatia but fixed core eff (underestimates varibility of FeO) """
# fig, axes = plt.subplots(1, 2, sharey=True)
# ylim = (1, 8)
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_80Fe/')
# fig, ax = plotpx.pop_scatter(dats, 'mgsi', 'mass_h2o_total', x_scale=1, y_scale=p.TO ** -1, c='b', alpha=0.2,
#                              xlabel='Mg/Si', ylabel='Mantle water capacity (Earth oceans)', fig=fig, ax=axes[0],
#                              save=False, show=False, ylim=ylim,
#                              # earth=earth,
#                              )
# dirs = [px.perplex_path_default + 'output/hypatia1M_1600K_' + s + 'Fe/' for s in ['70', '80', '90', '99']]
# fig, ax = plotpx.compare_pop_scatter(dirs, "wt_oxides['MgO']/dat.wt_oxides['FeO']", 'mass_h2o_total', x_scale=1,
#                                      y_scale=p.TO ** -1, c='r', alpha=0.2, lw=0,
#                                      xlabel='Mg/Fe', ylabel='Mantle water capacity (Earth oceans)', #fig=fig, ax=axes[1],
#                                      save=False, show=False,# ylim=ylim,
#                                      # earth=earth,
#                                      )

""" compare Mg/Si vs. Mg/Fe - scaling from solar """
# fig, axes = plt.subplots(1, 2, sharey=True)
# ylim = (1, 4)
# dats = rw.read_dir(px.perplex_path_default + 'output/MgSi_from_sun/')
# fig, ax = plotpx.pop_scatter(dats, 'mgsi', 'mass_h2o_total', y_scale=p.TO ** -1,
#                              c='b', alpha=0.5, xlabel='Mg/Si', ylabel='Mantle water capacity (Earth oceans)',
#                              save=False, show=False, ylim=ylim, fig=fig, ax=axes[0],
#                              # earth=earth,
#                              )
# dats = rw.read_dir(px.perplex_path_default + 'output/bulk_Fe_test/')  # MgFe_from_sun
# fig, ax = plotpx.pop_scatter(dats, "wt_oxides['MgO']/pl.wt_oxides['FeO']", 'mass_h2o_total', y_scale=p.TO ** -1,
#                              c='r', alpha=0.5, xlabel='Mg/Fe', ylabel='Mantle water capacity (Earth oceans)',
#                              save=False, show=False, ylim=ylim, fig=fig, ax=axes[1],  # xlim=(0.5, 8),
#                              # earth=earth,
#                              )

""" hist of star iron content """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_80Fe/')
# plotpx.pop_hist1D(dats, x_name, scale=1, earth=None, xlabel=None,
#                save=False, show=True, data_label=None,
#                xlim=None, bins=None, showmedian=True)


""" what is earth-composition gravity at min and max mass"""
# for m in [0.5]:
#     pl = rw.build_planet(M_p=m * p.M_E, test_CMF=0.325, test_oxides=px.wt_oxides_Earth,
#                          maxIter=30, tol=1e-4, n=800, Tp=1600, plot_all=False, get_saturation=False, verbose=False)
#     print('Mp', m, 'g', pl.gravity[-1], 'm/s^2')
# 0.1 M - 3.9167161840595828
# 2 - 13.704242549268521
# 3 - 16.525783688468746
# 4 M - 18.952818930786943
# 5 M - 21.13518499447704

""" histogram of comp """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia4M_1600K_88Fe/')
# # x_name = "df_all['X_wus'].iloc[-1] + pl.df_all['X_ppv'].iloc[-1] + pl.df_all['X_fapv'].iloc[-1]"
# x_name = "df_all['X_ppv'].iloc[-1]"
# plotpx.pop_hist1D(dats, x_name, scale=1, earth=None, xlabel=None, title=None, c_hist='k', ls='-',
#                filename=None, extension='.png', save=False, show=True, data_label=None, fig=None, ax=None,
#                xlim=None, labelsize=12, legsize=12, bins=100, showmedian=False, showsigma=True)


""" explain earth high water mass ? """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_88Fe/')
# params = ['c_h2o_mantle', 'mass_h2o_total', 'mgsi', "wt_oxides['MgO']/pl.wt_oxides['FeO']"]
# xlims = [(0, 0.004), (0, 1.5e22), None, None]
# xlabels = ['Mantle c_h2o', 'Mantle m_h2o (kg)', 'Mg/Si', 'Mantle Mg/Fe']
# fig, axes = plt.subplots(1, len(params))
# for ii, pm in enumerate(params):
#     fig, axes[ii] = plotpx.pop_hist1D(dats, pm, xlim=xlims[ii], bins=100, scale=1, earth=earth, ls_stats='--', c_hist='k', xlabel=xlabels[ii],
#                       histtype='step', showmedian=True, showsigma=True, save=False, show=False, annotate_n=False,
#                       fig=fig, ax=axes[ii])
# todo: make axis lims +/- 3 sigma to compare, and cut off very extreme low Mg/Si because skewing?


""" hist of FeO """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_70Fe/')
# # plotpx.pop_hist1D(dats, "wt_oxides['MgO']/(dat.wt_oxides['FeO'] + dat.wt_oxides['MgO'])", scale=1, earth=None, showmedian=True, save=False)
# plotpx.pop_hist1D(dats, "wt_oxides['CaO']/(dat.wt_oxides['Al2O3'])", scale=1, earth=None, showmedian=True, save=False, xlim=(0, 5), bins=100)

""" hist of p_CMB """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_65Fe/')
# plotpx.pop_hist1D(dats, 'pressure[pl.i_cmb]', earth=None, showmedian=True)

""" find highest p_cmb """
# plotpx.find_extrema('hypatia3,5M_1600K_88Fe', 'pressure[0]', get_min=False, get_max=True, n=1, output_base='output/', scale=1e-9)
# plotpx.find_extrema('hypatia1M_1600K_65Fe', "wt_oxides['FeO']", get_min=False, get_max=True, n=1, output_base='output/')  # max FeO is 040%


""" find n richest planets """
# plotpx.find_rich_planets(dir_name='hypatia3M', phase='st', n=10)
# plotpx.find_rich_planets(dir_name='hypatia4M_1600K_88Fe', phase='ppv', n=10, get_min=True)


plt.show()
