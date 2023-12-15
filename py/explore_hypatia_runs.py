import numpy as np
import parameters as p
import ask_hypatia as hyp
import perplexdata as px
import plot_perplex as plotpx
import saturation as sat
import main as rw
import matplotlib.pyplot as plt
import pandas as pd
from useful_and_bespoke import cornertext, colorize, colourbar, colourised_legend
from UM_mass_scaling import um_mass_scaling
import saturation as sat
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
import pickle as pkl
from scipy import stats

output_parent_path_aopp = '/home/g/guimond/Work/hypatia-compositions/'
perplex_path_aopp = '/home/g/guimond/Work/perple_x/'

# dat = rw.build_planet(M_p=1 * p.M_E, test_oxides=rw.update_MgSi(0.7244359600749921, px.wt_oxides_MD95),
#                       maxIter=30, tol=1e-4,  # n=800,
#                       Tp=1600, test_CMF=0.325,
#                       plot_all=False, get_saturation=True, verbose=True, clean=True,
#                       vertex_data='stx21ver', option_file='perplex_option_claire_mol', excluded_phases=[],
#                       name='1M_32CMF_27Mg_1600K_mol',
#                       output_parent_path=px.output_parent_px + 'MgSi_from_earth_fig1/',
#                       )
# dat.werami_garnet_composition()

""" get earth benchmark """
Tp = 1900
for m in [1] : #[0.1, 0.3, 0.5, 2, 2.5, 3, 4, 5]:
    earth = rw.build_planet(M_p=m * p.M_E, test_oxides=px.wt_oxides_MD95,
                            maxIter=30, tol=1e-4, #n=800,
                            Tp=Tp, test_CMF=0.325, #core_efficiency=0.88,
                            plot_all=False, get_saturation=False, verbose=True, clean=True,
                            vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                            perplex_path=perplex_path_aopp, output_parent_path=output_parent_path_aopp,
                            name='Earth_' + str(Tp) + 'K_' + str(m).replace('.', ',') + 'M',
                            )
# print('earth CMF', earth.CMF, 'core eff', earth.core_eff)
# earth.find_lower_mantle()
# print('earth mass um', earth.mass_um, 'kg')
# print('earth mgsi', earth.mgsi)
# earth.femg_star = 0.81


# Mp = 2.25
# pl_E = rw.build_planet(M_p=Mp * p.M_E, test_oxides=px.wt_oxides_MD95,
#                         maxIter=30, tol=1e-4, #n=800,
#                         Tp=Tp, test_CMF=0.325, #core_efficiency=0.88,
#                         plot_all=False, get_saturation=True, verbose=True, clean=True,
#                         vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                         # name='Earth_Si_test', x_Si_core=9,
#                         name='test1_' + str(Tp) + 'K',
#                         )
# rho = pl_E.M_p / (4/3 * np.pi * pl_E.R_p**3)
# print('Earth-like bulk density', rho, 'kg/m3')  # 6562 kg/m3
#
# wt_oxides_CaAl = {'SiO2': 27, 'MgO': 14, 'CaO': 16, 'Al2O3': 43, 'FeO': 0.1}
# pl_Ca = rw.build_planet(M_p=Mp * p.M_E, test_oxides=wt_oxides_CaAl,
#                         maxIter=30, tol=1e-4, #n=800,
#                         Tp=Tp, test_CMF=0.325, #core_efficiency=0.88,
#                         plot_all=False, get_saturation=True, verbose=True, clean=True,
#                         vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                         # name='Earth_Si_test', x_Si_core=9,
#                         name='test2_' + str(Tp) + 'K',
#                         )
# rho = pl_Ca.M_p / (4/3 * np.pi * pl_Ca.R_p**3)
# print('bulk density', rho, 'kg/m3')  # 6445 kg/m3


# earth.get_obm_water()
# m_w_obm = earth.mass_h2o_obm
# m_w_um = earth.mass_h2o_um
# m_w_mtz = m_w_um - m_w_obm
# print('earth obm mass', m_w_obm/p.TO, 'OM', 'mtz mass', m_w_mtz/p.TO, 'um mass', m_w_um/p.TO)
# i_lm = earth.find_lower_mantle()
# i_mtz = earth.find_transition_zone()

# sun = rw.build_planet(M_p=1 * p.M_E,
#                         maxIter=30, tol=1e-4, n=800, Tp=Tp, core_efficiency=0.88, star='sun',
#                         plot_all=False, get_saturation=True, verbose=True, clean=True,
#                         vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                         # name='Earth_Si_test', x_Si_core=9,
#                         name='sun_' + str(Tp) + 'K',
#                         )
# print('solar CMF', sun.CMF, 'core eff', sun.core_eff)
#

# earth.get_garnet_composition()
# plotpx.single_composition(earth, denominator='pressure', modality_type='water', comp_stacked=True, save=True,
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
# dats_1600 = [rw.read_name(output_path=px.output_parent_px + 'hypatia1M/', star=st, M_p=p.M_E,
#                      core_efficiency=0.8, Tp=1600) for st in stars]
# dats_1900 = [rw.read_name(output_path=px.output_parent_px + 'hypatia1M_1900K/', star=st, M_p=p.M_E,
#                      core_efficiency=0.8, Tp=1900) for st in stars]

# dat0 = rw.read_name(output_path=px.output_parent_px + 'hypatia1M_1900K_60Fe/', star='2MASS 19375133+4945541', M_p=1*p.M_E, core_efficiency=0.6, Tp=1900)
# dats = [rw.read_name(output_path=px.output_parent_px + 'hypatia0,1M_1600K_88Fe/', star=s, M_p=0.1*p.M_E,
#                      core_efficiency=0.88, Tp=1600) for s in ['2MASS 19394601-2544539',
# '2MASS 19141179+3833548', 'HIP 24186', 'HIP 79126',
#                                                               'HD 240210', 'HIP 1475', 'HIP 86087', 'HIP 80459',
#                                                               #'HIP 86287', 'HIP 114046', 'HIP 84460'
#                                                               ]]
# star = '2MASS19375133+4945541'#'HIP1692'
# dats = [#rw.read_name(output_path=px.output_parent_px + '/earthsize_planets_1600K_88Fe/', name='Kepler-68c')
#     rw.read_name(output_path=px.output_parent_px + 'MgSi_from_Earth/', name='1M_70Ceff_' + star + '_1600K'),
#         rw.read_name(output_path=px.output_parent_px + 'hypatia1M_1600K_88Fe/', name='1M_88Ceff_' + star + '_1600K'),
#     #     rw.read_name(output_path=px.output_parent_px + 'hypatia1M_1600K_99Fe/', name='1M_99Ceff_' + star + '_1600K')
#         ]
# dats = [rw.build_planet(M_p=0.1 * p.M_E,
#                         maxIter=30, tol=1e-4, n=1200, Tp=1900, core_efficiency=0.9999, star='HIP 24186',
#                         plot_all=False, get_saturation=True, verbose=True, clean=True, use_local_composition=True,
#                         vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                         )]
# # # # # m_mtl = dat0.cum_mass[-1] - dat0.cum_mass[dat0.i_cmb]
# # # # # print('mass ratio UM/mantle', dat0.mass_um / dat0.M_p)
# for dat in dats:
#     # print('p_mtz', dat.p_mtz * 1e-9, 'GPa')
#     # print('pmax', dat.pressure[dat.i_cmb + 1] * 1e-9, 'GPa', 'Tmax', dat.temperature[dat.i_cmb + 1], 'K')
#     m_w_tot = dat.mass_h2o_total / p.TO
#     # plotpx.single_composition(dat, denominator='pressure',# modality_type='water',
#     #                    xlabel=None, ylabel=None, cmap='tab20', labelsize=16,
#     #                           # plot_phases_order=['gt', 'cpx', 'opx','hpcpx',  'ol', 'wad', 'ring','pv', 'qtz', 'coes',
#     #                           #                                  'st', 'wus', 'capv', 'ppv', 'fapv']
#     #                           )
#     dat.find_lower_mantle()
#     # print(dat.p_lm * 1e-9, 'GPa')
#     fig, axes = plotpx.composition_subfig(dat, 'c_h2o', var_scale=1e6, var_log=False, vertical_pressure=False,
#                                            title='core eff = ' + str(dat.core_eff),
#                                           phase_order=['qtz', 'gt', 'cpx', 'ol', 'opx', 'hpcpx', 'wad', 'capv', 'ring', 'st', 'wus', 'pv', 'ppv', 'fapv'],
#                                            var_label='Water capacity\n(ppm)', save=True, show=False,
#                                            labelsize=12, legsize=10, ticksize=12, xpad=10, ypad=10,
#                                            # ax_ticks=[1e1, 1e3, 1e5],
#                                            annotation='total = {0:3.1f} earth oceans'.format(m_w_tot),
#                                            cmap_phases='tab20', linec='xkcd:navy', linew=2, ymin=1e0, ymax=50e2,# p_max=30,
#                                           )
# #     axes[1].axvline(dat.p_lm*1e-9, c='k', lw=2)

""" check individual mgsi """
# test1 = rw.build_planet(M_p=1 * p.M_E, test_oxides=rw.update_MgSi(0.72, px.wt_oxides_Earth),
#                         maxIter=30, tol=1e-4, n=800, Tp=Tp, core_efficiency=0.88,
#                         plot_all=False, get_saturation=True, verbose=True, clean=True,
#                         vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                         )
#
#
# test2 = rw.build_planet(M_p=1 * p.M_E, test_oxides=rw.update_MgSi(0.801, px.wt_oxides_Earth),
#                         maxIter=30, tol=1e-4, n=800, Tp=Tp, core_efficiency=0.88,
#                         plot_all=False, get_saturation=True, verbose=True, clean=True,
#                         vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                         )
# for mgsi in [1.4, 1.5, 1.6, 1.7]:
#     dat = rw.build_planet(M_p=1 * p.M_E, test_oxides=rw.update_MgSi(mgsi, px.wt_oxides_MD95),
#                         n=800, Tp=1600, core_efficiency=0.88,
#                         plot_all=False, get_saturation=True, verbose=True, clean=True,
#                         vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                         )
#     m_w_tot = dat.mass_h2o_total / p.TO
#     # plotpx.single_composition(dat, denominator='pressure', #modality_type='water',
#     #                    xlabel=None, ylabel=None, cmap='tab20', labelsize=16, p_max=35, title=str(dat.mgsi),
#     #                           plot_phases_order=['gt', 'cpx', 'opx','hpcpx',  'ol', 'wad', 'ring','pv',
#     #                                              # 'qtz', 'coes','st', 'seif',
#     #                                              'wus', 'capv', 'ppv', 'fapv'] ,
#     #                           )
#     fig, axes = plotpx.composition_subfig(dat, 'c_h2o', var_scale=1e6, var_log=False, vertical_pressure=False,
#                                            title='Mg/Si = ' + str(dat.mgsi),
#                                           plot_phases_order=['gt', 'cpx', 'opx', 'hpcpx', 'ol', 'wad', 'ring', 'pv',
#                                                              # 'qtz', 'coes','st', 'seif',
#                                                              'wus', 'capv', 'ppv', 'fapv'],
#                                            var_label='Water capacity\n(ppm)', save=True, show=False,
#                                            labelsize=12, legsize=10, ticksize=12, xpad=10, ypad=10,
#                                            # ax_ticks=[1e1, 1e3, 1e5],
#                                            # annotation='total = {0:3.1f} earth oceans'.format(m_w_tot),
#                                            cmap_phases='tab20', linec='xkcd:navy', linew=2, ymin=1e0, ymax=100e2, p_max=35,
#                                           )
# print(test1.get_phase_masses())
# print(test2.get_phase_masses())


""" try some crossplots - compositional effects """
# earth.get_phase_masses()

# UM water vs. mg/si
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_88Fe/')
# plotpx.pop_scatter(dats, 'mg_number', 'mass_um', x_scale=1, y_scale=p.TO**-1, c='g', alpha=0.4,
#                         # xlabel='Mg/Si', ylabel='Upper mantle water capacity (Earth oceans)', filename=None,
#                         earth=earth, save=True)

# stellar fe/mg
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_88Fe/')
# plotpx.pop_scatter(dats, 'femg_star', 'mass_um', x_scale=1, y_scale=p.TO**-1, c='g', alpha=0.4,
#                         # xlabel='Mg/Si', ylabel='Upper mantle water capacity (Earth oceans)', filename=None,
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

# plotpx.hist_saturation_subfig(denominator='um', masses=[0.1], xlim=(0, 1),
#                        bins=100, showmedian=True, show=False)

def find_sd(nsd, dir_name, prop, dats=None):
    # if np.size(q) > 1:
    #     raise NotImplementedError('find_percentile() - multiple percentiles at once not implemented. use scalar q.')
    if dats is None:
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
    mean = np.mean(x)
    std = np.std(x)
    pcen = (mean - nsd * std, mean + nsd * std)
    print(prop, pcen)
    return pcen


def find_percentile(q, dir_name, prop, dats=None):
    # if np.size(q) > 1:
    #     raise NotImplementedError('find_percentile() - multiple percentiles at once not implemented. use scalar q.')
    if dats is None:
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
    for qq, p in zip(q, pcen):
        i_near = abs(x - p).argmin()
        print(prop, qq, '%:', p, '@', names[i_near])
    return pcen


# folder = 'hypatia1M_1600K_80Fe'
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
# x_name = 'nH_star[-1]'
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1900K_90Fe/')
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
# # x_name = "data['X_wus'].iloc[-1] + pl.data['X_ppv'].iloc[-1] + pl.data['X_fapv'].iloc[-1]"
# x_name = "data['X_ppv'].iloc[-1]"
# plotpx.pop_hist1D(dats, x_name, scale=1, earth=None, xlabel=None, title=None, c_hist='k', ls='-',
#                filename=None, extension='.png', save=False, show=True, data_label=None, fig=None, ax=None,
#                xlim=None, labelsize=12, legsize=12, bins=100, showmedian=False, showsigma=True)


""" explain earth high water mass ? """

# folder = 'hypatia1M_1600K_88Fe'
# dats = rw.read_dir(px.perplex_path_default + 'output/' + folder + '/')
# params = ['mass_um', 'mass_h2o_um', "wt_oxides['MgO']/dat.wt_oxides['FeO']", 'p_mtz', 'p_lm'] #]#, 'mgsi', 'femg_star']
# # # xlims = [(0, 0.004), (0, 1.5e22), None, None]
# xlabels = ['UM mass (kg)', 'UM m_h2o (kg)', 'Mg/Fe', 'p mtz (Pa)', 'p LM (Pa)'] #Mantle Mg/Fe', 'Stellar Fe/Mg']
# for jj, pl in enumerate([earth]):#, sun)):
#     fig, axes = plt.subplots(len(params), 1)
#     for ii, pm in enumerate(params):
#         # get percentiles
#         # pmin, pmax = find_percentile([2.27, 97.73], folder, pm, dats=dats)
#         # pmin, pmax = find_sd(2, folder, pm, dats=dats)
#         fig, axes[ii] = plotpx.pop_hist1D(dats, pm,  #xmin=pmin, xmax=pmax,xlim=(pmin, pmax),
#                                           bins=50, scale=1, earth=pl, ls_stats='--', c_hist='k', alpha=0.5, xlabel=xlabels[ii],
#                           histtype='step', showmedian=True, showsigma=1, save=False, show=False, annotate_n=False,
#                           fig=fig, ax=axes[ii])
#     plt.tight_layout()
#     fig.savefig(plotpx.fig_path + ['earth', 'sun'][jj] + '_hists' + str(jj) + '.png')

def get_percentile_of_attr(target_dat=None, attr=None, x_target=None, dats=None, folder=None):
    from scipy import stats
    if dats is None:
        dats = rw.read_dir(px.perplex_path_default + 'output/' + folder + '/')
    x = []
    for dat in dats:
        try:
            x.append(eval('dat.' + attr))
        except (AttributeError, KeyError):
            pass  # blank data, didn't work because star not measured maybe
    if target_dat is not None:
        x_target = eval('target_dat.' + attr)
    p = stats.percentileofscore(x, x_target, 'rank')
    return p


# print(get_percentile_of_attr(attr='mgsi', x_target=0.81, dats=rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_88Fe/')))
#
# for ii, pm in enumerate(params):
#     print('earth', pm, get_percentile_of_attr(earth, pm, dats=dats), 'percentile')
#     print('sun', pm, get_percentile_of_attr(sun, pm, dats=dats), 'percentile')

""" hist of FeO """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_88Fe/')
# plotpx.pop_hist1D(dats, "femg_star", scale=1, earth=earth, showmedian=True, save=False)
# plotpx.pop_hist1D(dats, "wt_oxides['MgO']/(dat.wt_oxides['FeO'] + dat.wt_oxides['MgO'])", scale=1, earth=earth, showmedian=True, save=False)
# plotpx.pop_hist1D(dats, "wt_oxides['CaO']/(dat.wt_oxides['Al2O3'])", scale=1, earth=None, showmedian=True, save=False, xlim=(0, 5), bins=100)

""" hist of p_CMB """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_65Fe/')
# plotpx.pop_hist1D(dats, 'pressure[pl.i_cmb]', earth=None, showmedian=True)

""" hist of obm """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia3M_1600K_88Fe/')
# plotpx.pop_hist1D(dats, 'c_h2o_obm', earth=None, showmedian=True, showsigma=2)
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M_1600K_88Fe/')
# plotpx.pop_hist1D(dats, 'c_h2o_obm', earth=None, showmedian=True, showsigma=2)

""" hist of wmf """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia2M_1600K_88Fe/')
# plotpx.pop_hist1D(dats, 'mass_h2o_total', scale=1/(2*p.M_E) * 1e6, bins=100, range=(0, 500), earth=None, showmedian=True)

""" find extrema """
# plotpx.find_extrema('hypatia3,5M_1600K_88Fe', 'pressure[0]', get_min=False, get_max=True, n=1, output_base='output/', scale=1e-9)
# plotpx.find_extrema('hypatia1M_1600K_65Fe', "wt_oxides['FeO']", get_min=False, get_max=True, n=1, output_base='output/')  # max FeO is 040%
# plotpx.find_extrema('hypatia1M_1600K_99Fe', "mass_um", get_min=False, get_max=True, n=10, output_base='output/')
# plotpx.find_extrema('coreSi5_coreeff0,88', "wt_oxides['CaO']", get_min=False, get_max=True, n=10,
#                     perplex_path=perplex_path_aopp, output_base='output/')
# plotpx.find_extrema('coreSi5_coreeff0,88', "wt_oxides['Al2O3']", get_min=False, get_max=True, n=10,
#                     perplex_path=perplex_path_aopp, output_base='output/')

""" find n richest planets """
# plotpx.find_rich_planets(dir_name='hypatia3M', phase='st', n=10)
# plotpx.find_rich_planets(dir_name='hypatia4M_1600K_88Fe', phase='ppv', n=10, get_min=True)


""" update for unthought-of parameter """
# # # dirs = [px.perplex_path_default + 'output/hypatia' + s + 'M_1600K_88Fe/' for s in ['0,1', '0,5', '1','1,5', '2', '3']]  # todo other masses
# dirs = [px.perplex_path_default + 'output/hypatia' + s + 'M_1600K_88Fe/' for s in ('0,1', '0,3', '0,5', '1', '1,5', '2', '2,5', '3', '3,5', '4', '5')]
# # dirs = [px.perplex_path_default + 'output/hypatia1M_1600K_' + s + 'Fe/' for s in ['65', '70', '75', '80', '85', '90', '93', '95', '97', '99']]
# dirs = [px.perplex_path_default + 'output/hypatia0,1M_1600K_88Fe/']  # 'output/hypatia1M_1900K_88Fe/'
# for folder in dirs:
# dats = rw.update_dir(folder, px.PerplexData.write_star_composition, store=True)
# dats = rw.update_dir(folder, px.PerplexData.get_femg_star, store=True)
#     # dats = rw.update_dir(folder, px.PerplexData.get_obm_water, store=True)
#     dats = rw.update_dir(folder, px.PerplexData.find_lower_mantle, store=True)
#     # dats = rw.update_dir(folder, px.PerplexData.find_transition_zone, store=True)
#     dats = rw.update_dir(folder, px.PerplexData.get_um_mass, store=True)


""" depth to lm """
# dirs = [px.perplex_path_default + 'output/apollo/hypatia' + mp + 'M_1600K_88Fe_hires/' for mp in ['1', '2', '3', '4']]
# # for folder in dirs:
# #     dats = rw.update_dir(folder, px.PerplexData.find_lower_mantle, store=True)
#
# fig, ax = plotpx.compare_pop_scatter(dirs, "M_p", 'z_lm', x_scale=p.M_E**-1,
#                                      y_scale=1e-3,
#                                      c='r', alpha=0.2, lw=0,
#                                      xlabel='Mp', ylabel='z lm (km)', #fig=fig, ax=axes[1],
#                                      save=False, show=False,# ylim=ylim,
#                                      # earth=earth,
#                                      )


""" mass of ringwoodite as a cross-plot with Mg/Si for different mantle iron partitionings """

# dirs = [output_parent_path_aopp + 'hypatia1M_1600K_' + fe + 'Fe_hires/' for fe in ['70', '99']]
# # for folder in dirs:
# #     dats = rw.update_dir(folder, px.PerplexData.get_phase_masses, store=True)
# fig, ax = None, None
# for directory, colour, fe in zip(dirs, ['r', 'b'],['70', '99'] ):
#     dats = rw.read_dir(directory, verbose=True)
#     fig, ax = plotpx.pop_scatter(dats, "mgsi", "phase_mass['Ring']", x_scale=1, y_scale=(p.M_E**-1) * 1e2,
#                                  xlabel='Mg/Si', ylabel='mass ring (wt pt of planet)', range_min=None, range_max=None,
#                                  data_label=fe, fig=fig, ax=ax, xlim=None, ylim=None, save=False, show=False,
#                                  c=colour, annotate_n=True, ms=50, )

""" as above but ol/opx as a cross-plot with Fe/Si for different mantle iron partitionings - just at single pressure"""

# dirs = ['/home/g/guimond/Work/perple_x/output/coreSi0_coreeff' + fe + '/' for fe in ['0,1', '0,5', '0,88', '0,999']]
# # for folder in dirs:
# #     dats = rw.update_dir(folder, px.PerplexData.get_femg_si_mantle, store=True)
# fig, ax = None, None
# fig2, ax2 = None, None
# fig3, ax3 = None, None
# for directory, colour, fe in zip(dirs, ['r', 'b', 'g', 'k'], ['10', '50', '70', '99'] ):
#     dats = rw.read_dir(directory, verbose=True)
#     fig, ax = plotpx.pop_scatter(dats, "fesi_planet", "df_comp['O']/dat.df_comp['Opx']", x_scale=1, y_scale=1,
#                                  xlabel='Fe/Si', ylabel='Ol/Opx (mass)', range_min=None, range_max=None,
#                                  data_label=fe, fig=fig, ax=ax, xlim=None, ylim=(0, 80), save=False, show=False,
#                                  c=colour, annotate_n=True, ms=20, alpha=0.1 )
#     fig2, ax2 = plotpx.pop_scatter(dats, "mgsi", "df_comp['O']/dat.df_comp['Opx']", x_scale=1, y_scale=1,
#                                  xlabel='Mg/Si', ylabel='Ol/Opx (mass)', range_min=None, range_max=None,
#                                  data_label=fe, fig=fig2, ax=ax2, xlim=None, ylim=(0, 80), save=False, show=False,
#                                  c=colour, annotate_n=True, ms=20, alpha=0.1)
#     fig3, ax3 = plotpx.pop_scatter(dats, "femg_si_mantle", "df_comp['O']/dat.df_comp['Opx']", x_scale=1, y_scale=1,
#                                  xlabel='(Mg+Fe)/Si', ylabel='Ol/Opx (mass)', range_min=None, range_max=None,
#                                  data_label=fe, fig=fig3, ax=ax3, xlim=None, ylim=(0, 80), save=False, show=False,
#                                  c=colour, annotate_n=True, ms=20, alpha=0.1)

""" trends in (Ca+Al)/Si """
# dirs = ['/home/g/guimond/Work/perple_x/output/coreSi' + si + '_coreeff0,88/' for si in ['0', '5']]
# # for folder in dirs:
# #     dats = rw.update_dir(folder, px.PerplexData.get_minor_oxide_ratio, store=True)
# fig2, ax2 = None, None
# for directory, colour, fe in zip(dirs, ['g', 'k'], ['10%', '50%'] ):
#     dats = rw.read_dir(directory, verbose=True)
#     fig2, ax2 = plotpx.pop_scatter(dats, 'caal_si', "df_comp['Gt']", x_scale=1, y_scale=1,
#                                  xlabel='(Ca+Al)/(Si)', ylabel='Gt (wt%)', range_min=None, range_max=None,
#                                  data_label=fe, fig=fig2, ax=ax2, xlim=None, ylim=(0, 80), save=False, show=False,
#                                  c=colour, annotate_n=True, ms=20, alpha=0.1)

# folder = '/home/g/guimond/Work/hypatia-compositions/hypatia1M_1600K_88Fe_hires/'
# # dats = rw.update_dir(folder, px.PerplexData.get_um_mass, store=True)
# dats = rw.update_dir(folder, px.PerplexData.get_phase_masses, store=False, um_only=True)
# fig2, ax2 = plotpx.pop_scatter(dats, 'caal_si', "phase_mass['Gt']/dat.mass_um", x_scale=1, y_scale=1e2,
#                              xlabel='(Ca+Al)/(Si)', ylabel='Gt (wt% of UM)', range_min=None, range_max=None,
#                               xlim=None, #ylim=(0, 80),
#                                save=False, show=False,
#                              c='k', annotate_n=True, ms=20, alpha=0.1)
# hard to just look at total mass of garnet because also depends on mg/si (ol in transition zone)

""" mineralogy of low mg/si planets """
# directory = 'hypatia1M_1600K_70Fe_hires'
# star = 'HIP57087'
# plotpx.find_extrema(directory, "mgsi", get_min=True, get_max=False, n=20, output_base=output_parent_path_aopp, perplex_path='')
# with open(output_parent_path_aopp + directory + '/1M_70Ceff_' + star + '_1600K/dat.pkl', 'rb') as f:
#     dat = pkl.load(f)
#     plotpx.single_composition(dat,
#                        xlabel=None, ylabel=None, cmap='tab20', labelsize=16, p_max=35, title='Mg/Si = ' + str(dat.mgsi),
#                               plot_phases_order=['Gt', 'Cpx', 'Opx', 'C2/c',  'O', 'Wad', 'Ring','Pv',
#                                                  'qtz', 'coes','st', 'seif',
#                                                  'Wus', 'ca-pv', 'Ppv', 'fapv'] ,
#                               save=False, show=False)


""" enquire C/O ratio"""
def abundance_ratios_from_hypatia_tsv(path_to_tsv='/home/g/guimond/Work/hypatia-compositions/hypatia-04122023.tsv',
                                      el1=None, el2=None, xlim=None, percentile_of=None):
    df = pd.read_csv(path_to_tsv, sep='\t', header=0)
    df.replace('', np.nan, inplace=True)  # first convert empty string to nan
    df.dropna(subset=[el1, el2], inplace=True)

    lognH_1 = df[el1]
    lognH_2 = df[el2]

    sol_val1 = eval('p.' + el1.lower() + '_sol')
    sol_val2 = eval('p.' + el2.lower() + '_sol')

    ratio = (10 ** (lognH_1 + sol_val1) / 10 ** (lognH_2 + sol_val2)).to_numpy()

    plt.hist(ratio, histtype='step', bins=50, facecolor='k', edgecolor='k')
    plt.xlabel(el1 + '/' + el2)
    # plt.gca().set_xlim(xlim)

    if percentile_of is not None:
        p = stats.percentileofscore(ratio, percentile_of, 'rank')
        count = len([i for i in ratio if i >= percentile_of])
        print(percentile_of, 'is', p, 'percentile...', count, '/', len(ratio), '=', count/len(ratio)*100, '% have >= C/O')
    return ratio


# abundance_ratios_from_hypatia_tsv(el1='C', el2='O', xlim=(0, 1.4), percentile_of=0.9)

plt.show()
