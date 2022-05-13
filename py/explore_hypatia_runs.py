import numpy as np
import parameters as p
import ask_hypatia as hyp
import perplexdata as px
import plot_perplex as plotpx
import main as rw
import parameters as p
import matplotlib.pyplot as plt
from useful_and_bespoke import cornertext, colorize
import saturation as sat

""" get earth benchmark """
Tp = 1600
earth = rw.build_planet(M_p=1 * p.M_E, test_CMF=0.325, test_oxides=px.wt_oxides_Earth,
                        maxIter=30, tol=1e-4, n=800, Tp=Tp, #core_efficiency=0.8,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        # name='Earth_' + str(Tp) + 'K',
                        )
earth.get_garnet_composition()
# plotpx.single_composition(earth, which='pressure', modality_type='water', comp_stacked=True, save=True,
#                           show=True, cmap='tab20', labelsize=16, plot_phases_order=None, p_max=None, make_legend=True)

""" check a few cases and plot instantly to verify """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia2/', subsample=10)
# for dat in dats:
#     m_w_tot = sat.total_water_frac(dat.df_all) * dat.M_p / p.TO
#     fig, axes = plotpx.single_phase_subfig(dat, 'c_h2o', var_scale=1e6, var_log=True, vertical_pressure=False,
#                                            title=dat.name,
#                                            var_label='Water capacity\n(ppm)', save=False, show=True, phase_order=None,
#                                            labelsize=12, legsize=10, ticksize=12, xpad=10, ypad=10,
#                                            # ax_ticks=[1e1, 1e3, 1e5],
#                                            annotation='total = {0:3.1f} earth oceans'.format(m_w_tot),
#                                            cmap_phases='tab20', linec='xkcd:navy', linew=2, ymin=1e0, ymax=1e5)

""" try some crossplots """
# earth.get_phase_masses()

# UM water vs. mg/si
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia2/')
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

""" plot variable across masses with composition continuity """
# dirs = [px.perplex_path_default + 'output/' + s + '/' for s in ('hypatia1M', 'hypatia2M', 'hypatia4M')]
# # plotpx.compare_pop_scatter(dirs, 'M_p', 'mass_um', x_scale=p.M_E ** -1, y_scale=1, xlog=False, ylog=False,
# #                            ylabel='Upper mantle mass (kg)', xlabel='$M_p$ ($M_\oplus$)', title='all Mg/Si',
# #                            filename='pops_M_p_mass_um', save=True, show=True, extension='.png',
# #                            labelsize=12, lw=0.2, vmin=0.2, vmax=2.5, c='rainbow', earth=earth, v_name='mgsi', alpha=1)
# plotpx.compare_pop_fillbetween(dirs, 'M_p', 'mass_um', x_scale=p.M_E ** -1, y_scale=1, xlog=False, ylog=False,
#                            ylabel='Upper mantle mass (kg)', xlabel='$M_p$ ($M_\oplus$)', sigma=2,
#                             save=True, show=True, extension='.png',
#                            labelsize=12, earth=earth)
# plotpx.compare_pop_fillbetween(dirs, 'M_p', 'mass_h2o_um', x_scale=p.M_E ** -1, y_scale=p.TO ** -1, xlog=False, ylog=False,
#                            ylabel='UM water capacity (Earth oceans)', xlabel='$M_p$ ($M_\oplus$)', sigma=2,
#                             save=True, show=True, extension='.png',
#                            labelsize=12, earth=earth)

""" upper mantle water mass - copy back to water dist hypatia"""
# fig, axes = plt.subplots(3, 1)
# labelsize = 14
# masses = [1, 2, 4]
# c = colorize(masses, 'rainbow')[0]
# for ii, Mp in enumerate(masses):
#     planets = rw.read_dir(px.perplex_path_default + 'output/hypatia' + str(Mp) + 'M/')
#     fig, ax = plotpx.pop_hist1D(planets, 'mass_h2o_um', scale=p.TO ** -1, earth=None,
#                                 xlabel='', c_hist=c[ii], ls='-',
#                                 filename='sat_um_exohosts_4M', save=False, show=False, data_label=str(Mp) + ' $M_\oplus$',
#                                 alpha=0.7, xlim=(0.5, 3.5),
#                                 histtype='step', fig=fig, ax=axes[ii], labelsize=labelsize)
# axes[-1].set_xlabel('UM water capacity (Earth oceans)')
# fig.suptitle('Upper mantle wantle water capacities, $T_p$ = 1600')

# with lower mantle
# plotpx.pop_hist1D(planets, ['mass_h2o_um', 'mass_h2o_total'], scale=TO ** -1, earth=earth, xlabel='Water capacity (Earth oceans)',
#                   title='Mantle water capacities, ' + str(Mp) + ' $M_\oplus$, $T_p$ = ' + str(Tp),
#                     c_hist='k', ls=['-', '--'],
#                   filename='sat_um_exohosts', save=True, show=True, data_label=['Upper mantle', 'Whole mantle'], alpha=0.7,
#                   histtype='step')


""" find n richest planets """
# plotpx.find_rich_planets(dir_name='hypatia2M', phase='st', n=30)

plt.show()
