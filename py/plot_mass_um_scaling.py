import numpy as np
import parameters as p
from perplexdata import wt_oxides_Earth, wt_oxides_MD95, perplex_path_default
import plot_perplex as plotpx
import main as rw
import parameters as p
import matplotlib.pyplot as plt
from useful_and_bespoke import cornertext, colorize, colourbar, colourised_legend
from UM_mass_scaling import um_mass_scaling
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec

""" get earth benchmark """
Tp = 1600
earth = rw.build_planet(M_p=1 * p.M_E, test_CMF=0.325, test_oxides=wt_oxides_MD95,
                        maxIter=30, tol=1e-4, n=1200,
                        Tp=Tp,  # core_efficiency=0.8,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Earth_' + str(Tp) + 'K',
                        )
print('earth mgsi', earth.mgsi)
print('earth mg#', earth.wt_oxides['MgO'] / (earth.wt_oxides['MgO'] + earth.wt_oxides['FeO']))
print('earth core eff', earth.core_eff_from_cmf())

sun = rw.build_planet(M_p=1 * p.M_E,star='sun',
                        maxIter=30, tol=1e-4,# n=800,
                        Tp=Tp,  core_efficiency=0.88,
                        plot_all=False, get_saturation=True, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Sun300_' + str(Tp) + 'K',
                        )
print('sun CMF', sun.CMF)

""" plot variable across masses with composition continuity """
plt.rc('text', usetex=True)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
dirs = [perplex_path_default + 'output/apollo/hypatia' + s + 'M_1600K_88Fe_hires/' for s in ('0,1', '0,3', '0,5', '1', '1,5', '2', '2,5', '3', '4', '5')]
fig, ax = plotpx.compare_pop_fillbetween(dirs, 'M_p', 'mass_um', x_scale=p.M_E ** -1, y_scale=p.M_E ** -1, xlog=False, ylog=False,
                           ylabel='Upper mantle mass ($M_\oplus$)', xlabel='Planet mass ($M_\oplus$)', sigma=2,
                            save=True, show=False, extension='.pdf', xlim=(0.1, 5), ylim=(None, 0.3), #(0.4, 1.8),
                                         show_n=False,
                           labelsize=16, legsize=14, ticksize=14, earth=earth, sun=sun, earth_real=None, #1.06e24,
                                         show_scaling_fn=um_mass_scaling, scalinglabel=r'Constant-$\rho$ scaling',
                                         c='xkcd:peach',dpi=400)


# dirs = [perplex_path_default + 'output/apollo/hypatia' + s + 'M_1600K_99Fe_hires/' for s in ('0,1', '0,3', '0,5', '1', '1,5', '2', '2,5', '3', '4', '5')]
# fig, ax = plotpx.compare_pop_fillbetween(dirs, 'M_p', 'mass_um', x_scale=p.M_E ** -1, y_scale=1e-24, xlog=False, ylog=False,
#                            ylabel='Upper mantle mass ($10^{24}$ kg)', xlabel='Planet mass ($M_\oplus$)', sigma=2,
#                             save=True, show=False, extension='.png', xlim=(0.1, 5), ylim=(0.4, 1.8), show_n=False,
#                            labelsize=16, legsize=14, ticksize=14, earth=earth, sun=sun, earth_real=None, #1.06e24,
#                                          show_scaling_fn=um_mass_scaling, scalinglabel=r'Constant-$\rho$ scaling',
#                                          filename='mass_um_0fe')
#
# plt.show()


"""add different Fe to same axes"""
# dirs = [perplex_path_default + 'output/hypatia' + s + 'M_1900K/'  # 'M_1600K_50Fe
#         for s in ('0,1', '0,3', '0,5', '1', '2', '3', '4', '5')]
# fig, ax = plotpx.compare_pop_fillbetween(dirs, 'M_p', 'mass_um', x_scale=p.M_E ** -1, y_scale=1e-24, xlog=False, ylog=False,
#                            ylabel='Upper mantle mass ($10^{24}$ kg)', xlabel='Planet mass ($M_\oplus$)', sigma=2,
#                             save=True, show=False, extension='.png', xlim=(0.1, 5), show_n=False,
#                            labelsize=12, earth=None, fig=fig, ax=ax)

# dirs = [perplex_path_default + 'output/hypatia' + s + 'M_1600K_99Fe/' for s in ('0,1', '0,3', '0,5', '1', '2', '3', '4', '5')]
# fig, ax = plotpx.compare_pop_fillbetween(dirs, 'M_p', 'mass_um', x_scale=p.M_E ** -1, y_scale=1e-24, xlog=False, ylog=False,
#                            ylabel='Upper mantle mass ($10^{24}$ kg)', xlabel='Planet mass ($M_\oplus$)', sigma=2,
#                             save=True, show=False, extension='.png', xlim=(0.1, 5), show_n=False,
#                            labelsize=12, earth=None, show_scaling_fn=um_mass_scaling, fig=fig, ax=ax)


""" just water mass """
# plotpx.compare_pop_fillbetween(dirs, 'M_p', 'mass_h2o_um', x_scale=p.M_E ** -1, y_scale=p.TO ** -1, xlog=False, ylog=False,
#                            ylabel='UM water capacity (Earth oceans)', xlabel='$M_p$ ($M_\oplus$)', sigma=2,
#                             save=True, show=True, extension='.png',
#                            labelsize=12, earth=earth)



""" scatter plot of all runs """

# # plotpx.compare_pop_scatter(dirs, 'M_p', 'mass_um', x_scale=p.M_E ** -1, y_scale=1, xlog=False, ylog=False,
# #                            ylabel='Upper mantle mass (kg)', xlabel='$M_p$ ($M_\oplus$)', title='all Mg/Si',
# #                            filename='pops_M_p_mass_um', save=True, show=True, extension='.png',
# #                            labelsize=12, lw=0.2, vmin=0.2, vmax=2.5, c='rainbow', earth=earth, v_name='mgsi', alpha=1)

plt.show()
