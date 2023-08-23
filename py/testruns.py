import numpy as np
import perplexdata as px
import parameters as p
import main as rw
import eos
import os
import matplotlib.pyplot as plt
import plot_perplex as plotpx
import matplotlib.colors

# tol = 1e-2  # quick
# n = 800  # quick
tol = 1e-4  # better
n = 'auto'  # better

# get earth-scaled density
perplex_path = '/home/claire/Works/perple_x/'
for ii, Mp in enumerate([1.31, 1.37, 1.04, 0.69, 1.32, 4.91, 4.36, 0.39, 6.38, 4.13]):
    if ii != 8:
        pl = rw.build_planet(M_p=Mp * p.M_E, name=str(ii), test_CMF=0.325, Tp=1600, star='sun',  tol=1e-4, n='auto',
                             plot_all=False, get_saturation=False, verbose=True, clean=True,
                             perplex_path=perplex_path,
                             output_parent_path=perplex_path + 'output/random_coreeff/',
                             )
        print('earth-scaled bulk density', pl.M_p / (np.pi * 4/3 * pl.R_p**3))



""" one-off variations """
Mp = 1
Tp = 1600

# """ test random stars """
# pl = rw.build_planet(M_p=Mp * p.M_E, core_efficiency=0.88, Tp=Tp, star='HIP 21850',  # tol=1e-4, n='auto',
#                       plot_all=False, get_saturation=True, verbose=True, clean=False,
#                       )


wt_oxides_MD95 = {'SiO2': 45.0, 'MgO': 37.8, 'CaO': 3.55, 'Al2O3': 4.45, 'FeO': 8.05}  # McDonough & Sun 1995 BSE
# McD S case mol fractions
# test_oxides = {'MgO':40.60822798 ,
# 'SiO2': 43.24067322 ,
# 'FeO': 8.05,
# 'CaO': 3.59487446 ,
# 'Al2O3': 4.50622434}  # Mg/Si = 1.4
test_oxides = {'MgO': 38.27887621 , 'SiO2': 45.57002469 ,
'FeO': 8.05,
'CaO': 3.5948746,
'Al2O3': 4.5062245}  # BSE







#test_oxides = rw.update_MgSi(1.4, px.wt_oxides_MD95)
# pl = rw.build_planet(M_p=Mp * p.M_E, core_efficiency=0.88, Tp=Tp,
#                       test_oxides=test_oxides,
#                       maxIter=30, tol=tol, n=1600,  # tol=1e-4, n='auto',
#                       plot_all=False, get_saturation=True, verbose=True, clean=True,
#                       name='MgSi_1,1_stx11',
#                       option_file='perplex_option_claire_mol', vertex_data='stx11ver', excluded_phases=['q']
#                       )

# test stellar iron effect on mineralogy
# oxides = ['MgO', 'SiO2', 'CaO', 'Al2O3', 'FeO']

""" fixed core eff, vary stellar Fe"""
# pl0 = rw.build_planet(M_p=Mp * p.M_E, core_efficiency=0.88, Tp=Tp,
#                       test_nH_star=[p.mg_sol, p.si_sol, p.ca_sol, p.al_sol, p.fe_sol],
#                       maxIter=30, tol=tol, n=n,  # tol=1e-4, n='auto',
#                       plot_all=False, get_saturation=True, verbose=True, clean=True,
#                       name='solar_Fe',
#                       plot_kwargs={'comp_stacked': True}  # , 'plot_phases': True}
#                       )
# pl1 = rw.build_planet(M_p=Mp * p.M_E, core_efficiency=0.88, Tp=Tp,
#                       test_nH_star=[p.mg_sol, p.si_sol, p.ca_sol, p.al_sol, p.fe_sol + 0.2],
#                       maxIter=30, tol=tol, n=n,  # tol=1e-4, n='auto',
#                       plot_all=False, get_saturation=True, verbose=True, clean=True,
#                       name='solar_Fe_plus',
#                       plot_kwargs={'comp_stacked': True}  # , 'plot_phases': True}
#                       )
# pl2 = rw.build_planet(M_p=Mp * p.M_E, core_efficiency=0.88, Tp=Tp,
#                       test_nH_star=[p.mg_sol, p.si_sol, p.ca_sol, p.al_sol, p.fe_sol - 0.2],
#                       maxIter=30, tol=tol, n=n,  # tol=1e-4, n='auto',
#                       plot_all=False, get_saturation=True, verbose=True, clean=True,
#                       name='solar_Fe_minus',
#                       plot_kwargs={'comp_stacked': True}  # , 'plot_phases': True}
#                       )


""" fixed CMF, vary stellar Fe"""
# test_CMF = 0.325
# # pl0 = rw.build_planet(M_p=Mp * p.M_E, test_CMF=test_CMF, Tp=Tp,
# #                       test_nH_star=[p.mg_sol, p.si_sol, p.ca_sol, p.al_sol, p.fe_sol],
# #                       maxIter=30, tol=tol, n=n,  # tol=1e-4, n='auto',
# #                       plot_all=False, get_saturation=True, verbose=True, clean=True,
# #                       name='solar_Fe',
# #                       plot_kwargs={'comp_stacked': True}  # , 'plot_phases': True}
# #                       )
# # pl1 = rw.build_planet(M_p=Mp * p.M_E, test_CMF=test_CMF, Tp=Tp,
# #                       test_nH_star=[p.mg_sol, p.si_sol, p.ca_sol, p.al_sol, p.fe_sol + 0.2],
# #                       maxIter=30, tol=tol, n=n,  # tol=1e-4, n='auto',
# #                       plot_all=False, get_saturation=True, verbose=True, clean=True,
# #                       name='solar_Fe_plus',
# #                       plot_kwargs={'comp_stacked': True}  # , 'plot_phases': True}
# #                       )
# # pl2 = rw.build_planet(M_p=Mp * p.M_E, test_CMF=test_CMF, Tp=Tp,
# #                       test_nH_star=[p.mg_sol, p.si_sol, p.ca_sol, p.al_sol, p.fe_sol - 0.2],
# #                       maxIter=30, tol=tol, n=n,  # tol=1e-4, n='auto',
# #                       plot_all=False, get_saturation=True, verbose=True, clean=True,
# #                       name='solar_Fe_minus',
# #                       plot_kwargs={'comp_stacked': True}  # , 'plot_phases': True}
# #                       )


""" fixed stellar Fe, vary CMF"""
# test_nH_star=[p.mg_sol, p.si_sol, p.ca_sol, p.al_sol, p.fe_sol]
# pl0 = rw.build_planet(M_p=Mp * p.M_E, test_CMF=0.2, Tp=Tp,
#                       test_nH_star=test_nH_star,
#                       maxIter=30, tol=tol, n=n,  # tol=1e-4, n='auto',
#                       plot_all=False, get_saturation=True, verbose=True, clean=True,
#                       name='20',
#                       plot_kwargs={'comp_stacked': True}  # , 'plot_phases': True}
#                       )
# pl1 = rw.build_planet(M_p=Mp * p.M_E, test_CMF=0.3, Tp=Tp,
#                       test_nH_star=test_nH_star,
#                       maxIter=30, tol=tol, n=n,  # tol=1e-4, n='auto',
#                       plot_all=False, get_saturation=True, verbose=True, clean=True,
#                       name='30',
#                       plot_kwargs={'comp_stacked': True}  # , 'plot_phases': True}
#                       )
# pl2 = rw.build_planet(M_p=Mp * p.M_E, test_CMF=0.4, Tp=Tp,
#                       test_nH_star=test_nH_star,
#                       maxIter=30, tol=tol, n=n,  # tol=1e-4, n='auto',
#                       plot_all=False, get_saturation=True, verbose=True, clean=True,
#                       name='40',
#                       plot_kwargs={'comp_stacked': True}  # , 'plot_phases': True}
#                       )
#
# from plot_mgsi_modality import demo_composition_gridspec
#
# sample = [rw.read_name(output_path=px.output_parent_px + 'test_bulk_Fe_coreeff/',
#                        name='solar_Fe' + m) for m in ['_minus', '', '_plus']]
# # sample = [rw.read_name(output_path=px.output_parent_px + 'test_CMF/',
# #                        name=m) for m in ['20', '30', '40']]
# for dat in sample:
#     print('\n', dat.name)
#     print('core eff', dat.core_eff_from_cmf())
#     print('mantle mass', dat.M_p * (1 - dat.CMF))
#     print('total water (kg):', dat.mass_h2o_total)
#     print('total water (% mantle)', dat.mass_h2o_total / (dat.M_p * (1 - dat.CMF)) * 100)
#     # plotpx.single_structure(dat, fig_path=plotpx.fig_path)
# demo_composition_gridspec(dats=sample, comp_var=None, labelsize=12, legsize=10, p_max=150,  # cmap=None,
#                           comp_labels=['low stellar Fe', 'solar', 'high stellar Fe'],
#                           save=True, figtitle='phases_vs_bulkFe_coreeff', log_profile=False, cumulative=True,
#                           phases_order=['gt', 'cpx', 'opx', 'hpcpx', 'ol', 'wad', 'ring', 'pv', 'qtz', 'coes',
#                                         'st', 'wus', 'dvm', 'ppv'], profile_var='mass_h2o(kg)', ylim=(0, 3e21),
#                           prof_ax_label='Cumulative\nmass ' + r'H$_2$O (kg)', show_cum_mass=True, profile_scale=1,
#                           cmap=matplotlib.colors.ListedColormap([plt.cm.tab20b(14),  # gt
#                                                                  plt.cm.tab20b(10),  # cpx
#                                                                  plt.cm.tab20b(11),  # opx
#                                                                  plt.cm.tab20b(9),  # hpcpx
#                                                                  plt.cm.tab20b(7),  # ol
#                                                                  plt.cm.tab20b(5),  # wad
#                                                                  plt.cm.tab20b(4),  # ring
#                                                                  plt.cm.tab20c(19),  # pv
#                                                                  plt.cm.tab20b(19),  # qtz
#                                                                  plt.cm.tab20b(18),  # coes
#                                                                  plt.cm.tab20b(17),  # st
#                                                                  plt.cm.tab20b(3),  # wus
#                                                                  plt.cm.tab20b(0),  # capv
#                                                                  plt.cm.tab20b(1),  # ppv
#                                                                  ])
#                           )
# plt.show()

# CMF = 0.3
# wt_oxides = px.wt_oxides_Earth  # px.update_MgSi(1.27, px.wt_oxides_Earth)
# rw.build_planet(M_p=Mp*p.M_E, test_CMF=CMF, test_oxides=wt_oxides, Tp=Tp,
#                 maxIter=30, tol=tol, n=n, #tol=1e-4, n='auto',
#                 plot_all=True, get_saturation=False, verbose=True, clean=True,
#                 vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                 name='res_test_' + str(n),
#                 plot_kwargs={'comp_stacked': True, 'p_max': 50} #, 'plot_phases': True}
#                 )

# """ test some extremes of Mg/Si from Hypatia """
# rw.build_planet(M_p=Mp*p.M_E, star='HIP 102042', Tp=Tp, core_efficiency=0.8,
#                 maxIter=30, tol=tol, n=n, #tol=1e-4, n='auto',
#                 plot_all=True, verbose=True, clean=True,
#                 vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[], get_saturation=True,
#                 # name='3M_1600_stx11',
#                 plot_kwargs={'comp_stacked': True, 'p_max': 150} #, 'plot_phases': True}
#                 )

# rw.build_planet(M_p=Mp*p.M_E, star='HIP 91300', Tp=Tp, core_efficiency=0.8,
#                 maxIter=30, tol=tol, n=n, #tol=1e-4, n='auto',
#                 plot_all=True, get_saturation=False, verbose=True, clean=True,
#                 vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                 # name='3M_1600_stx11',
#                 plot_kwargs={'comp_stacked': True} #, 'plot_phases': True}
#                 )

""" check out some stars """
# star = 'HIP 27384'
# for m in [1, 2, 4]:
#     dat = rw.build_planet(M_p=m*p.M_E, star=star, Tp=Tp, core_efficiency=0.8,
#                     maxIter=30, tol=tol, n=n,  #tol=1e-4, n='auto',
#                     plot_all=True, verbose=True, clean=True,star = 'HIP 27384'
# for m in [1, 2, 4]:
#     dat = rw.build_planet(M_p=m*p.M_E, star=star, Tp=Tp, core_efficiency=0.8,
#                     maxIter=30, tol=tol, n=n,  #tol=1e-4, n='auto',
#                     plot_all=True, verbose=True, clean=True,
#                     vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[], get_saturation=True,
#                     output_parent_path=px.perplex_path_default + 'output/hypatia' + str(m) + 'M/',
#                     plot_kwargs={'comp_stacked': True, 'p_max': 150}
#                     )
#     print(m, 'mass_um', dat.mass_um, 'mtl res', len(dat.data))
#                     vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[], get_saturation=True,
#                     output_parent_path=px.perplex_path_default + 'output/hypatia' + str(m) + 'M/',
#                     plot_kwargs={'comp_stacked': True, 'p_max': 150}
#                     )
#     print(m, 'mass_um', dat.mass_um, 'mtl res', len(dat.data))

# """ Dorn+ 2018 Ca-Al rich planet """
# wt_oxides_CaAl = {'SiO2': 27, 'MgO': 14, 'CaO': 16.0, 'Al2O3': 43.0,  'FeO': 0.01, 'Na2O': 0.01}
# rw.build_planet(name='3M_CaAl_stx11_1800', M_p=3*p.M_E, test_CMF=0, test_oxides=wt_oxides_CaAl, Tp=1800,
#                 vertex_data='stx11ver', option_file='perplex_option_claire', excluded_phases=['q'],
#                 rho_m0=1.15*p.rho_E, clean=True,
#                 plot_all=True, get_saturation=False,
#                 maxIter=20, verbose=True, n=n, tol=tol, #n=1000, tol=1e-4,
#                 plot_kwargs={'comp_stacked': True}
#                 # plot_phases: ['Aki', 'O', 'Gt', 'Ring', 'Wus', 'Pv', 'CF', 'ca-ov', 'st', 'Ppv', 'Wus', 'seif',]
#                 )1

""" like Umemoto """
# CMF = 0.3
# Tp = 1600
# wt_oxides = {'SiO2': 26.99, 'MgO': 13.99, 'CaO': 16.0, 'Al2O3': 43.0,  'FeO': 1e-5, 'Na2O': 0.01}
# wt_oxides = rw.update_MgSi(1.3, wt_oxides)
# rw.build_planet(M_p=5*p.M_E, test_CMF=CMF, test_oxides=wt_oxides, Tp=Tp,
#                 maxIter=30, tol=tol, n=n,
#                 plot_all=True, comp_stacked=False, get_saturation=False,
#                 vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                 )

""" like Unterborn and Panero """
# masses = [1, 2, 3, 4]
# CMF = 0.25
# Tp = 1600
# wt_oxides = {'SiO2': 50, 'MgO': 50}
# for M in masses:
#     rw.build_planet(M_p=M*p.M_E, test_CMF=CMF, oxides=['SiO2', 'MgO'], test_oxides=wt_oxides, Tp=Tp,
#                     maxIter=30, tol=tol, n=n,
#                     plot_all=True, get_saturation=False,
#                     vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                     name='unterborn_' + str(M) + 'M'
#                     )


"""test Hakim EoS"""
# # pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, psurf, Tp)
# _, density, alpha, cp = eos.EOS_all(234.399, 1700, 4)
# print('rho Bouchet', density)
# _, density, alpha, cp = eos.EOS_all(234.4, 1700, 4)
# print('rho Hakim', density)
