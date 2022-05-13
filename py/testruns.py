import numpy as np
import perplexdata as px
import parameters as p
import main as rw
import eos
import os

# tol = 1e-2  # quick
# n = 800  # quick
tol = 1e-4  # better
n = 'auto'  # better


""" one-off variations """
Mp = 1
Tp = 1600

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
star = 'HIP 27384'
for m in [1, 2, 4]:
    dat = rw.build_planet(M_p=m*p.M_E, star=star, Tp=Tp, core_efficiency=0.8,
                    maxIter=30, tol=tol, n=n,  #tol=1e-4, n='auto',
                    plot_all=True, verbose=True, clean=True,
                    vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[], get_saturation=True,
                    output_parent_path=px.perplex_path_default + 'output/hypatia' + str(m) + 'M/',
                    plot_kwargs={'comp_stacked': True, 'p_max': 150}
                    )
    print(m, 'mass_um', dat.mass_um, 'mtl res', len(dat.df_all))

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
