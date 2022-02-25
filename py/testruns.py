import numpy as np
import perple_x as px
import parameters as p
import eos
import os


""" one-off variations """
Mp = 5
CMF = 0.33
Tp = 1600
wt_oxides = px.wt_oxides_Earth  # px.update_MgSi(1.27, px.wt_oxides_Earth)
px.build_planet(M_p=Mp*p.M_E, test_CMF=CMF, test_oxides=wt_oxides, Tp=Tp,
                maxIter=30, tol=1e-4, n='auto',
                plot=True, comp_stacked=False, get_saturation=False,
                vertex_data='stx11ver', option_file='perplex_option_new', excluded_phases=[],
                name='test_5M_stx11')

""" Dorn+ 2018 Ca-Al rich planet """
# wt_oxides_CaAl = {'SiO2': 26.99, 'MgO': 13.99, 'CaO': 16.0, 'Al2O3': 43.0,  'FeO': 0.01, 'Na2O': 0.01}
# # wt_oxides_CaAl = {'SiO2': 27.0, 'MgO': 14, 'CaO': 16.0, 'Al2O3': 43.0}
# px.build_planet(name='2M_CaAl', M_p=2*p.M_E, test_CMF=0, test_oxides=wt_oxides_CaAl, Tp=1200,
#                 vertex_data='stx11ver', comp_stacked=True, rho_m0=1.15*p.rho_E,
#                 plot=True, get_saturation=False,
#                 maxIter=20, verbose=True, n=1000,
#                 # excluded_phases=['Opx', 'odi', 'en', 'fs', 'ts', 'an', 'coe']
#                 )

# px.build_planet(name='5M_10CMF', M_p=5*p.M_E, test_CMF=0.1, test_oxides=px.wt_oxides_Earth,
#                 plot=True, store_all=True, get_saturation=True, overwrite=True,
#                 maxIter=20)


"""test Hakim EoS"""
# # pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, psurf, Tp)
# _, density, alpha, cp = eos.EOS_all(234.399, 1700, 4)
# print('rho Bouchet', density)
# _, density, alpha, cp = eos.EOS_all(234.4, 1700, 4)
# print('rho Hakim', density)
