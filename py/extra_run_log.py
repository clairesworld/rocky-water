import numpy as np
import parameters as p
import ask_hypatia as hyp
import perplexdata as px
import pickle as pkl
import plot_perplex as plotpx
import main as rw
from bulk_composition import stellar_mantle
from saturation import TO
import matplotlib.pyplot as plt

""" get every star from hypatia (run once) - todo only with measured mg?"""
# hyp.retrieve_star_names(exo_hosts=False, writeto='all_hypatia_names.txt')
# hyp.retrieve_star_names(exo_hosts=True, writeto='host_names.txt')

""" run over names list """
n_sample = -1


"""run over all stars at 1 M_E and different core partitioning"""
# Tp = 1900
# Mp = 1
# for core_eff in [0.5, 0.6, 0.7, 0.9, 0.999]:
#     if isinstance(Mp, float):
#         mass_str = str(Mp).replace('.', ',')
#     elif isinstance(Mp, int):
#         mass_str = str(Mp)
#     directory = 'output/hypatia' + mass_str + 'M_' + str(Tp) + 'K_' + str(int(core_eff * 100)) + 'Fe/'
#     planet_dict = {'M_p': Mp, 'Tp': Tp, 'core_efficiency': core_eff,
#                    'maxIter': 30, 'tol': 1e-4,  # 'n': 1600,
#                    'get_saturation': True, 'verbose': True,
#                    'output_parent_path': px.perplex_path_default + directory}
#     # # get water capacity across planets
#     planets = rw.planets_from_hypatia(n_sample,
#                                       use_local_composition=True, **planet_dict)
#

"""vary Mg/Si with otherwise Earth compositon and CMF for fig. 4"""
Tp = 1600
# test_CMF = 0.325
core_eff = 0.8831461545794602
for Mp in [0.1, 0.5, 2, 3]:
    directory = 'output/MgSi_from_earth_M' + str(Mp).replace('.', ',') + '/'
    for ms in np.linspace(0.6, 2, num=12, endpoint=True):
    # for ms in [0.7244359600749921, 1.0715193052376069, 1.4125375446227555]:  # or, hypatia 2sigma range
        oxides = rw.update_MgSi(MgSi=ms, oxides=px.wt_oxides_MD95)
        planet_dict = {'M_p': Mp*p.M_E, 'Tp': Tp, #'test_CMF': test_CMF,
                       'core_efficiency': core_eff,
                       'test_oxides': oxides,
                       'maxIter': 30, 'tol': 1e-4, #'n': 1600,
                       'get_saturation': True, 'verbose': True, 'star': None,
                       'output_parent_path': px.perplex_path_default + directory}
        rw.build_planet(plot_all=False, **planet_dict)


"""start with solar composition, vary FeO in mantle via tuning core efficiency (so also change CMF)"""
# Tp = 1600
# Mp = 1
# directory = 'output/MgFe_from_sun/'
# for core_eff in [0.7, 0.75, 0.88, 0.9, 0.95, 0.999]:
#     planet_dict = {'M_p': Mp*p.M_E, 'Tp': Tp, 'core_efficiency': core_eff,
#                    'maxIter': 30, 'tol': 1e-4, #'n': 1600,
#                    'get_saturation': True, 'verbose': True, 'star': 'sun',
#                    'output_parent_path': px.perplex_path_default + directory}
#     pl = rw.build_planet(plot_all=False, **planet_dict)

# start with solar composition and vary Mg/Si for fixed core eff - should be the same as varying from Earth comp
# Tp = 1600
# Mp = 1
# core_eff = 0.8831461545794602
# nH_star = [eval('p.' + el + '_sol') for el in ['si', 'mg', 'ca', 'al', 'fe']]
# oxides = bulk_composition(['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO'], nH_star, core_eff)
# directory = 'output/MgSi_from_sun/'
# for ms in np.linspace(1.7, 3, num=5, endpoint=True):
#     oxides_new = rw.update_MgSi(MgSi=ms, oxides=oxides)
#     planet_dict = {'M_p': Mp*p.M_E, 'Tp': Tp, 'core_efficnecy': core_eff, 'test_oxides': oxides_new,
#                    'maxIter': 30, 'tol': 1e-4,
#                    'get_saturation': True, 'verbose': True, 'star': 'sun',
#                    'output_parent_path': px.perplex_path_default + directory}
#     pl = rw.build_planet(plot_all=False, **planet_dict)

"""keep Earth compositon and change only Mg/Fe ratio in mtl (same core size)"""
# Tp = 1600
# Mp = 1
# cmf = 0.325
# directory = 'output/MgFe_from_earth/'
# for ms in [50, 200]:  #np.linspace(1.4, 7, num=9, endpoint=True):  # 7 is core_eff > 0.9, 1.4 is core_eff ~ 0.5
#     oxides_new = rw.update_MgFe(MgFe=ms, oxides=px.wt_oxides_Earth)
#     planet_dict = {'M_p': Mp*p.M_E, 'Tp': Tp, 'test_CMF': cmf, 'test_oxides': oxides_new,
#                    'maxIter': 30, 'tol': 1e-4,
#                    'get_saturation': True, 'verbose': True, 'star': 'sun',
#                    'output_parent_path': px.perplex_path_default + directory}
#     pl = rw.build_planet(plot_all=False, **planet_dict)


"""vary only stellar iron in hypatia range (solar +- 0.5)"""
# wt_min = {'SiO2': 53.278903580618, 'MgO': 37.423827767238684, 'CaO': 3.210589842392148, 'Al2O3': 4.029027382228369, 'FeO': 2.0576514275227846}
# wt_max = {'SiO2': 44.953944292463035, 'MgO': 31.57626294080103, 'CaO': 2.7089273093327835, 'Al2O3': 3.399481977316648, 'FeO': 17.361383480086502}
# names = ['Fe_min', 'Fe_max']
# Tp = 1600
# core_eff = 0.88
# Mp = 1
# directory = 'output/bulk_Fe_test/'
# for wt, nm in zip([wt_min, wt_max], names):
#     planet_dict = {'M_p': Mp*p.M_E, 'Tp': Tp, 'test_oxides': wt, 'core_efficiency': core_eff,
#                    'maxIter': 30, 'tol': 1e-4,
#                    'get_saturation': True, 'verbose': True, 'star': None, 'name': nm,
#                    'output_parent_path': px.perplex_path_default + directory}
#     pl = rw.build_planet(plot_all=False, **planet_dict)

"""Mars 1600 K"""
# mars_oxides = px.wt_oxides_Mars
# mars_oxides['FeO'] = 13.7  # update Khan+ 2022
# planet_dict = {'M_p': 6.9e23, 'Tp': 1600, 'test_CMF': 0.25, 'test_oxides': mars_oxides,
#                'maxIter': 30, 'tol': 1e-4,
#                'get_saturation': True, 'verbose': True, 'star': 'sun',
#                'output_parent_path': px.perplex_path_default, 'name': 'Mars_1600K'}
# pl = rw.build_planet(plot_all=False, **planet_dict)
# print(pl.core_eff_from_cmf())

""" initial tests of Fe-dependent solidus for 1 ME """
# n_sample = 20
# core_eff = 0.60
# Mp = 1
# directory = 'output/test_compositional_solidus_60Fe/'
# planet_dict = {'M_p': Mp, 'Tp': 'Tm', 'core_efficiency': core_eff, 'profile': 'warm',
#                'maxIter': 30, 'tol': 1e-4, 'get_saturation': True, 'verbose': True,
#                'output_parent_path': px.perplex_path_default + directory}
# planets = rw.planets_from_hypatia(n_sample,
#                                   use_local_composition=True, **planet_dict)

""" trying fake Si partitioning """
# n_sample = -1
# core_eff = 0.88
# Tp = 1600
# Mp = 1
# for x_Si_core in [1, 6, 9]:
#     directory = 'output/test_Si_core' + str(x_Si_core) + '/'
#     planet_dict = {'M_p': Mp, 'Tp': Tp, 'core_efficiency': core_eff, 'x_Si_core': x_Si_core,
#                    'maxIter': 30, 'tol': 1e-4, 'get_saturation': True, 'verbose': True,
#                    'output_parent_path': px.perplex_path_default + directory}
#     planets = rw.planets_from_hypatia(n_sample,
#                                       use_local_composition=True, solve_interior=False, **planet_dict)
#
#     # make hist of mg/si - above runs in like a second
#     plotpx.pop_hist1D(planets, 'mgsi', save=False, show=False, showmedian=True, showsigma=True, bins=100, alpha=0.3)
# plt.show()
#
