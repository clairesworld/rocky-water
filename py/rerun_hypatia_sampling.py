import numpy as np
import parameters as p
import ask_hypatia as hyp
import perplexdata as px
import pickle as pkl
import plot_perplex as plotpx
import main as rw
from saturation import TO

""" get every star from hypatia (run once) - todo only with measured mg?"""
# hyp.retrieve_star_names(exo_hosts=False, writeto='all_hypatia_names.txt')
# hyp.retrieve_star_names(exo_hosts=True, writeto='host_names.txt')

""" run over names list """
# set run parameters
n_sample = -1
# directory = 'output/hypatia5M_n1600/'
directory = 'output/hypatia2M/'

# set planet parameters
Mp = 2
Tp = 1600
planet_dict = {'M_p': Mp, 'Tp': Tp, 'core_efficiency': 0.8,
               'maxIter': 30, 'tol': 1e-4, #'n': 1600,
               'get_saturation': True, 'verbose': True,
               'output_parent_path': px.perplex_path_default + directory}

# # get water capacity across planets
planets = rw.planets_from_hypatia(n_sample, **planet_dict)

# # or, load from pickle
# planets = rw.read_dir(px.perplex_path_default + 'output/hypatia2M/')

#########################################################################################################################
# already ran the below but leaving it here for ref
# dirs = [px.perplex_path_default + 'output/' + s + '/' for s in ('hypatia1M', 'hypatia2M', 'hypatia4M')]
# for dir_ in dirs:
#     rw.update_dir(dir_, px.PerplexData.get_um_mass, store=True)
#     rw.update_dir(dir_, px.PerplexData.get_mgsi, store=True)