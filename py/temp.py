import numpy as np
import parameters as p
import ask_hypatia as hyp
import perplexdata as px
import pickle as pkl
import plot_perplex as plotpx
import main as rw
from saturation import TO

""" run over names list """
# set run parameters
n_sample = -1


# set planet parameters
Tp = 1600
# core_eff = 0.8831461545794602
Mp = 1

for core_eff in [0.98]:
    if isinstance(Mp, float):
        mass_str = str(Mp).replace('.', ',')
    elif isinstance(Mp, int):
        mass_str = str(Mp)
    directory = 'output/hypatia' + mass_str + 'M_' + str(Tp) + 'K_' + str(int(core_eff * 100)) + 'Fe/'
    planet_dict = {'M_p': Mp, 'Tp': Tp, 'core_efficiency': core_eff,
                   'maxIter': 30, 'tol': 1e-4,  # 'n': 1600,
                   'get_saturation': True, 'verbose': True,
                   'output_parent_path': px.perplex_path_default + directory}

    # # get water capacity across planets
    planets = rw.planets_from_hypatia(n_sample,
                                      restart='2MASS 19592683+4549384',
                                      use_local_composition=True, **planet_dict)



# when this is done also do 0.65