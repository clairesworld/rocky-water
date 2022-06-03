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

# set planet parameters
Tp = 1600
core_eff = 0.8831461545794602

for Mp in [3]:
    if isinstance(Mp, float):
        mass_str = str(Mp).replace('.', ',')
    elif isinstance(Mp, int):
        mass_str = str(Mp)
    directory = 'output/hypatia' + mass_str + 'M_' + str(Tp) + 'K_' + str(int(core_eff * 100)) + 'Fe_hires/'
    planet_dict = {'M_p': Mp, 'Tp': Tp, 'core_efficiency': core_eff,
                   'maxIter': 30, 'tol': 1e-4,  'n': 1600,
                   'get_saturation': True, 'verbose': True,
                   'output_parent_path': px.perplex_path_default + directory}

    # # get water capacity across planets
    planets = rw.planets_from_hypatia(n_sample,
                                      # stopafter='2MASS 19461589+4406211',
                                      use_local_composition=True, **planet_dict)




# # or, load from pickle
# planets = rw.read_dir(px.perplex_path_default + 'output/hypatia2M/')

#########################################################################################################################
# already ran the below but leaving it here for ref
# dirs = [px.perplex_path_default + 'output/' + s + '/' for s in ('hypatia1M', 'hypatia2M', 'hypatia4M')]
# for dir_ in dirs:
#     rw.update_dir(dir_, px.PerplexData.get_um_mass, store=True)
#     rw.update_dir(dir_, px.PerplexData.get_mgsi, store=True)
