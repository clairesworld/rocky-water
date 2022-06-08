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
perplex_path = '/home/cmg76/Works/perple_x/'  # will have to edit this in .dat after copying over...
n_sample = -1

# set planet parameters
Tp = 1600
core_eff = 0.8831461545794602
n = 1200

# run at higher res over masses (primarily to get upper mantle)
for Mp in [0.1, 0.3, 0.5, 1, 1.5, 2]:
    if isinstance(Mp, float):
        mass_str = str(Mp).replace('.', ',')
    elif isinstance(Mp, int):
        mass_str = str(Mp)
    directory = 'output/hypatia' + mass_str + 'M_' + str(Tp) + 'K_' + str(int(core_eff * 100)) + 'Fe_hires/'
    planet_dict = {'M_p': Mp, 'Tp': Tp, 'core_efficiency': core_eff,
                   'maxIter': 30, 'tol': 1e-4, 'n': n,
                   'get_saturation': True, 'verbose': True,
                   }

    # # get water capacity across planets
    planets = rw.planets_from_hypatia(n_sample,
                                      # stopafter='2MASS 19461589+4406211',
                                      use_local_composition=False,
                                      perplex_path=perplex_path,
                                      output_parent_path=perplex_path + directory,
                                      **planet_dict)

# # or, load from pickle
# planets = rw.read_dir(px.perplex_path_default + 'output/hypatia2M/')

#########################################################################################################################
# already ran the below but leaving it here for ref
# dirs = [px.perplex_path_default + 'output/' + s + '/' for s in ('hypatia1M', 'hypatia2M', 'hypatia4M')]
# for dir_ in dirs:
#     rw.update_dir(dir_, px.PerplexData.get_um_mass, store=True)
#     rw.update_dir(dir_, px.PerplexData.get_mgsi, store=True)
