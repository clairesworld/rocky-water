# import numpy as np
# import parameters as p
# import ask_hypatia as hyp
# from perplexdata import perplex_path_default
# import pickle as pkl
# import plot_perplex as plotpx
import main as rw
# from saturation import TO

""" get every star from hypatia (run once) - todo only with measured mg?"""
# hyp.retrieve_star_names(exo_hosts=False, writeto='all_hypatia_names.txt')
# hyp.retrieve_star_names(exo_hosts=True, writeto='host_names.txt')

""" run over names list """
def run_all_masses(masses=None, Tp=1600, core_eff=0.8831461545794602, n_sample=-1, n='auto', restart=None,
                   perplex_path='/raid1/cmg76/perple_x/', tail='_hires'  # apollo hires defaults
                   ):
    # run at higher res over masses (primarily to get upper mantle)

    if masses is None:
        masses = [0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5]
    for Mp in masses:
        if isinstance(Mp, float):
            mass_str = str(Mp).replace('.', ',')
        elif isinstance(Mp, int):
            mass_str = str(Mp)
        directory = 'output/hypatia' + mass_str + 'M_' + str(Tp) + 'K_' + str(int(core_eff * 100)) + 'Fe' + tail + '/'
        planet_dict = {'M_p': Mp, 'Tp': Tp, 'core_efficiency': core_eff,
                       'maxIter': 30, 'tol': 1e-4, 'n': n,
                       'get_saturation': True, 'verbose': True,
                       }

        # get water capacity across planets
        planets = rw.planets_from_hypatia(n_sample,
                                          restart=restart,
                                          use_local_composition=True,
                                          perplex_path=perplex_path,
                                          output_parent_path=perplex_path + directory,
                                          **planet_dict)
        restart = None  # reset

def run_Fe_partitioning(core_effs=None, Tp=1600, Mp=1, n_sample=-1, n='auto', restart=None,
                   perplex_path='/raid1/cmg76/perple_x/', tail='_hires'  # apollo hires defaults
                   ):
    # run at higher res over masses (primarily to get upper mantle)

    if core_effs is None:
        core_effs = [0.65, 0.75, 0.85, 0.9, 0.93, 0.97, 0.98]
    for core_eff in core_effs:
        if isinstance(Mp, float):
            mass_str = str(Mp).replace('.', ',')
        elif isinstance(Mp, int):
            mass_str = str(Mp)
        directory = 'output/hypatia' + mass_str + 'M_' + str(Tp) + 'K_' + str(int(core_eff * 100)) + 'Fe' + tail + '/'
        planet_dict = {'M_p': Mp, 'Tp': Tp, 'core_efficiency': core_eff,
                       'maxIter': 30, 'tol': 1e-4, 'n': n,
                       'get_saturation': True, 'verbose': True,
                       }

        # get water capacity across planets
        planets = rw.planets_from_hypatia(n_sample,
                                          restart=restart,
                                          use_local_composition=True,
                                          perplex_path=perplex_path,
                                          output_parent_path=perplex_path + directory,
                                          **planet_dict)
        restart = None  # reset


# set run parameters
# perplex_path = '/raid1/cmg76/perple_x/'  # will have to edit this in .dat after copying over...
# perplex_path = perplex_path_default
# tail = '_hires'

# set planet parameters
# Tp = 1600
# core_eff = 0.8831461545794602
# n = 'auto'  # 1200

# run_all_masses(#masses=None,
#                Tp=Tp, core_eff=core_eff,
#             )

# # or, load from pickle
# planets = rw.read_dir(px.perplex_path_default + 'output/hypatia2M/')

#########################################################################################################################
# already ran the below but leaving it here for ref
# dirs = [px.perplex_path_default + 'output/' + s + '/' for s in ('hypatia1M', 'hypatia2M', 'hypatia4M')]
# for dir_ in dirs:
#     rw.update_dir(dir_, px.PerplexData.get_um_mass, store=True)
#     rw.update_dir(dir_, px.PerplexData.get_mgsi, store=True)
