import perplexdata as px
import main as rw
# from saturation import TO


""" run over names list """
# set run parameters
perplex_path = '/home/cmg76/Works/perple_x/'  # will have to edit this in .dat after copying over...
# perplex_path = px.perplex_path_default
n_sample = -1

# set planet parameters
Tp = 1900
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
                                      use_local_composition=True,
                                      perplex_path=perplex_path,
                                      output_parent_path=perplex_path + directory,
                                      **planet_dict)
