import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
import numpy as np
import meltsfugacitydata as mfug
import perplexfugacitydata as pfug
import py.bulk_composition as bulk
import py.main as rw
import os

output_parent_path_px = '/home/claire/Works/min-fo2/perplex_output/hypatia/hypatia_88coreeff_3ferric_ext/'
output_parent_path_mlts = '/home/claire/Works/min-fo2/alphamelts_output/earth-tea23/hypatia_88coreeff_3ferric_ext/'
exclude_names = []


def get_fe2o3(model, ph, output_parent_path):
    # get directory names in folder
    subfolders = rw.get_run_dirs(output_path=output_parent_path)
    x = []
    n_plag = 0
    if subfolders:  # nonzero
        for ii, sub in enumerate(subfolders):
            name = os.path.basename(sub)
            # print(sub + name + '_results.csv')
            # print(os.path.exists(sub + name + '_results.csv'))

            if (len(os.listdir(sub)) > 1) and (
                    os.path.exists(sub + '/' + name + '_results1373.csv')):  # 1 if contains nH_star e.g.
                if name not in exclude_names:
                    if model == 'perplex':
                        p_of_interest = 1
                        T_of_interest = 1373
                        dat = pfug.init_from_results(name, output_parent_path=output_parent_path)
                    elif model == 'melts':
                        p_of_interest = 10000
                        T_of_interest = 1373.15
                        dat = mfug.init_from_results(name, output_parent_path=output_parent_path, verbose=False,
                                                     )
                    if dat is not None:

                        # load fe2o3 composition
                        d = dat.get_phase_composition_dict(p_of_interest=p_of_interest,
                                                           T_of_interest=T_of_interest,
                                                           component='Fe2O3',
                                                           phases=[ph],
                                                           to_absolute_abundance=False,
                                                           verbose=False)
                        if d is not None:
                            fe2o3_in_phase = d[ph]
                            x.append(fe2o3_in_phase)

                        # plag?
                        if 'X_Plag' in dat.data.columns:
                            n_plag += 1
    print('n plag', n_plag)
    return x


phase = 'Opx'

x_px = get_fe2o3(model='perplex', ph=phase, output_parent_path=output_parent_path_px)
print(x_px)

x_mlt = get_fe2o3(model='melts', ph=phase, output_parent_path=output_parent_path_mlts)
print(x_mlt)

fig = plt.figure()
plt.gca().hist(x_px, bins=50)
plt.gca().set_xlabel('wt.% Fe2O3 in ' + phase +', perplex')

fig = plt.figure()
plt.gca().hist(x_mlt, bins=50)
plt.gca().set_xlabel('wt.% Fe2O3 in ' + phase +', melts')

plt.show()
