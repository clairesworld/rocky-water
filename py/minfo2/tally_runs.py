import numpy as np
import pandas as pd
import os
import pathlib
import subprocess
# import py.perplexdata as px
# import py.bulk_composition as bulk
import py.main as rw
import py.parameters as p
import perplexfugacitydata as pfug
import meltsfugacitydata as mfug
import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
from py.useful_and_bespoke import iterable_not_string

"""
count number of runs that worked
"""

def completion_df(output_parent_path, T_of_interest=1373.0, col_out=['logfo2'], verbose=False, **kwargs):
    if not iterable_not_string(col_out):
        col_out = [col_out]

    # get directory names in folder
    subfolders = rw.get_run_dirs(output_path=output_parent_path)
    n_stars = len(subfolders)

    df_master = pd.DataFrame(columns=['name', 'p_missing', 'p_GPa'] + col_out, index=range(n_stars))
    if subfolders:  # nonzero
        # loop over runs
        for row, sub in enumerate(subfolders):
            failed = False
            name = os.path.basename(sub)
            df_master.loc[row].name = name
            # print(sub + '/' + name + '_results.csv', os.path.exists(sub + '/' + name + '_results.csv'))
            if (len(os.listdir(sub)) > 0) and os.path.exists(sub + '/' + name + '_results.csv'):  # 1 if contains nH_star e.g.
                df = pd.read_csv(sub + '/' + name + '_results.csv', sep='\t')
                # find T_of_iterest
                if (T_of_interest == df['T(K)']).any() or (T_of_interest + 0.15 == df['T(K)']).any():

                    df_master.loc[row].p_missing = df['P(bar)'].isnull().sum()
                    df_master.loc[row].p_GPa = df['P(bar)'].to_numpy() * 1e-4
                    for c in col_out:
                        df_master.loc[row, c] = df[c].to_numpy()
                else:
                    failed = True

                # if model == 'melts':
                #     dat = mfug.init_from_results(name=name, output_parent_path=output_parent_path, verbose=False)
                # elif model == 'perplex':
                #     dat = pfug.init_from_results(name=name, output_parent_path=output_parent_path)

                # if dat is not None:
                #     # loop over components to check (subplots)
                #     for row, component in enumerate(col_out):
                #
                #         # search for component in oxides
                #         if component in dat.wt_oxides.keys():
                #             x = dat.wt_oxides[component]
                #         # search for component in phase comp
                #         elif 'X_' + component in row.index:
                #             x = row['X_' + component]
                #         else:
                #             if verbose:
                #                 print(name, ':', component, 'not found in', dat.wt_oxides.keys(), 'or', row.index)

                    # else:
                    #     print('problem initialising data object for', name)
            else:
                failed = True
            if failed:
                df_master.loc[row].p_missing = 999
                df_master.loc[row].p_GPa = np.nan
                for c in col_out:
                    df_master.loc[row].c = np.nan
                if verbose:
                    print(sub, 'is empty')

    return df_master

def tally_runs(output_parent_path, **kwargs):
    df = completion_df(output_parent_path, **kwargs)
    good = [df.loc[ii].name for ii in range(len(df)) if df.loc[ii].p_missing == 0]
    bad = [df.loc[ii].name for ii in range(len(df)) if df.loc[ii].name not in good]
    print('good:', len(good), '/', len(df))
    print('missing:', len(bad), '/', len(df))
    return good, bad


X_ferric = 3
coreeff_list = [80, 85, 88, 99]
dirs = [fo2plt.output_parent_mlt_earth + 'hypatia_' + str(ce) + 'coreeff_' + str(X_ferric) + 'ferric_ext/' for ce in coreeff_list]

for d in dirs:
    print('-----\n', d)git
    good, bad = tally_runs(d)