import pandas as pd
import numpy as np
import main as rw
import bulk_composition as bulk

ce = 0.88
Xf = 0.03
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO', 'TiO2', 'Na2O', 'Cr2O3']
opp = '/home/claire/Works/perple_x/output/comps/'

planet_kwargs = {'core_efficiency': ce, 'X_ferric': Xf}
pl_list = rw.planets_from_hypatia(n_sample=-1, M_p=1, Tp=999, oxide_list=oxide_list, oxides=oxide_list,
                                  names_file='/home/claire/Works/rocky-water/py/host_names.txt',
                                  use_local_composition=True, existing_dir='comps/',
                                  output_parent_path=opp,
                                  run_vertex=False, solve_interior=False, plot_all=False, skip_stars=[],
                                  **planet_kwargs)

cols = ['name'] + [ox + '(wt%)' for ox in oxide_list]
cols[cols.index('FeO(wt%)')] = 'FeO*(wt%)'
cols.extend(['FeO(wt%)', 'Fe2O3(wt%)', 'Fe_c/Fe_bulk(mol%)', 'Fe3+/Fe(mol%)'])
cols.extend(['log10(' + ox[:2] + '/H)' for ox in oxide_list])
df = pd.DataFrame(columns=cols, index=range(len(pl_list)))
for row, pl in enumerate(pl_list):
    df.name.iloc[row] = pl.star.replace(' ', '')
    try:
        df['FeO*(wt%)'].iloc[row] = pl.wt_oxides['FeO']  # full FeO complement
        # print('FeO*', pl.wt_oxides['FeO'])

        # split up ferric/ferrous
        wt_oxides_split = bulk.insert_fe2o3(pl.wt_oxides, Xf, mass_tot=100, verbose=False)
        for ox in oxide_list:
            df[ox + '(wt%)'].iloc[row] = wt_oxides_split[ox]
        df['Fe2O3(wt%)'].iloc[row] = wt_oxides_split['Fe2O3']
        # print('FeO*', wt_oxides_split['FeO'] + wt_oxides_split['Fe2O3'])

        # add solar
        for ox in oxide_list:
            df['log10(' + ox[:2] + '/H)'].iloc[row] = pl.nH_star[oxide_list.index(ox)]  # log10(N_X/N_H)
    except TypeError as e:
        print('ERROR: pl', pl.star, pl.wt_oxides)

    df['Fe_c/Fe_bulk(mol%)'].iloc[row] = ce * 100
    df['Fe3+/Fe(mol%)'].iloc[row] = Xf * 100

df.to_csv(opp + 'hypatia_compositions.csv', sep=',')


