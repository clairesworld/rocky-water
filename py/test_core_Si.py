import bulk_composition as bulk
import numpy as np
import parameters as p
import main as rw
import pickle as pkl
import time
import perplexdata as px
import ask_hypatia as hyp
import pandas as pd
import os
import matplotlib.pyplot as plt

perplex_path = '/home/g/guimond/Work/perple_x/'
oxide_list = ['MgO', 'SiO2', 'CaO', 'Al2O3', 'FeO']  # 'Na2O

def core_eff_from_cmf(CMF, M_p, wt_oxides, core_light_wt=None):  # todo
    M_c = CMF * M_p  # core mass in kg
    M_m = M_p - M_c  # mantle mass in kg
    n_Fe_mtl = wt_oxides['FeO'] / 100 * M_m / p.M_FeO  # n_Fe = n_FeO
    n_Fe_core = M_c / p.M_Fe  # pure Fe core
    core_eff = n_Fe_core / (n_Fe_core + n_Fe_mtl)
    return core_eff


# test for sol

nH_star = [p.mg_sol, p.si_sol, p.ca_sol, p.al_sol, p.fe_sol]
core_eff = 0.85
core_Si_wtpt = 5  # need 10 to match cmf and mg/si at 0.9 Fe_eff, but something like 4-5 max is favoured for earth (hirose)

wt_oxides = bulk.stellar_mantle_corelightelements(oxide_list, nH_star, core_eff, core_Si_wtpt=core_Si_wtpt)
cmf = bulk.core_mass_fraction(wt_oxides, core_eff, core_light_wtpt={'Si': core_Si_wtpt, 'Ni': 85 - core_Si_wtpt})  # Ni is everything else (from Hirose+ review - sum to 12.5%)
#
# print('\n\nno core Si\n')
# wt_oxides = bulk.stellar_mantle_corelightelements(oxide_list, nH_star, core_eff, core_Si_wtpt=0)
# cmf = bulk.core_mass_fraction(wt_oxides, core_eff, core_light_wtpt=None)
# cmf2 = bulk.Fecore_mass_fraction(wt_oxides, core_eff)
# print('old cmf function', cmf2)


def count_qtz_planets(output_path=''):
    import matplotlib.pyplot as plt
    n_qtz = 0
    dats = rw.read_dir(output_path, subsample=None, verbose=True)
    mgsi = []
    for dat in dats:
        # print(dat.df_comp.head())
        try:
            x_coes = dat.df_comp['coes']
            n_qtz += 1
            mgsi.append(dat.mgsi)
        except KeyError:
            try:
                x_qtz = dat.df_comp['qtz']
                n_qtz += 1
                mgsi.append(dat.mgsi)
            except KeyError:
                continue
    print(n_qtz, '/', len(dats), '=', n_qtz/len(dats)*100, '% of planets stabilise pure SiO2')
    print('    with max Mg/Si =', np.max(mgsi))

    plt.hist(mgsi, bins=20)
    plt.xlabel('Mg/Si')
    plt.show()
    return n_qtz


def count_pyrolite_planets(output_path=''):
    n_pyr = 0
    dats = rw.read_dir(output_path, subsample=None, verbose=True)
    mgsi_olfree = []
    mgsi_opxfree = []
    femg = []
    wt_FeO = []
    for dat in dats:
        # print(dat.df_comp.head())
        if ('O' in dat.df_comp.columns) and ('Opx' in dat.df_comp.columns):
            n_pyr += 1
        elif ('Opx' in dat.df_comp.columns):  # ol-free
            mgsi_olfree.append(dat.mgsi)
        elif ('O' in dat.df_comp.columns):  # opx-free
            mgsi_opxfree.append(dat.mgsi)
            print(dat.df_comp.head())
        femg.append(bulk.mass_ratio_to_mol_ratio(dat.wt_oxides['FeO'], dat.wt_oxides['MgO'], s_i='FeO', s_j='MgO'))
        wt_FeO.append(dat.wt_oxides['FeO'])
    print('\n', n_pyr, '/', len(dats), '= {:.1f}% of planets stabilise ol+opx'.format(n_pyr/len(dats)*100))

    print('\nWhen do ol and opx disappear?')
    print(' opx-free limits', np.quantile(mgsi_opxfree, [0.16, 0.84] ))  # 2 sigma [0.0227, 0.9773]
    print(' ol-free limits', np.quantile(mgsi_olfree, [0.16, 0.84]))
    print(' min Mg/Si for opx-free:', np.min(mgsi_opxfree))
    print(' max Mg/Si for ol-free:', np.max(mgsi_olfree))
    print(' average Fe/Mg = {:.3f}'.format(np.mean(femg)), 'average mantle FeO = {:.2f} wt%'.format(np.mean(wt_FeO)))
    fig, axes = plt.subplots(1, 2)
    for set, axnum, label in zip([mgsi_olfree, mgsi_opxfree], [0,1], ['Mg/Si, ol-free', 'Mg/Si, opx-free']):
        axes[axnum].hist(set, histtype='step', bins=50)
        axes[axnum].set_xlabel(label)
        axes[axnum].axvline(np.median(set), c='k')
        axes[axnum].axvline(np.mean(set), c='k', ls='--')
        axes[axnum].set_title(os.path.basename(output_path) + ' n = ' + str(len(set)))
    return n_pyr


# count_qtz_planets(output_path=perplex_path + 'output/coreSi_0/')

""" when does ol or opx disappear? """

# count_pyrolite_planets(output_path=perplex_path + 'output/coreSi0_coreeff0,5/')
count_pyrolite_planets(output_path=perplex_path + 'output/coreSi0_coreeff0,88/')
# count_pyrolite_planets(output_path=perplex_path + 'output/coreSi0_coreeff0,999/')
# core_eff = 0.1 --> get lots of pure FeO phases at 4 GPa, perplex is probably breaking


""" average Mg/Si increase over bulk """

def mgsi_increase_from_core(core_eff=0.88, core_Si_wtpt=5, oxide_list=oxide_list,
                              path_to_tsv='/home/g/guimond/Work/hypatia-compositions/hypatia-04122023.tsv'):
    # read stellar abundances
    df = pd.read_csv(path_to_tsv, sep='\t', header=0)

    # remove rows with no measurements
    df.replace('', np.nan, inplace=True)  # first convert empty string to nan
    df.dropna(subset=[ox[:2] for ox in oxide_list], inplace=True)

    mgsi_increase_all = []
    mgsi_bulk_all = []
    mgsi_partitioned_all = []
    for ii, star in enumerate([name.rstrip() for name in df.Name.to_numpy()]):
        nH_star = []
        for metal in oxide_list:
            el = metal[:2]
            sol_val = eval('p.' + el.lower() + '_sol')
            nH_star.append(df.iloc[ii][el] + sol_val)

        # get bulk mantle oxides
        wt_oxides_purecore = bulk.stellar_mantle_corelightelements(oxide_list=oxide_list, nH_star=nH_star, core_eff=core_eff,
                                                     core_Si_wtpt=0)
        mgsi_bulk = (wt_oxides_purecore['MgO'] / p.M_MgO) / (wt_oxides_purecore['SiO2'] / p.M_SiO2)

        wt_oxides_partitioned = bulk.stellar_mantle_corelightelements(oxide_list=oxide_list, nH_star=nH_star, core_eff=core_eff,
                                                     core_Si_wtpt=core_Si_wtpt)
        mgsi_partitioned = (wt_oxides_partitioned['MgO'] / p.M_MgO) / (wt_oxides_partitioned['SiO2'] / p.M_SiO2)

        mgsi_increase_all.append(mgsi_partitioned - mgsi_bulk)
        mgsi_bulk_all.append(mgsi_bulk)
        mgsi_partitioned_all.append(mgsi_partitioned)

    print('\n\n\n------------------')
    print('average Mg/Si difference', np.mean(mgsi_increase_all))
    print('bulk average', np.mean(mgsi_bulk_all))
    print('partitioned average', np.mean(mgsi_partitioned_all))

# mgsi_increase_from_core(core_eff=0.5, core_Si_wtpt=5)

plt.show()
