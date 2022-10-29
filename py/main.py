import numpy as np
from py.parameters import M_E, M_Fe, M_FeO, M_MgO, M_SiO2, G, R_E, rho_E
import os
import py.saturation as sat
import time
import py.perplexdata as px
import py.ask_hypatia as hyp
import py.plot_perplex as plotpx
import pickle as pkl
import random


def rename_phases(phases_px):
    """ perplex phases may be different than what you need for plotting """
    new_phases = []
    for ph in phases_px:
        if ph == 'O':
            ph1 = 'ol'
        elif ph == 'C2/c':
            ph1 = 'hpcpx'
        elif ph == 'ca-pv':
            ph1 = 'dvm'
        else:
            ph1 = ph.lower()
        new_phases.append(ph1)
    return new_phases


def update_MgSi(MgSi=None, oxides=px.wt_oxides_Earth):
    """ given an oxide composition in wt%, change MgO and SiO2 to desired molar Mg/Si ratio """
    # assume 100 g
    oxides = {k: v * 100.0/sum(oxides.values()) for k, v in oxides.items()}

    # update m_MgO and m_SiO2, preserving their total mass
    print('assume 100 g')
    print('m total =', sum(oxides.values()), 'g')
    m_tot = oxides['MgO'] + oxides['SiO2']  # total mass (%) of MgO and SiO2 - conserve this
    print('m_MgO + m_SiO2 =', m_tot, 'g')

    # convert to moles
    n_MgO = oxides['MgO'] / M_MgO
    n_SiO2 = oxides['SiO2'] / M_SiO2
    MgSi_old = n_MgO / n_SiO2  # n_Mg = n_MgO and n_Si = n_SiO2
    if MgSi is None:
        MgSi = MgSi_old  # no update
    print('n total =', n_MgO + n_SiO2)
    print('Mg/Si', MgSi_old)

    m_ratio_new = MgSi * (M_MgO / M_SiO2)  # m = nM | new mass ratio
    x = m_ratio_new
    y = m_tot  # conserved
    m_MgO_new = x * y / (1 + x)
    m_SiO2_new = y - m_MgO_new  # m_tot = m_MgO + m_SiO2
    oxides['MgO'] = m_MgO_new
    oxides['SiO2'] = m_SiO2_new

    print('-------\nCHECK:')
    m_tot = oxides['MgO'] + oxides['SiO2']  # total mass (%) of MgO and SiO2 - conserve this
    print('m_MgO + m_SiO2 =', m_tot, 'g')
    n_MgO = oxides['MgO'] / M_MgO  # convert to moles
    n_SiO2 = oxides['SiO2'] / M_SiO2
    MgSi_old = n_MgO / n_SiO2  # n_Mg = n_MgO etc
    print('n total =', n_MgO + n_SiO2)
    print('Mg/Si', MgSi_old, '\n\n\n\n\n\n\n\n')

    return oxides

# update_MgSi(0.7)






def update_MgFe(MgFe=None, oxides=px.wt_oxides_Earth):
    # vary Mg/Fe, use original otherwise
    # update m_MgO and m_FeO, preserving total mass
    m_tot = oxides['MgO'] + oxides['FeO']  # total mass (%) of MgO and FeO - conserve this
    n_MgO = oxides['MgO'] / M_MgO  # convert to moles
    n_FeO = oxides['FeO'] / M_FeO
    n_MgFe_old = n_MgO / n_FeO  # n_Mg = n_MgO etc
    if MgFe is None:
        MgFe = n_MgFe_old  # no update

    m_ratio_new = MgFe * (M_MgO / M_FeO)  # m = nM, in % | new mass ratio
    x = m_ratio_new
    y = m_tot
    m_MgO_new = x * y / (1 + x)
    m_FeO_new = m_tot - m_MgO_new
    oxides['MgO'] = m_MgO_new
    oxides['FeO'] = m_FeO_new
    return oxides


def get_name(M_p=None, star=None, core_efficiency=None, Tp=None, test_CMF=None, test_oxides=None, suffix=None,
             **kwargs):
    mass_str = str(int(np.round(M_p / M_E))) + 'M'
    if test_CMF is not None:
        cmf_str = str(int(test_CMF * 100)) + 'CMF'
    else:
        cmf_str = str(int(core_efficiency * 100)) + 'Ceff'
    if test_oxides is not None:
        comp_str = str(int(np.round(test_oxides['MgO']))) + 'Mg'
    elif star is not None:
        comp_str = star.replace(' ', '')
    else:
        raise NotImplementedError('no bulk composition scenario input (need star or test_oxides)')
    if 'profile' in kwargs:
        if kwargs['profile'] != 'adiabat':
            temp_str = kwargs['profile']
    elif Tp is not None:
        temp_str = str(int(Tp)) + 'K'
    else:
        raise NotImplementedError('no temperature scenario input (need Tp or profile name)')
    name = mass_str + '_' + cmf_str + '_' + comp_str + '_' + temp_str
    if suffix is not None:
        name = name + '_' + suffix
    name = name.replace('.', ',')  # dots will crash file i/o
    return name


def update_saturation(dat):  # weird choice that this isn't in class file but ok
    """
    gets saturation from composition dataframe
    :param dat:
    :return:
    """
    df_all = dat.df_all
    df_all = sat.mineral_water_contents(df_all['P(bar)'].to_numpy() * 1e5, df_all['T(K)'].to_numpy(), df=df_all)

    # save whole-mantle concentration and mass
    dat.c_h2o_mantle = sat.total_water_frac(df_all)  # weight fraction wrt mantle (NOT PLANET!)
    dat.mass_h2o_total = sat.total_water_mass(df_all)  # total in kg
    print('\n\ntotal water in mantle:', dat.c_h2o_mantle * 1e2, 'wt% = ', dat.mass_h2o_total / sat.TO, 'OM')

    # also isolate upper mantle
    i_um_base = dat.find_lower_mantle() - 1  # base of upper mantle
    dat.mass_h2o_um = sat.total_water_mass(df_all, i_min=0, i_max=i_um_base)

    # also isolate olivine-bearing mantle
    dat.get_obm_water()

    dat.df_all = df_all
    return dat


def build_planet(name=None, get_saturation=True, plot_all=False, plot_kwargs=None, solve_interior=True,
                 **kwargs):
    """ kwargs include
     overwrite (bool) : auto overwrite build etc files,
     head (bool) : to print DataFrame headers when reading them
     n (int) : interior structure radial resolution,
     tol (float) : iteration tolerance for R_p in fraction
     maxIter (int): max iterations
     verbose (bool) : lots of places """
    if plot_kwargs is None:
        plot_kwargs = {}
    time_start = time.time()

    # name object systematically
    if name is None:
        name = get_name(**kwargs)

    print('>>>>>> STARTED', name)

    # create class object
    dat = px.PerplexData(name=name, **kwargs)
    okay = True

    if os.listdir(dat.output_path) and os.path.isfile(dat.output_path + 'dat.pkl'):
        # check if data exists in directpry and load if so
        print(dat.output_path,
              'already contains data, trying to load... [to overwrite, please delete directory manually]')

        # todo cleaner storage - but can write functions to get everything you're not pickling. pickle for tn
        with open(dat.output_path + 'dat.pkl', "rb") as pfile:
            dat_loaded = pkl.load(pfile)

        # check if loaded but output_path is inconsistent (manual directory renaming?)
        if dat.output_path != dat_loaded.output_path:
            dat_loaded.output_path = dat.output_path
        dat = dat_loaded
    else:
        # run --> this is the time-consuming part
        # print('kwargs build_planet', kwargs)
        okay = dat.setup_interior(solve_interior=solve_interior, **kwargs)

    if not okay and solve_interior:
        print('planet is not okay - returning None')
        return None

    elif solve_interior:  # i.e .if meant to solve
        # rename and scale columns (for use with saturation calcs but being consistent regardless)
        df_all = dat.df_comp.copy()
        renamed_phases = rename_phases(dat.phases_px)
        # prefix phase wt composition with 'X_'
        for old, new in zip(dat.phases_px, renamed_phases):
            try:
                col_name = 'X_' + new
                df_all.rename(columns={old: col_name}, inplace=True)
                df_all[col_name] = df_all[col_name] * 1e-2  # also convert from % to frac
            except KeyError as e:
                print('something very wrong:', e)
                print('df columns:', df_all.columns)

        # put important mantle data into single dataframe for convenience - note these are all SI units
        df_all['mass(kg)'] = dat.mass[dat.i_cmb + 1:][::-1]  # invert and just store mantle (surface down)
        df_all['z(m)'] = dat.R_p - dat.radius[dat.i_cmb + 1:][::-1]
        df_all['rho'] = dat.density[dat.i_cmb + 1:][::-1]
        df_all['alpha'] = dat.alpha[dat.i_cmb + 1:][::-1]
        df_all['cp'] = dat.cp[dat.i_cmb + 1:][::-1]
        dat.df_all = df_all

        # print('total mass Ol:', np.sum(df_all['X_ol'] * df_all['mass(kg)']))
        # print('total mass Ring:', np.sum(df_all['X_ring'] * df_all['mass(kg)']))
        # print('total mass Wad:', np.sum(df_all['X_wad'] * df_all['mass(kg)']))

        # calculate saturation water content profile per mineral phase
        if get_saturation:
            dat = update_saturation(dat)

        # calculate some other simple things - only here because you thought of this later but should be done inside class
        if not hasattr(dat, 'mass_um') or dat.mass_um is None:  # needs to be here because df_all not created in class
            dat.get_um_mass()
        if not hasattr(dat, 'mgsi') or dat.mgsi is None:
            dat.get_mgsi()
        if not hasattr(dat, 'mg_number') or dat.mg_number is None:
            dat.get_mg_number()

        # write df to file
        df_all.to_csv(path_or_buf=dat.output_path + dat.name + '_profiles.dat', sep='\t', na_rep='NaN', index=True,
                      index_label=None)

        # move perplex files to their own folder
        for fend in ['_comp.tab', '_adiabat.dat', '.dat', '_WERAMI_options.txt', '_VERTEX_options.txt']:
            file = dat.name + fend
            if os.path.exists(dat.perplex_path + file):
                os.rename(dat.perplex_path + file, dat.output_path + file)

    else:
        print('mgsi', dat.mgsi)
        return dat  # just calculated bulk composition

    # for now also save as pickle for laziness
    with open(dat.output_path + "dat.pkl", "wb") as pfile:
        pkl.dump(dat, pfile)

    time_end = time.time()
    print('\n>>>>>> COMPLETED', dat.name, 'in', time_end - time_start, 's\n\n')

    if plot_all:
        # plot
        plotpx.single_structure(dat, fig_path=dat.output_path, **plot_kwargs)
        plotpx.single_composition(dat, fig_path=dat.output_path, **plot_kwargs)
        if get_saturation:
            plotpx.profile(dat, parameter='c_h2o', independent_ax='p', ax_label='H$_2$O storage capacity (ppm)',
                           scale=1e6, log=True, fig_path=dat.output_path, **plot_kwargs)
    return dat


def build_multi_planets(loop_var_name, loop_var_vals, names=None, **kwargs):
    if names is None:
        # name after loop variable by default
        names = [loop_var_name.replace(" ", "") + str(s) for s in loop_var_vals]
    dat_list = []
    for ii, val in enumerate(loop_var_vals):
        kwargs[loop_var_name] = val  # add to dict, others input or otherwise use default
        dat = build_planet(name=names[ii], **kwargs)
        dat_list.append(dat)
    return dat_list


def get_run_dirs(output_path, prevent_blank_dir=True, verbose=False):
    try:
        subfolders = [f.path for f in os.scandir(output_path) if f.is_dir()]
    except FileNotFoundError:
        # directory doesn't exist
        if prevent_blank_dir:
            raise FileNotFoundError(output_path, 'does not exist')
        if verbose:
            print('warning:', output_path, 'does not exist')
        return []
    except PermissionError:
        print('probably have a for loop reading characters rather than file strings')
    return subfolders


def read_dir(output_path, subsample=None, verbose=False, prevent_blank_dir=True, include_names=None):
    """ read all """
    subfolders = get_run_dirs(output_path, prevent_blank_dir=prevent_blank_dir, verbose=verbose)
    if subsample is not None:
        # a random subset of the directory
        subfolders = random.sample(subfolders, subsample)
    dats = []
    count = 0
    for dir in subfolders:
        try:
            with open(dir + '/dat.pkl', "rb") as pfile:
                dat = pkl.load(pfile)
            dats.append(dat)
            count += 1
        except FileNotFoundError:
            if verbose:
                print('warning:', dir, 'is empty')
        except PermissionError as e:
            print(dir, e)
    if verbose:
        print('loaded', count, 'pickle files from', output_path)
    return dats


def read_name(output_path, name=None, star=None, M_p=None, core_efficiency=None, Tp=None, test_CMF=None, test_oxides=None, suffix=None):
    """ find single folder in output path with name """
    if name is None:
        name = get_name(M_p=M_p, star=star, core_efficiency=core_efficiency, Tp=Tp, test_CMF=test_CMF,
                         test_oxides=test_oxides, suffix=suffix)
    with open(output_path + name + '/dat.pkl', "rb") as pfile:
        dat = pkl.load(pfile)
    return dat


def update_dir(output_path, func, store=False, **func_args):
    """ if you need to redo some analysis or add something """
    dats = read_dir(output_path)
    new_dats = []
    for dat in dats:
        func(dat, **func_args)
        new_dats.append(dat)
        if store:
            # for now save as pickle for laziness
            with open(output_path + dat.name + "/dat.pkl", "wb") as pfile:
                pkl.dump(dat, pfile)
    return new_dats


def planets_from_hypatia(n_sample=-1, M_p=1, names_file='host_names.txt', plot_all=False, restart=None,
                         stopafter=None, **kwargs):
    """ restart is string name of last star that worked"""

    sample_names = hyp.random_star(n_sample, names_file=names_file)
    planets = []

    if restart is not None and n_sample == -1:
        # must be looping all stars for this to work, otfherwise can't work as will be random names
        ii_start = sample_names.index(restart)
    else:
        ii_start = 0
    if stopafter is not None and n_sample == -1:
        ii_last = sample_names.index(stopafter)
        subsample_names = sample_names[ii_start:ii_last + 1]
    else:
        subsample_names = sample_names[ii_start:]
    for ii, star in enumerate(subsample_names):
        print(ii + 1, '/', len(subsample_names))
        # print('kwargs planets_From_hypatia', kwargs)
        pl = build_planet(star=star, M_p=M_p * M_E, plot_all=plot_all, **kwargs)
        pl.write_star_composition(fname='nH_star.txt',
                                  path=pl.output_path)  # save parameter nH_star (abundances) to file
        planets.append(pl)
    return planets
