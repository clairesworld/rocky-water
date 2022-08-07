import numpy as np
import pandas as pd
import os
import pathlib
import subprocess
import perplexdata as px
import bulk_composition as bulk
import main as rw
from scipy import interpolate
import parameters as p

# import ask_hypatia as hyp

Rb = 8.31446261815324
output_parent_default = '/home/claire/Works/min-fo2/perplex_output/'
output_parent_apollo = '/raid1/cmg76/perple_x/output/rocky-fo2/'
wt_oxides_DMM = {'SiO2': 44.71, 'MgO': 38.73, 'CaO': 3.17, 'Al2O3': 3.98, 'FeO': 8.18}  # Workman & Hart depleted mantle
solution_phases_default = ['Cpx(JH)', 'O(JH)', 'Sp(JH)', 'Grt(JH)', 'Opx(JH)']

wt_oxides_DMM_ext = {'SiO2': 44.71, 'MgO': 38.73, 'CaO': 3.17, 'Al2O3': 3.98, 'FeO': 8.17,
                     'Na2O': 0.13, 'Cr2O3': 0.57}  # Workman & Hart depleted mantle, extended

class PerplexFugacity(px.PerplexData):

    def __init__(self, name='test', star=None, solution_phases=solution_phases_default,
                 output_parent_path=output_parent_default, wt_oxides=wt_oxides_DMM, X_ferric=0.03,
                 perplex_path=px.perplex_path_default,
                 **kwargs):
        # super().__init__(output_parent_path=output_parent_path, **kwargs)

        self.name = name
        self.star = star

        if 'O2' not in wt_oxides:
            # get O2 mantle fraction
            self.wt_oxides = bulk.o2_from_wt_comp(wt_oxides, X_ferric=X_ferric)
        else:
            self.wt_oxides = wt_oxides
        self.oxide_list = [k for k in self.wt_oxides.keys()]
        self.get_mgsi()
        self.get_mg_number()
        self.solution_phases = solution_phases
        self.output_path = output_parent_path + self.name + '/'
        self.perplex_path = perplex_path

    def command_vertex_grid(self, build_file_end=''):
        """ string for vertex command file - only entry is project name in operational mode 5 """
        s = self.name + build_file_end
        return s

    def command_werami_mu(self, points_file=px.perplex_path_default + 'Earth_1600K_adiabat.dat', build_file_end='',
                          **kwargs):
        """ string for werami command file to get chemical potentials """

        s = self.name + build_file_end + '\n'  # Enter the project name (the name assigned in BUILD)
        s = s + '4\n'  # Select operational mode: 4 - as in 3 [properties along a 1d path], but input from file
        s = s + '2\n'  # Path will be described by:  2 - a file containing a list of x-y points
        s = s + points_file + '\n'  # Enter the file name:
        s = s + '1\n'  # every nth plot will be plotted, enter n:
        s = s + '23\n'  # Select a property: 23 - Dependent potentials (J/mol, bar, K)
        s = s + str(self.oxide_list.index('O2') + 1) + '\n'  # Enter a component:  [here index starts at 1]
        s = s + '0\n'  # Select an additional property or enter 0 to finish:
        s = s + '0\n'  # Select operational mode: 0 - EXIT
        s = s + 'EOF'
        return s

    def command_werami_composition_grid(self, points_file=px.perplex_path_default + 'Earth_1600K_adiabat.dat',
                                        build_file_end='', **kwargs):
        """ string for werami command file to get chemical potentials """

        s = self.name + build_file_end + '\n'  # Enter the project name (the name assigned in BUILD)
        s = s + '4\n'  # Select operational mode: 4 - as in 3 [properties along a 1d path], but input from file
        s = s + '2\n'  # Path will be described by:  2 - a file containing a list of x-y points
        s = s + points_file + '\n'  # Enter the file name:
        s = s + '1\n'  # every nth plot will be plotted, enter n:
        s = s + '25\n'  # Select a property: 25 - Modes of all phases
        s = s + 'n\n'  # Output cumulative modes (y/n)?
        s = s + 'n\n'  # Include fluid in computation of aggregate( or modal) properties(y / n)?
        s = s + '0\n'  # Select operational mode: 0 - EXIT
        return s

    def command_frendly_mu(self, T_min, T_max, build_file_end='', vertex_data='hp622ver', **kwargs):
        """ string for frendly command file to get chemical potentials """

        s = vertex_data + '.dat\n'  # Enter thermodynamic data file name [default = hp02ver.dat]:
        s = s + 'claire\n'  # What is your name?
        s = s + 'y\n'  # Can you say "therm-o-dy-nam-ics" (y/n)?
        s = s + '2\n'  # Choose from the following options:  2 - [default] thermodynamic properties for a phase or reaction relative to the reference state.
        s = s + 'n\n'  # List database phases (y/n)?
        s = s + 'n\n'  # Calculate thermodynamic properties for a reaction (y/n)?
        s = s + 'O2\n'  # Calculate thermodynamic properties for phase:
        s = s + '1\n'  # Enter activity [default is 1] of O2
        s = s + 'y\n'  # Write a properties table (Y/N)?
        s = s + 'n\n'  # Compute properties along a path (Y/N)?
        s = s + 'y\n'  # Make a 1-dimensional (e.g., isobaric) table (y/n)?
        s = s + '2\n'  # Select the independent (x) variable: 2 - T(K)
        s = s + str(T_min) + ' ' + str(T_max) + ' 1\n'  # Enter minimum, maximum, and increment for T(K)    :
        s = s + '1\n'  # Specify the value for: P(bar)
        s = s + self.name + build_file_end + '_standard\n'  # Enter the output file name [without the .plt/.tab suffix, default is my_project]:
        s = s + 'O2 potential\n'  # Enter calculation title:
        s = s + '4'  # Choose from the following options: 4 - quit.
        return s

    def command_frendly_mu_qfm(self, T_min, T_max, p0=1, vertex_data='hp622ver_qfm', fout=None, **kwargs):
        """ string for frendly command file to get delta G of FMQ reaction, p0 in bar """

        if fout is None:
            fout = 'qfm_frendly_' + str(p).replace('.', ',') + 'bar'

        s = vertex_data + '.dat\n'  # Enter thermodynamic data file name [default = hp02ver.dat]:
        s = s + 'claire\n'  # What is your name?
        s = s + 'y\n'  # Can you say "therm-o-dy-nam-ics" (y/n)?
        s = s + '2\n'  # Choose from the following options:  2 - [default] thermodynamic properties for a phase or reaction relative to the reference state.
        s = s + 'n\n'  # List database phases (y/n)?
        s = s + 'y\n'  # Calculate thermodynamic properties for a reaction (y/n)?
        s = s + '4\n'  # How many phases or species in the reaction?
        s = s + 'fa\n'  # Enter phase or species number  1 in your reaction:
        s = s + '+3\n'  # Enter reaction coefficient for: fa       products (+), reactants (-):
        s = s + '1\n'  # Enter activity [default is 1] of fa
        s = s + 'O2\n'  # Enter phase or species number  2 in your reaction:
        s = s + '+1\n'  # Enter reaction coefficient for: O2       products (+), reactants (-):
        s = s + '1\n'  # Enter activity [default is 1] of O2
        s = s + 'mt\n'  # Enter phase or species number  3 in your reaction:
        s = s + '-2\n'  # Enter reaction coefficient for: mt       products (+), reactants (-):
        s = s + '1\n'  # Enter activity [default is 1] of mt
        s = s + 'q\n'  # Enter phase or species number  4 in your reaction:
        s = s + '-3\n'  # Enter reaction coefficient for: q       products (+), reactants (-):
        s = s + '1\n'  # Enter activity [default is 1] of q
        s = s + 'y\n'  # Write a properties table (Y/N)?
        s = s + 'n\n'  # Compute properties along a path (Y/N)?
        s = s + 'y\n'  # Make a 1-dimensional (e.g., isobaric) table (y/n)?
        s = s + '2\n'  # Select the independent (x) variable: 2 - T(K)
        s = s + str(T_min) + ' ' + str(T_max) + ' 1\n'  # Enter minimum, maximum, and increment for T(K)    :
        s = s + str(p0) + '\n'  # Specify the value for: P(bar)
        s = s + fout + '\n'  # Enter the output file name [without the .plt/.tab suffix, default is my_project]:
        s = s + 'QFM\n'  # Enter calculation title:
        s = s + '4'  # Choose from the following options: 4 - quit.
        return s

    def run_frendly(self, T_min, T_max, frendly_command_end='_frendly_command.txt', build_file_end='',
                    suppress_output=True, clean=True, verbose=True, frendly_command_text_fn=None, frendly_kwargs=None,
                    **kwargs):

        if frendly_kwargs is None:
            frendly_kwargs = {}
        cwd = os.getcwd()
        os.chdir(self.perplex_path)
        if suppress_output:
            stderr, stdout = subprocess.DEVNULL, subprocess.DEVNULL
        else:
            stderr, stdout = None, None

        # create frendly command file
        frendly_command_file = self.name + build_file_end + frendly_command_end
        with open(self.perplex_path + frendly_command_file, 'w') as file:
            s = frendly_command_text_fn(T_min, T_max, build_file_end=build_file_end, **frendly_kwargs)
            file.write(s)

        # run frendly
        output = subprocess.run('./frendly < ' + frendly_command_file, shell=True,
                                stdout=stdout, stderr=stderr)
        if verbose:
            print('  ', output)

        # delete extra files
        if clean:
            os.remove(frendly_command_file)
        os.chdir(cwd)  # return to original dir

    def fo2_calc(self, p_min=1000, p_max=40000, T_min=1600, T_max=1603, isotherm=None, verbose=True, build_file_end='',
                 vertex_data='hp622ver', run=True, compare_buffer=None, check_comp=False, points_file=None,
                 mu0_file=None, save=True, **kwargs):
        """
        Parameters
        ----------
        p_min :
        p_max :
        T_min :
        T_max :
        isotherm : Numeric or None
            Temperature value for doing calculations for an isothermal section, None if not applicable
        verbose :
        build_file_end :
        vertex_data :
        run :
        compare_buffer :
        check_comp :
        points_file :
        mu0_file : Path or string
            Given relative to perplex_path
        save : bool
            If True, auto-save important results to csv

        Returns
        -------

        """
        if verbose:
            print(self.__dict__)

        # decide which T, p points to calculate
        if isotherm is not None:
            if points_file is not None:
                raise Exception('ERROR: cannot take both points_file and isotherm, set one of these to None')
            try:
                # run werami over an isothermal section
                points_file = self.write_isothermal_points(T0=isotherm, p_min=p_min, p_max=p_max, delta_p=1000)
            except TypeError as e:
                print('ERROR: trying to write an isotherm section but isotherm is not a valid numeric')
                raise e
        elif points_file is None:
            raise NotImplementedError('ERROR: need either input points_file or input isotherm')

        if run:
            # create build file
            self.write_build(title='Planet', p_min=p_min, p_max=p_max,  # adiabat_file=points_file,
                             verbose=verbose, overwrite=True, vertex_data=vertex_data, option_file='perplex_option_fo2',
                             excluded_phases=None, use_solutions=True, build_file_end=build_file_end,
                             calculation_type='5',  # '5' for grid
                             T_min=T_min, T_max=T_max, **kwargs)

            # run vertex and werami to get mu at T, P of interest
            # strategy is to run each case once in vertex over grid and keep all the files
            self.run_perplex(werami_command_end='_werami_command_mu.txt',
                             werami_command_text_fn=self.command_werami_mu,
                             vertex_command_text_fn=self.command_vertex_grid,
                             output_file_end='_mu.tab', build_file_end=build_file_end, verbose=verbose, clean=True,
                             werami_kwargs={'points_file': points_file}, run_vertex='auto',
                             store_vertex_output=True,
                             **kwargs)

            if check_comp:
                self.run_perplex(werami_command_end='_werami_command_comp.txt',
                                 werami_command_text_fn=self.command_werami_composition_grid,
                                 vertex_command_text_fn=self.command_vertex_grid,
                                 output_file_end='_comp.tab', build_file_end=build_file_end, verbose=verbose,
                                 clean=True, werami_kwargs={'points_file': points_file}, run_vertex='auto',
                                 **kwargs)

        # extract mu and T of interest from output data
        df_w = self.read_chem_potential_werami(verbose=False, return_df=True)
        mu = df_w['mu_O2(J/mol)'].to_numpy()
        T = df_w['T(K)'].to_numpy()
        P = df_w['P(bar)'].to_numpy()

        # prepare all data
        df_save = df_w.copy()

        if check_comp:
            df_c = self.read_werami(fend='_comp.tab')

            # concat with full results dataframe
            px_phases = [col for col in df_c.columns if col not in ['P(bar)', 'T(K)']]
            print('px_phases', px_phases)
            print({phase: 'X_' + phase for phase in px_phases})
            df_c.rename(columns={phase: 'X_' + phase for phase in px_phases}, inplace=True)
            df_save = pd.concat([df_w, df_c], axis=1)

        if mu0_file is None and run:
            # run frendly to get mu_0 at T of interest, P = 1 bar
            self.run_frendly(T.min(), T.max(), build_file_end=build_file_end,
                             frendly_kwargs={'vertex_data': vertex_data},
                             frendly_command_text_fn=self.command_frendly_mu, verbose=verbose, clean=True, **kwargs)
            mu0_file = self.name + build_file_end + '_standard.tab'

        # extract mu_0 and T of interest from output data
        mu_0, T_mu_0 = self.read_o2_potential_frendly(T0=T, fin=mu0_file, verbose=False)
        # print('mu =', mu, '\nmu_0 =', mu_0, '\nT_mu_0 =', T_mu_0)

        # calculate log10(fo2) from standard state potential
        logfo2 = mu_to_logfo2(mu, mu_0, T)  # mu is P-dependent here, given a range of P
        df_save['logfo2'] = logfo2

        if compare_buffer == 'qfm':
            # run frendly to get fo2 of qfm buffer at T, p
            # todo: this is accurate at 1 bar but fo2 gets way too low at 2 GPa - do you need to subtract from standard state?
            # p0 = 1  # p[0]
            # fout = 'qfm_frendly_' + str(p0).replace('.', ',') + 'bar'
            # self.run_frendly(T.min(), T.max(),
            #                  frendly_command_text_fn=self.command_frendly_mu_qfm, verbose=verbose, clean=False,
            #                  frendly_kwargs={'vertex_data': 'hp622ver_qfm', 'p': p0, 'fout': fout},
            #                  **kwargs)
            # df_rx = pd.read_csv(self.perplex_path + fout + '.tab', skiprows=12,
            #                     index_col=None,
            #                     sep=r"\s+")
            # deltaG = df_rx['g(J/mol)']
            # logfo2_buffer = -deltaG[0] / (np.log(10) * Rb * T.min())

            logfo2_buffer = read_qfm_os(T, P, verbose=False)
            # print('logfo2_qfm', logfo2_buffer)
            # print('logfo2 = QFM + ', logfo2 - logfo2_buffer, 'at', P, 'bar,', T, 'K')
            df_save['logfo2_qfm'] = logfo2_buffer
            df_save['delta_qfm'] = logfo2 - logfo2_buffer

        # store mega df
        df_save = df_save.loc[:, ~df_save.columns.duplicated()].copy()
        if save:
            df_save.to_csv(self.output_path + self.name + '_results.csv', sep="\t")

        # store werami/frendly i/o files - note must be after reading-in is finished
        for fend in ['_WERAMI_options.txt', '_VERTEX_options.txt', '_mu.tab', '_comp.tab', '_standard.tab', '.dat']:
            if os.path.isfile(self.perplex_path + self.name + build_file_end + fend):
                os.rename(self.perplex_path + self.name + build_file_end + fend,
                          self.output_path + self.name + build_file_end + fend)

        return logfo2

    def read_o2_potential_frendly(self, build_file_end='', T0=None, fin=None, verbose=False, return_df=False):
        if fin is None:
            fin = self.name + build_file_end + '_standard.tab'
        try:
            df = pd.read_csv(self.perplex_path + fin, skiprows=12,
                             index_col=None,
                             sep=r"\s+")
        except FileNotFoundError as e:
            try:
                df = pd.read_csv(self.output_path + fin, skiprows=12,
                                 index_col=None,
                                 sep=r"\s+")
            except FileNotFoundError as e:
                print(e)
                raise Exception('Missing frendly output files')
        # df.rename(columns={"g(J/mol)": "mu_O2(J/mol)"}, inplace=True)
        if verbose:
            print(df.head())
        if T0 is not None:
            # return row for specific T
            try:
                df = df[df['T(K)'].isin(T0)]  # note this should work for an array of T0 as well
            except:
                df = df[df['T(K)'] == T0]
        mu_0 = df['g(J/mol)'].to_numpy()
        T = df['T(K)'].to_numpy()
        if return_df:
            return df
        else:
            return np.squeeze(mu_0), np.squeeze(T)

    def read_werami(self, fend=''):
        axis = 1
        try:
            df_w = pd.read_csv(self.perplex_path + self.name + fend, skiprows=8, index_col=None, sep=r"\s+")
        except FileNotFoundError as e:
            try:
                df_w = pd.read_csv(self.output_path + self.name + fend, skiprows=8, index_col=None, sep=r"\s+")
            except FileNotFoundError as e:
                print(e)
                raise Exception('Missing vertex/werami output files')
        return df_w

    def read_chem_potential_werami(self, build_file_end='', verbose=False, return_df=False):
        df = self.read_werami(fend=build_file_end + '_mu.tab')
        df.rename(columns={"mu[O2],J/mol": "mu_O2(J/mol)"}, inplace=True)
        if verbose:
            print(df.head())
        if return_df:
            return df
        else:
            return df['mu_O2(J/mol)'].to_numpy(), df['T(K)'].to_numpy()

    def write_isothermal_points(self, T0, p_min, p_max, delta_p=1000, fname=None):
        """ write sample points in p, T with p in bar, T in K """
        if fname is None:
            fname = str(T0).replace(',', '.') + 'K_isotherm.xy'
        if not os.path.isfile(self.perplex_path + 'points_files/' + fname):
            P = np.arange(p_min, p_max, step=delta_p)
            T = T0 * np.ones_like(P)
            self.write_adiabat(P, T, fout=self.perplex_path + 'points_files/' + fname, verbose=False, overwrite=True)
        return 'points_files/' + fname


def init_from_build(name, **kwargs):
    return PerplexFugacity(**read_from_build(name, **kwargs))


def read_from_build(name, output_parent_path=output_parent_default, verbose=False, **kwargs):
    import re
    # parse the parameters file into a python dictionary

    fin = output_parent_path + name + '/' + name + '.dat'
    if verbose:
        print("Reading build parameters from", fin)

    build_file = open(fin).readlines()

    # read thermodynamic data (always 1st line?)
    vertex_data = build_file[0].split()[0]

    # read thermodynamic components
    wt_oxides = {}
    idx = [x for x in range(len(build_file)) if 'begin thermodynamic component list' in build_file[x]][
        0]  # start of list
    line = build_file[idx + 1]
    while 'end thermodynamic component list' not in line:
        parts = line.split()  # split arbitrary whitespace
        wt_oxides[parts[0]] = float(parts[2])
        idx += 1
        line = build_file[idx + 1]

    # read excluded phases
    excluded_phases = []
    idx = [x for x in range(len(build_file)) if 'begin excluded phase list' in build_file[x]][0]  # start of list
    line = build_file[idx + 1]
    while 'end excluded phase list' not in line:
        excluded_phases.append(line.rstrip())
        idx += 1
        line = build_file[idx + 1]

    # read solution phases
    solution_phases = []
    idx = [x for x in range(len(build_file)) if 'begin solution phase list' in build_file[x]][0]  # start of list
    line = build_file[idx + 1]
    while 'end solution phase list' not in line:
        solution_phases.append(line.rstrip())
        idx += 1
        line = build_file[idx + 1]

    # read min, max (p, T)
    idx = [x for x in range(len(build_file)) if 'max p' in build_file[x]][0]  # start of list
    line_max = build_file[idx].strip().split()  # split arbitrary whitespace
    line_min = build_file[idx + 1].strip().split()  # split arbitrary whitespace
    p_max, T_max = line_max[0], line_max[1]
    p_min, T_min = line_min[0], line_min[1]

    # compile in args dict
    d = {'vertex_data': vertex_data, 'wt_oxides': wt_oxides, 'excluded_phases': excluded_phases,
         'solution_phases': solution_phases, 'p_min': p_min, 'p_max': p_max, 'T_min': T_min, 'T_max': T_max}

    # TODO parse star from name using possible prefixes
    if verbose:
        print('d', d)
    return d

    # d = {}
    # keys = []
    #
    # def nested_set(dic, keys, value):
    #     for k in keys[:-1]:
    #         dic = dic[k]
    #     dic[keys[-1]] = value
    #
    # def num(s):
    #     try:
    #         return int(s)
    #     except ValueError:
    #         try:
    #             return float(s)
    #         except ValueError:
    #             return s

    # for ii, line in enumerate(build_file):
    #     line, _, comments = line.partition('|')  # separate comments
    #     line = line.strip()
    #     if line:
    #         parts = line.split()  # splot by whitespace
    #         if parts[1] == 'thermodynamic data file':
    #             vertex_data = parts[0]
    #
    #         if parts[1] == 'thermodynamic' and parts[2] == 'component':
    #

    # if re.match("begin(.*)", line):
    #     sub_name = line[5:].strip()
    #     keys.append(sub_name)
    #     nested_set(d, keys, {})
    #
    # if re.match("end(.*)", line):
    #     keys.pop()
    #
    # if re.match("set(.*)", line):
    #     i = line.index('=')
    #     key = line[4:i].strip()
    #     value = line[i + 1:].strip()
    #     keys.append(key)
    #     nested_set(d, keys, num(value))
    #     keys.pop()
    # return d


def mu_to_logfo2(mu, mu_0, T):
    return np.log(np.exp((mu - mu_0) / (Rb * T))) / np.log(10)


def read_qfm_os(T, P, path=px.perplex_path_default, fin='data_tables/fmqNl_fo2_oli.dat', verbose=False):
    """ read in oli's qfm calculations, T in K, p in bar """

    def do(df, T0, p0):
        # first check if exact values already there
        exists = None
        if p0 in df['P(bar)'].unique():
            exists = 'P'
            df = df[df['P(bar)'] == p0]

        if T0 in df['T(K)'].unique():
            if exists:  # both T and P already there
                return np.squeeze(df[(df['P(bar)'] == p0) and (df['T(K)'] == T0)]['logfo2'])
            exists = 'T'
            df = df[df['T(K)'] == T0]

        # too many points to interpolate so find 4 neares
        if exists != 'P':
            # get nearest pressures
            tmp1 = p0 - df['P(bar)']  # want smallest nonnegative number so p1 < p0
            idx1 = tmp1.mask(tmp1 < 0).argmin()
            tmp2 = df['P(bar)'] - p0  # want smallest nonnegative number so p0 < p2
            idx2 = tmp2.mask(tmp2 < 0).argmin() + nT
            df = df.iloc[idx1:idx2 + 1]

        if exists != 'T':
            # get nearest temperatures
            tmp1 = T0 - df['T(K)']  # want smallest nonnegative number so p1 < p0
            idx1 = tmp1.mask(tmp1 < 0).argmin()
            T1 = df['T(K)'].iloc[idx1]
            tmp2 = df['T(K)'] - T0  # want smallest nonnegative number so p0 < p2
            idx2 = tmp2.mask(tmp2 < 0).argmin()
            T2 = df['T(K)'].iloc[idx2]
            df = df[df['T(K)'].isin([T1, T2])]

        # interpolate in this range
        if verbose:
            print('starting interpolation over df length', len(df), '| p_min', df['P(bar)'].iloc[0], 'p_max',
                  df['P(bar)'].iloc[-1], 'p0', p0,
                  '| T_min', df['T(K)'].iloc[0], 'T_max', df['T(K)'].iloc[-1], 'T0', T0)

        x_in = df['P(bar)'].to_numpy()
        y_in = df['T(K)'].to_numpy()
        z_in = df['logfo2'].to_numpy()
        # f = interpolate.griddata((x_in, y_in),

        if exists == 'P':  # interpolate over temperatures
            f = interpolate.interp1d(x=y_in, y=z_in)
            return np.squeeze(f(T0))
        elif exists == 'T':  # interpolate over pressures
            f = interpolate.interp1d(x=x_in, y=z_in)
            return np.squeeze(f(p0))
        else:  # interpolate 2D
            f = interpolate.interp2d(x=x_in, y=y_in, z=z_in)
            return np.squeeze(f(p0, T0))

    df = pd.read_csv(path + fin, delimiter=r"\s+", index_col=False, header=None,
                     names=['P(bar)', 'T(K)', 'logfo2', 'thing1', 'thing2', 'Rb'])
    df = df.drop(columns=['thing1', 'thing2', 'Rb'])
    nT = 2023 - 1323  # number of temperature points

    T_is_mult = np.size(T) > 1
    P_is_mult = np.size(P) > 1

    if T_is_mult:
        if (T[:-1] == T[1:]).all():  # if all elements equal
            T = T[0]
            T_is_mult = False
    if P_is_mult:
        if (P[:-1] == P[1:]).all():  # if all elements equal
            P = P[0]
            P_is_mult = False

    if T_is_mult:  # for input multiple values
        if P_is_mult:
            raise NotImplementedError('ERROR: cannot interpolate buffer fo2 for gridded P, T of interest')
        else:
            logfo2 = []
            for TT in T:  # loop over each T for single p
                logfo2.append(do(df, TT, P))
    elif P_is_mult:  # copy above
        logfo2 = []
        for pp in P:  # loop over each p for single T
            logfo2.append(do(df, T, pp))
    else:
        logfo2 = do(df, T, P)

    return logfo2


def fo2_from_hypatia(p_min, p_max, n_sample=5, core_efficiency=0.88, planet_kwargs={},
                     output_parent_path=output_parent_default, **kwargs):
    planet_kwargs.update({'core_efficiency': core_efficiency})
    pl_list = rw.planets_from_hypatia(n_sample=n_sample, names_file='host_names.txt', plot_all=False, restart=None,
                                      get_saturation=False, solve_interior=False,
                                      stopafter=None, output_parent_path=output_parent_path,
                                      **planet_kwargs, **kwargs)
    for pl in pl_list:
        fo2_from_oxides(pl.name, p_min, p_max, pl=pl, output_parent_path=output_parent_path, **kwargs)


def fo2_from_local(output_parent_path=output_parent_default, **kwargs):
    # perform fo2 calculations on local (existing) vertex data
    subfolders = rw.get_run_dirs(output_path=output_parent_path)
    if subfolders:
        for sub in subfolders:
            name = os.path.basename(sub)
            # star = name.split('_')[2]

            d = read_from_build(name=name, output_parent_path=output_parent_path, verbose=False)
            dat = PerplexFugacity(name=name, output_parent_path=output_parent_path, **d, **kwargs)
            logfo2 = dat.fo2_calc(run=False, **d, **kwargs)
            print('log fo2 of system:', logfo2)
    else:
        print('no local output found in', output_parent_path)


def fo2_from_oxides(name, p_min, p_max, T_min=1373, T_max=1900, pl=None,
                    core_efficiency=0.88, test_oxides=None, star=None, verbose=True, planet_kwargs={},
                    output_parent_path=output_parent_default, **kwargs):
    """ p_min and p_max in bar """

    if pl is None:
        # make planet object without solving interior structure, but getting bulk composition

        if star is not None and test_oxides is not None:
            raise NotImplementedError('Are you sure you want to input both test_oxides and star? Only need one.')
        planet_kwargs.update({'test_oxides': test_oxides, 'core_efficiency': core_efficiency})
        pl = rw.build_planet(name=name, star=star, get_saturation=False, solve_interior=False, verbose=verbose,
                             output_parent_path=output_parent_path, **planet_kwargs)

    print('Starting fo2 calc for planet of', pl.star, 'with', [(v, k) for k, v in pl.wt_oxides.items()])
    dat = PerplexFugacity(name=name, wt_oxides=pl.wt_oxides, verbose=verbose, output_parent_path=output_parent_path,
                          **kwargs)

    logfo2 = dat.fo2_calc(p_min=p_min, p_max=p_max, T_min=T_min, T_max=T_max,
                          verbose=verbose, **kwargs)

    print('log fo2 of system:', logfo2)
