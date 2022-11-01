import numpy as np
import pandas as pd
import py.perplexdata as px
import py.bulk_composition as bulk
import os
import subprocess
import py.main as rw
import py.minfo2.perplexfugacitydata as pf
import py.parameters as p

output_parent_default = '/home/claire/Works/min-fo2/alphamelts_output/'
alphamelts_path_default = '/home/claire/Works/alphamelts/'


# solution_phases_default = ['olivine', 'spinel', 'garnet', 'orthopyroxene', 'clinopyroxene']  # opx and cpx automatic?

class MeltsFugacityData:

    def __init__(self, name='test', star=None,  # solution_phases=solution_phases_default,
                 output_parent_path=output_parent_default, wt_oxides=pf.wt_oxides_DMM, X_ferric=0.03,
                 alphamelts_path=alphamelts_path_default,
                 T_final=None, pressures_of_interest=None,
                 **kwargs):

        """
        because pMELTS is bad at starting from subsolidus, need to execute isobarically for each desired
        pressure of interest. therefore output is structured by :

                        individual composition
                                 |
                       --------------------
                       |                  |
                      p[0]     ...       p[n]
        """

        if (T_final is None) or (pressures_of_interest is None):
            raise Exception('ERROR: Must assign a T_final and pressures_of_interest list directly.')

        try:
            pressures_of_interest[0]
        except TypeError:
            pressures_of_interest = [pressures_of_interest]  # must be list

        self.name = name
        self.star = star

        if 'Fe2O3' not in wt_oxides:
            # get O2 mantle fraction
            self.wt_oxides = bulk.insert_fe2o3(wt_oxides, X_ferric, mass_tot=100)
        else:
            self.wt_oxides = wt_oxides
        self.oxide_list = [k for k in self.wt_oxides.keys()]
        # self.solution_phases = solution_phases
        self.get_mgsi()
        self.get_mg_number()

        # initialise dataframe for storing important things
        self.df_all = pd.DataFrame(columns=['P(bar)', 'T(K)', 'logfo2'], index=range(len(pressures_of_interest)))
        self.alphamelts_path = alphamelts_path

        # set main output path for this composition and subpaths for pressure runs
        self.output_path = output_parent_path + self.name + '/'
        self.output_p_paths = [self.output_path + str(pp).replace('.', ',') + 'bar/' for pp in pressures_of_interest]
        self.pressures_of_interest = pressures_of_interest
        self.T_final = T_final

        print('\ninitialising MeltsFugacityData data with:')
        print('        T =', self.T_final, 'K')
        print('        p =', self.pressures_of_interest, 'bar')
        print('\n')

    def get_mgsi(self, **kwargs):
        # if self.nH_star:
        #     self.mgsi = 10 ** self.nH_star[self.oxide_list.index('MgO')] / 10 ** self.nH_star[self.oxide_list.index('SiO2')]  # molar ratio
        # else:
        n_MgO = self.wt_oxides['MgO'] / p.M_MgO  # convert to moles
        n_SiO2 = self.wt_oxides['SiO2'] / p.M_SiO2
        self.mgsi = n_MgO / n_SiO2  # n_Mg = n_MgO etc
        return self.mgsi

    def get_mg_number(self):
        # Mg number of mantle
        n_MgO = self.wt_oxides['MgO'] / p.M_MgO  # convert to moles
        n_FeO = self.wt_oxides['FeO'] / p.M_FeO
        self.mg_number = n_MgO / (n_MgO + n_FeO) * 100
        return self.mg_number

    def batchfile_text_fo2calc(self, melts_file):
        """ contents to write to alphamelts batchfile for fo2 calc """
        s = '1\n'  # 1. Read MELTS file to set composition of system
        s = s + melts_file + '\n'  # MELTS filename:
        s = s + '4\n'  # 4. Execute (follow path, mineral isograd or melt contour)
        s = s + '1\n'  # Superliquidus (1) or subsolidus (0) initial guess ?
        s = s + '0\n'
        return s

    def meltsfile_text_fo2calc(self, p_of_interest, T_ini=2000, **kwargs):
        """ contents to write to alphamelts meltsfile for fo2 calc
        Parameters
        ----------
        p_of_interest :
            Desired pressure of interest in bar
        T_ini : float or int
            Initial temperature for pMELTS execution in K, should be superliquidus

        Returns
        -------
        s : string
            string to write to .melts file
        """

        s = 'Title: ' + self.name + '\n'
        for k, v in self.wt_oxides.items():
            s = s + 'Initial Composition: ' + k + ' ' + str(v) + '\n'
        s = s + 'Initial Temperature: ' + str(T_ini) + '\n'
        s = s + 'Initial Pressure: ' + str(p_of_interest) + '\n'
        return s

    def run_alphamelts_all_p(self, **kwargs):
        for newpath in [self.output_path] + self.output_p_paths:
            if not os.path.exists(newpath):
                # create director(ies) if doesn't exist yet
                print('creating directory', newpath)
                os.makedirs(newpath)
        for (p_of_interest, ppath) in zip(self.pressures_of_interest, self.output_p_paths):
            self.run_alphamelts_at_p(p_of_interest=p_of_interest, output_p_path=ppath,
                                     batch_text_fn=self.batchfile_text_fo2calc,
                                     melts_text_fn=self.meltsfile_text_fo2calc, **kwargs)

    def run_alphamelts_at_p(self, p_of_interest=None, output_p_path=None, suppress_output=True, clean=True,
                            verbose=True,
                            batch_text_fn=None, melts_text_fn=None,
                            env_file=None,
                            melts_kwargs=None, overwrite=False, verify_on_path=False, **kwargs):

        if output_p_path is None:
            output_p_path = self.output_path + str(p_of_interest).replace('.', ',') + 'GPa/'
        if melts_kwargs is None:
            melts_kwargs = {}
        if env_file is None:
            env_file = self.alphamelts_path + 'examples/alphamelts_env_isobaric.txt'

        # check if run already exists
        if (not os.path.isfile(output_p_path + 'System_main_tbl.txt')) or overwrite:

            os.chdir(self.alphamelts_path)

            if suppress_output:
                stderr, stdout = subprocess.DEVNULL, subprocess.DEVNULL
            else:
                stderr, stdout = None, None

            melts_file = output_p_path + self.name + '.melts'
            batch_file = output_p_path + self.name + '.in'

            # create batch and melts files
            with open(batch_file, 'w') as file:
                file.write(batch_text_fn(melts_file))
            with open(melts_file, 'w') as file:
                file.write(melts_text_fn(p_of_interest, **melts_kwargs))

            # make sure installation is on path
            if verify_on_path:
                subprocess.call(['sh', './verify_path.sh'])  # should be in alphamelts directory

            # run alphamelts
            cmd = 'run_alphamelts.command -p "{}" -f "{}" -m "{}" -b "{}"'.format(output_p_path, env_file,
                                                                                  melts_file, batch_file)
            output = subprocess.run(cmd, shell=True, stdout=stdout, stderr=stderr)
            if verbose:
                print(output)

            if clean:
                # get rid of unneeded files
                for fend in ['alphaMELTS_tbl.txt', 'Bulk_comp_tbl.txt', 'Liquid_comp_tbl.txt', 'Solid_comp_tbl.txt',
                             'Trace_main_tbl.txt']:
                    if os.path.isfile(output_p_path + fend):
                        os.remove(output_p_path + fend)
                    else:
                        print('cannot find file: ', output_p_path + fend)
                        raise FileNotFoundError(
                            'ERROR: alphamelts did not complete, try running again with suppress_output=False')

            # return to original dir
            os.chdir(os.path.dirname(os.path.abspath(__file__)))
        else:
            print('Run', self.name, 'at p =', p_of_interest,
                  'bar already exists! To execute, delete files or set overwrite=True\n---------------------------------------------')

    def read_melts_TP(self, check_isothermal=True, reload_TP=False):

        # check if already loaded
        if (not self.df_all['P(bar)'].isnull().values.any()) or (not self.df_all['T(K)'].isnull().values.any()):
            if not reload_TP:
                print('already loaded T,P from alphaMELTS! to overwrite, use reload_TP=True')
                return None

        fname = 'System_main_tbl.txt'
        for row, path in enumerate(self.output_p_paths):
            # load output csv from pMELTS
            output_file = path + fname
            df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                             dtype=np.float64).tail(1)

            # append P and T
            if np.isnan(self.df_all.loc[row, 'P(bar)']):
                self.df_all.loc[row, 'P(bar)'] = df['Pressure'].iloc[-1]
            elif self.df_all.loc[row, 'P(bar)'] != df['Pressure'].iloc[-1]:
                error_str = 'Error: pressure mismatch!\nLoaded {} from {}\nHave {} in self.df_all'.format(
                    df['Pressure'].iloc[-1], output_file, self.df_all['P(bar)'].iloc[row])
                raise RuntimeError(error_str)
            if np.isnan(self.df_all['T(K)'].iloc[row]):
                self.df_all.loc[row, 'T(K)'] = df['Temperature'].iloc[-1]  # K
            elif self.df_all.loc[row, 'T(K)'] != df['Temperature'].iloc[-1]:
                error_str = 'Error: pressure mismatch!\nLoaded {} from {}\nHave {} in self.df_all'.format(
                    df['Temperature'].iloc[-1], output_file, self.df_all['T(K)'].iloc[row])
                raise RuntimeError(error_str)

        if check_isothermal:
            # check if all runs did indeed complete to desired T_final
            print('TODO: check isothermal')
            pass

    def read_melts_phases(self, which='mass', **kwargs):
        """ make csv of phase proportions in same format as perple_x - isothermal x-section for several pressures """

        if which == 'mass':
            fname = 'Phase_mass_tbl.txt'
        elif which == 'volume':
            fname = 'Phase_vol_tbl.txt'

        self.read_melts_TP(**kwargs)

        for row, path in enumerate(self.output_p_paths):
            # load output csv from pMELTS
            output_file = path + fname
            df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                             dtype=np.float64).tail(1)
            # print('loaded\n', df.head())

            # append phases to self df
            m_tot = df['mass'].iloc[-1]
            for ph in [col for col in df.columns if col.endswith('_0')]:
                if ph != 'liquid_0':
                    try:
                        self.df_all.loc[row, ph] = df[ph].iloc[
                                                       -1] / m_tot * 100  # renormalise to 100 g total mass, only need last row
                    except KeyError:
                        self.df_all[ph] = np.nan  # add column
                        self.df_all.loc[row, ph] = df[
                                                       ph].iloc[
                                                       -1] / m_tot * 100  # renormalise to 100 g total mass, only need last row

        print('done loading phases!')
        print(self.df_all.head())

    def read_melts_fo2(self, **kwargs):
        """ make csv of logfo2 - isothermal x-section for several pressures """

        self.read_melts_TP(**kwargs)

        fname = 'System_main_tbl.txt'
        for row, path in enumerate(self.output_p_paths):
            # load output csv from pMELTS
            output_file = path + fname
            df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                             dtype=np.float64).tail(1)
            # print('loaded\n', df.head())

            # append fo2
            self.df_all['logfo2'].iloc[row] = df['logfO2(absolute)'].iloc[-1]

            # append P and T
            if np.isnan(self.df_all['P(bar)'].iloc[row]):
                self.df_all['P(bar)'].iloc[row] = df['Pressure'].iloc[-1]
            elif self.df_all['P(bar)'].iloc[row] != df['Pressure'].iloc[-1]:
                error_str = 'Error: pressure mismatch!\nLoaded {} from {}\nHave {} in self.df_all'.format(
                    df['Pressure'].iloc[-1], output_file, self.df_all['P(bar)'].iloc[row])
                raise RuntimeError(error_str)
            if np.isnan(self.df_all['T(K)'].iloc[row]):
                self.df_all['T(K)'].iloc[row] = df['Temperature'].iloc[-1]  # K
            elif self.df_all['T(K)'].iloc[row] != df['Temperature'].iloc[-1]:
                error_str = 'Error: pressure mismatch!\nLoaded {} from {}\nHave {} in self.df_all'.format(
                    df['Temperature'].iloc[-1], output_file, self.df_all['T(K)'].iloc[row])
                raise RuntimeError(error_str)
        print('...done loading fO2!')

    def fo2_calc(self, compare_buffer=None, save=True, perplex_path=px.perplex_path_default, **kwargs):

        try:
            # run alphamelts
            self.run_alphamelts_all_p(**kwargs)

            # retrieve fo2
            self.read_melts_fo2()

            # retrieve phases
            self.read_melts_phases(which='mass')

            if compare_buffer == 'qfm':
                try:
                    logfo2_buffer = pf.read_qfm_os(self.df_all['T(K)'].to_numpy(), self.df_all['P(bar)'].to_numpy(),
                                                   verbose=False, perplex_path=perplex_path)
                    # print('logfo2_qfm', logfo2_buffer)
                    # print('logfo2 = QFM + ', logfo2 - logfo2_buffer, 'at', P, 'bar,', T, 'K')
                    self.df_all['logfo2_qfm'] = logfo2_buffer
                    self.df_all['delta_qfm'] = self.df_all.logfo2 - logfo2_buffer
                except NotImplementedError:
                    # probably didn't get to 1100 C
                    pass
            okay = True
        except FileNotFoundError as e:
            print(e)
            # alphamelts did not work at all for whatever reason
            okay = False

        # store mega df
        df_save = self.df_all.loc[:, ~self.df_all.columns.duplicated()].copy()
        if save:
            df_save.to_csv(self.output_path + self.name + '_results.csv', sep="\t")
            print('saved to', self.output_path + self.name + '_results.csv')
        return okay

    def find_common_T_final_from_results(self):
        # for already-ran melts data, get the last temperature that all pressures cooled to
        fname = 'Phase_mass_tbl.txt'
        T_min = self.T_final
        mass_melt_min = 0
        p_count = 0
        for ii, path in enumerate(self.output_p_paths):
            # load output csv from pMELTS
            output_file = path + fname
            df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                             dtype=np.float64).tail(1)
            p_count += 1
            T_last = df['Temperature'].item()
            # print('T_last', T_last, '(type =', type(T_last), ') at', path)
            if T_last >= T_min:
                T_min = T_last  # want highest value of min T
                mass_melt_min = df['liquid_0']  # also retrieve final melt mass fraction
        self.T_min = T_min
        self.mass_melt_min = mass_melt_min
        self.n_pressures = p_count



def init_from_results(name, output_parent_path=output_parent_default, alphamelts_path=alphamelts_path_default, T_final=1373, verbose=True, **kwargs):

    parts = name.split('_')
    star = parts[2]
    X_ferric = parts[4]
    wt_oxides = {}

    subfolders = [f.name for f in os.scandir(output_parent_path + name + '/') if f.is_dir()]

    if len(subfolders) > 0:
        pressures_of_interest = [float(s.replace(',', '.').replace('bar', '')) for s in subfolders]  # TODO parse directories

        # parse melts file
        try:
            melts_file_contents = open(output_parent_path + name + '/' + subfolders[0] + '/' + name + '.melts').readlines()
        except FileNotFoundError as e:
            print(e, 'skipping...')
            return None

        for line in melts_file_contents:
            if 'Initial Composition:' in line:
                parts = line.split()
                wt_oxides[parts[2]] = float(parts[3])

        dat = MeltsFugacityData(name=name, star=star, X_ferric=X_ferric, wt_oxides=wt_oxides, output_parent_path=output_parent_path,
                                alphamelts_path=alphamelts_path, pressures_of_interest=pressures_of_interest, T_final=T_final)

        # load results csv
        dat.df_all = pd.read_csv(dat.output_path + name + '_results.csv', sep='\t')
        if verbose:
            print('loaded df\n', dat.df_all.head())

        return dat
    else:
        print('no runs files at', name)
        return None


def fo2_from_hypatia(pressures_of_interest, n_sample=-1, core_efficiency=0.88, planet_kwargs={},
                     output_parent_path=output_parent_default, **kwargs):
    planet_kwargs.update({'core_efficiency': core_efficiency, 'solve_interior': False})
    pl_list = rw.planets_from_hypatia(n_sample=n_sample, plot_all=False,
                                      get_saturation=False, Tp=999,
                                      stopafter=None, output_parent_path=output_parent_path,
                                      **planet_kwargs, **kwargs)
    print('\nfinished generating compositions\n')
    bad = []
    for pl in pl_list:
        print(pl.wt_oxides)
        okay = fo2_from_oxides(pl.name, pressures_of_interest, pl=pl, output_parent_path=output_parent_path, **kwargs)
        if not okay:
            bad.append(pl.name)
    print('bad cases:', bad)
    return pl_list


def fo2_from_oxides(name, pressures_of_interest, pl=None,
                    core_efficiency=0.88, test_oxides=None, star=None, verbose=True, planet_kwargs={},
                    output_parent_path=output_parent_default, oxide_list=None, **kwargs):
    """ p_min and p_max in bar """

    if pl is None:
        # make planet object without solving interior structure, but getting bulk composition

        if star is not None and test_oxides is not None:
            raise NotImplementedError('Are you sure you want to input both test_oxides and star? Only need one.')
        planet_kwargs.update({'test_oxides': test_oxides, 'core_efficiency': core_efficiency})
        pl = rw.build_planet(name=name, star=star, get_saturation=False, solve_interior=False, verbose=verbose,
                             output_parent_path=output_parent_path, oxides=oxide_list, **planet_kwargs)

    try:
        tmp = pl.wt_oxides.items()
    except AttributeError as e:
        print(e)
        print(name, 'has no wt_oxides composition')
        return False
    print('\n\n\n\n----------------------------------------\nStarting fo2 calc for planet at', pl.star)
    dat = MeltsFugacityData(name=name, pressures_of_interest=pressures_of_interest,
                            wt_oxides=pl.wt_oxides, verbose=verbose, output_parent_path=output_parent_path,
                            **kwargs)
    okay = dat.fo2_calc(**kwargs)
    return okay


def common_Tmin(output_parent_path, **kwargs):
    names = [f.name for f in os.scandir(output_parent_path + '/') if f.is_dir()]
    df = pd.DataFrame(columns=['name', 'T_last', 'mass_melt_min', 'n_pressures'], index=range(len(names)))
    for row, sub in enumerate(names):
        dat = init_from_results(sub, output_parent_path=output_parent_path, **kwargs)
        if dat:  # existing runs
            dat.find_common_T_final_from_results()
            df.at[row, 'name'] = dat.name
            df.at[row, 'T_min'] = dat.T_min
            df.at[row, 'mass_melt_min'] = dat.mass_melt_min
            df.at[row, 'n_pressures'] = dat.n_pressures

    print(df)





