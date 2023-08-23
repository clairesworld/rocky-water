import numpy as np
import pandas as pd
import py.perplexdata as px
import py.bulk_composition as bulk
import os
import subprocess
import py.main as rw
import py.minfo2.perplexfugacitydata as pf
import py.parameters as p
from py.useful_and_bespoke import find_nearest_idx
from inspect import currentframe, getframeinfo
from pandas.core.common import SettingWithCopyError

pd.options.mode.chained_assignment = 'raise'

output_parent_default = '/home/claire/Works/min-fo2/alphamelts_output/'
alphamelts_path_default = '/home/claire/Works/alphamelts/'
output_parent_mlt_earth = '/home/claire/Works/min-fo2/alphamelts_output/earth-tea23/'

map_to_px_phase = {'olivine': 'Ol', 'orthopyroxene': 'Opx', 'clinopyroxene': 'Cpx', 'spinel': 'Sp', 'garnet': 'Gt',
                   'feldspar': 'Plag', 'quartz': 'q', 'coesite': 'coe'}


# solution_phases_default = ['olivine', 'spinel', 'garnet', 'orthopyroxene', 'clinopyroxene']  # opx and cpx automatic?

class MeltsFugacityData:

    def __init__(self, name='test', star=None,  # solution_phases=solution_phases_default,
                 output_parent_path=output_parent_default, wt_oxides=pf.wt_oxides_DMM, X_ferric=0.03, core_eff=None,
                 alphamelts_path=alphamelts_path_default,
                 T_final=None,
                 pressures_of_interest=None, # in bar
                 verbose=False,
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
        self.X_ferric = X_ferric
        self.core_eff = core_eff

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
        self.data = pd.DataFrame(columns=['P(bar)', 'T(K)', 'logfo2'], index=range(len(pressures_of_interest)))
        self.alphamelts_path = alphamelts_path
        self.logfo2_1GPa = None
        self.logfo2_4GPa = None
        self.delta_qfm_1GPa = None
        self.delta_qfm_4GPa = None

        # set main output path for this composition and subpaths for pressure runs
        self.output_path = output_parent_path + self.name + '/'
        self.output_p_paths = [self.output_path + str(int(pp)) + 'bar/' for pp in pressures_of_interest]
        self.pressures_of_interest = pressures_of_interest
        self.T_final = T_final


        if verbose:
            print('\ninitialising MeltsFugacityData data', name, 'with:')
            print('        T =', self.T_final, 'K')
            print('        p =', self.pressures_of_interest, 'bar')
            print('        X_ferric =', self.X_ferric)
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

    def print_comp(self):
        print('wt.% oxides\n-----------')
        for k in self.wt_oxides:
            print("{0:<7}".format(k), "{:5.2f}%".format(self.wt_oxides[k]))
        print('Mg/Si', self.mgsi)

        X_fer_back = bulk.test_ferric_ratio(self.wt_oxides['FeO'], self.wt_oxides['Fe2O3'])
        print('Fe3+/Fe', X_fer_back)

        # insert O2 like Perple_x - already gets printed
        test_wt_oxides = self.wt_oxides.copy()
        m_Fe2O3 = test_wt_oxides.pop('Fe2O3')
        test_wt_oxides['FeO'] = test_wt_oxides['FeO'] + m_Fe2O3 # add Fe2O3 back in
        wt_dict_o2 = bulk.insert_o2_from_wt_comp(test_wt_oxides, X_ferric=X_fer_back)

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

    def run_alphamelts_at_p(self, p_of_interest=None, output_p_path=None, suppress_output=False, clean=True,
                            verbose=True,
                            batch_text_fn=None, melts_text_fn=None,
                            env_file=None, dry_setup=False,
                            melts_kwargs=None, overwrite=False, verify_on_path=False, **kwargs):

        if output_p_path is None:
            output_p_path = self.output_path + str(int(p_of_interest)) + 'GPa/'
        if melts_kwargs is None:
            melts_kwargs = {}
        if env_file is None:
            env_file = self.alphamelts_path + 'examples/alphamelts_env_isobaric.txt'

        # check if run already exists
        if (not os.path.isfile(output_p_path + 'System_main_tbl.txt')) or overwrite:

            # os.chdir(self.alphamelts_path)
            os.chdir(output_p_path)

            if suppress_output:
                stderr, stdout = subprocess.DEVNULL, subprocess.DEVNULL
            else:
                stderr, stdout = None, None

            # melts_file = output_p_path + self.name + '.melts'
            # batch_file = output_p_path + self.name + '.in'
            # melts_file = self.name + '.melts'
            # batch_file = self.name + '.in'
            melts_file = 'pl.melts'
            batch_file = 'batch.in'

            # create batch and melts files
            with open(batch_file, 'w') as file:
                file.write(batch_text_fn(melts_file))
            with open(melts_file, 'w') as file:
                file.write(melts_text_fn(p_of_interest, **melts_kwargs))

            # run alphamelts
            if not dry_setup:
                # make sure installation is on path
                if verify_on_path:
                    subprocess.call(['sh', './verify_path.sh'])  # should be in alphamelts directory

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
                                'ERROR: alphamelts did not complete, try running again with suppress_output=False:', self.name)

            # return to original dir
            os.chdir(os.path.dirname(os.path.abspath(__file__)))
        else:
            print('Run', self.name, 'at p =', p_of_interest,
                  'bar already exists! To execute, delete files or set overwrite=True\n---------------------------------------------')

    def read_melts_TP(self, T_of_interest=1373.15, reload_TP=False, verbose=False, **kwargs):

        # check if already loaded
        if (not self.data['P(bar)'].isnull().values.any()) or (not self.data['T(K)'].isnull().values.any()):
            if not reload_TP:
                if verbose:
                    print('                already loaded T,P from alphaMELTS! to overwrite, use reload_TP=True')
                return True

        fname = 'System_main_tbl.txt'
        for row, path in enumerate(self.output_p_paths):
            # load output csv from pMELTS
            output_file = path + fname
            df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                             dtype=np.float64)

            # find idx of T of interest
            try:
                idx = df.loc[df['Temperature'] == T_of_interest].index[0]
                new_T = df.loc[idx, 'Temperature']
                new_P = int(df.loc[idx, 'Pressure'])  # floor to integer bar to avoid fp errors
            except IndexError as e:
                print('    System_main_tbl: T_of_interest', T_of_interest, 'K not found at',  self.pressures_of_interest[row], 'bar ( min:', df.Temperature.iloc[-1], ')', self.name)
                new_T = np.nan
                new_P = np.nan

            # append P and T
            if np.isnan(self.data.loc[row, 'P(bar)']):
                self.data.loc[row, 'P(bar)'] = new_P
            elif int(self.data.loc[row, 'P(bar)']) != new_P:
                error_str = 'Error: pressure mismatch!\nLoaded {} from {}\nHave {} in self.data'.format(
                    new_P, output_file, self.data.loc[row, 'P(bar)'])
                raise RuntimeError(error_str)
            if np.isnan(self.data.loc[row, 'T(K)']):
                self.data.loc[row, 'T(K)'] = new_T  # K
            elif self.data.loc[row, 'T(K)'] != new_T:
                error_str = 'Error: pressure mismatch!\nLoaded {} from {}\nHave {} in self.data'.format(
                    new_T, output_file, self.data.loc[row, 'T(K)'])
                raise RuntimeError(error_str)
        return True

    def read_melts_phases(self, T_of_interest=1373.15, which='mass', verbose=False, **kwargs):
        import re
        """ make csv of phase proportions in same format as perple_x - isothermal x-section for several pressures """

        if which == 'mass':
            fname = 'Phase_mass_tbl.txt'
        elif which == 'volume':
            fname = 'Phase_vol_tbl.txt'

        okay = self.read_melts_TP(T_of_interest=T_of_interest, verbose=verbose, **kwargs)
        if not okay:
            raise Exception(self.name, 'read_melts_TP crashed')
        if okay:
            for row, path in enumerate(self.output_p_paths):
                try:
                    # load output csv from pMELTS
                    output_file = path + fname
                    df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                                     dtype=np.float64)
                    # print('loaded\n', df.head())

                    try:
                        # find idx of T of interest
                        idx = df.loc[df['Temperature'] == T_of_interest].index[0]
                    except IndexError:
                        if verbose:
                            print(T_of_interest, 'K not found at', self.pressures_of_interest[row], 'bar')
                        # T_of_interest not found at pressure_of_interest
                        continue  # skip to next pressure

                    # append phases to self df
                    m_tot = df['mass'].loc[idx]
                    for ph in [col for col in df.columns if col not in ['Pressure', 'Temperature', 'mass', 'volume']]:
                        if ph not in ['liquid_0', 'alloy-liquid_0']:  # ignoring melt phases
                            ph2 = re.sub('_0', '', ph)
                            try:
                                label = 'X_' + map_to_px_phase[ph2]
                            except KeyError as e:
                                if ph2.endswith('_1'):
                                    if df.loc[idx, ph] != 0:
                                        # add to existing column
                                        ph3 = re.sub('_1', '', ph)
                                        label = 'X_' + map_to_px_phase[ph3]
                                        addition = self.data.loc[row, label] + df.loc[idx, ph]
                                        self.data.loc[row, label] = addition / m_tot * 100  # renormalise to 100 g total mass
                                        # raise Exception(self.name, 'repeated phase', ph2, '(nonzero):', path)
                                    else:
                                        continue
                                else:
                                    print('missing', ph2, 'in map_to_px_phase dictionary:', path)
                                    if self.pressures_of_interest[row] not in (10e3, 40e3):
                                        continue  # don't care
                                    else:
                                        raise e
                            try:
                                self.data.loc[row, label] = df.loc[idx, ph] / m_tot * 100  # renormalise to 100 g total mass
                            except KeyError:
                                self.data[label] = np.nan  # add column
                                # print('cols', self.data.columns)
                                self.data.loc[row, label] = df.loc[idx, ph] / m_tot * 100  # renormalise to 100 g total mass
                except SettingWithCopyError:
                    print('handling..')
                    frameinfo = getframeinfo(currentframe())
                    print(frameinfo.lineno)

            if verbose:
                print('             ...done loading phases!')
            # print(self.data.head())

        # also load Fe2O3 content
        # print('(checking indices of self.data')
        # print(self.data)
        for row, pp in enumerate(self.pressures_of_interest):
            d_Fe3 = self.read_phase_main_components(pp, T_of_interest, component='Fe2O3', phases=map_to_px_phase.keys(),
                                                    absolute_abundance=False, verbose=False, p_index=row)
            # print('row', row, 'pp', pp)
            if d_Fe3:
                for key in d_Fe3:
                    label = 'X_Fe3_' + key
                    if d_Fe3[key] > 0:
                        try:
                            self.data.loc[row, label] = d_Fe3[key]  # Fe3 composition of phase
                        except KeyError:
                            self.data[label] = np.nan  # add column
                            self.data.loc[row, label] = d_Fe3[key]  # Fe3 composition of phase

    def read_melts_fo2(self, T_of_interest=1373.15, verbose=False, **kwargs):
        """ make csv of logfo2 - isothermal x-section for several pressures """

        okay = self.read_melts_TP(T_of_interest=T_of_interest, verbose=verbose, **kwargs)
        if not okay:
            raise Exception(self.name, 'read_melts_TP crashed in read_melts_fo2()')
        if okay:
            fname = 'System_main_tbl.txt'
            for row, path in enumerate(self.output_p_paths):
                try:
                    # load output csv from pMELTS
                    output_file = path + fname
                    df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                                     dtype=np.float64)
                    # print('loaded\n', df.head())

                    # find idx of T of interest
                    try:
                        idx = df.loc[df['Temperature'] == T_of_interest].index[0]

                        # append fo2
                        self.data.loc[row, 'logfo2'] = df['logfO2(absolute)'].loc[idx]
                        # print('idx', idx, 'T', df.Temperature.loc[idx], self.data['T(K)'].iloc[row])

                    except IndexError:
                        # T_of_interest not found
                        self.data.loc[row, 'logfo2'] = np.nan

                except SettingWithCopyError:
                    print('handling..')
                    frameinfo = getframeinfo(currentframe())
                    print('line', frameinfo.lineno)

            if verbose:
                print('             ...done loading fO2!')

    def fo2_calc(self, compare_buffer=None, save=True, perplex_path=px.perplex_path_default, run_alphamelts=True,
                 verbose=False, T_of_interest=1373.15, **kwargs):

        try:
            # run alphamelts -- note will check if run exists at run_alphamelts_at_p()

            if run_alphamelts:
                self.run_alphamelts_all_p(T_of_interest=T_of_interest, **kwargs)

            if ('dry_setup' not in kwargs) or not (kwargs['dry_setup']):
                try:
                    # retrieve fo2
                    self.read_melts_fo2(T_of_interest=T_of_interest, verbose=verbose, **kwargs)

                    # retrieve phases
                    self.read_melts_phases(which='mass', T_of_interest=T_of_interest, verbose=verbose, **kwargs)

                    # if verbose:
                    #     print(self.name, 'finished reading')
                    #     print(self.data.head())
                except SettingWithCopyError:
                    print('handling..')
                    frameinfo = getframeinfo(currentframe())
                    print(frameinfo.lineno)

                if compare_buffer == 'qfm':
                    try:
                        print('getting QFM')
                        # use isothermal T but make sure this is ok
                        tmp = self.data.copy(deep=True)['T(K)']
                        tmp.dropna(inplace=True)
                        tmpa = tmp.to_numpy()
                        try:
                            if (tmpa[0] == tmpa).all():
                                T0 = tmpa[0]
                            else:
                                print(self.name, '\n', tmp)
                                raise NotImplementedError('ERROR: multiple values of T in melts results?')
                            logfo2_buffer = pf.read_qfm_os(T0, self.data['P(bar)'].to_numpy(),
                                                           verbose=False, perplex_path=perplex_path)
                            self.data['logfo2_qfm'] = logfo2_buffer
                            self.data['delta_qfm'] = self.data.logfo2 - logfo2_buffer
                            # print('added qfm\n', self.data['delta_qfm'])
                        except IndexError:
                            # tmpa is size 0 - all nan?
                            print('no valid temperatures for QFM interpolation (all NaN?)', self.name)
                            print(self.data['T(K)'])
                            self.data['logfo2_qfm'] = np.nan
                            self.data['delta_qfm'] = np.nan

                    except NotImplementedError as e:
                        # probably didn't get to 1100 C <-- idk about this anymore 5/03
                        print(self.name, e)
                        print('T', self.data['T(K)'].to_numpy(), '\nP', self.data['P(bar)'].to_numpy())
                        raise e
                        # pass
            else:
                save = False  # don't save empty results csv for dry run with no calculations
            okay = True
        except FileNotFoundError as e:
            if ('dry_setup' not in kwargs) or not (kwargs['dry_setup']):
                print(e)
                # alphamelts did not work at all for whatever reason
                okay = False
            else:
                okay = True
                print('Dry run, ignoring FileNotFoundError:', e)

        # store mega df
        df_save = self.data.loc[:, ~self.data.columns.duplicated()].copy()
        # print('df saved\n', df_save.head())
        if save:
            df_save.to_csv(self.output_path + self.name + '_results' + str(int(T_of_interest)) + '.csv', sep="\t")
            print('>>> saved to', self.output_path + self.name + '_results' + str(int(T_of_interest)) + '.csv')
        return okay

    def find_common_T_final_from_results(self, include_p=None, **kwargs):
        # for already-ran melts data, get the last temperature that *all* pressures cooled to
        fname = 'Phase_mass_tbl.txt'
        T_min = self.T_final
        mass_melt_min = 0
        p_count = 0
        for ii, path in enumerate(self.output_p_paths):
            if (not include_p) or (include_p and (self.pressures_of_interest[ii] in include_p)):
                # load output csv from pMELTS
                output_file = path + fname

                try:
                    df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                                     dtype=np.float64).tail(1)
                except FileNotFoundError:
                    # this run didn't complete
                    print(self.name, 'did not complete at', self.pressures_of_interest[ii], 'bar: no Phase_mass_tbl.txt found in', path)
                    self.T_min = np.nan
                    self.mass_melt_min = np.nan
                    self.n_pressures = 0
                    return None

                p_count += 1
                T_last = df['Temperature'].item()
                # print('T_last', T_last, '(type =', type(T_last), ') at', path)
                if T_last >= T_min:
                    T_min = T_last  # want highest value of min T
                    mass_melt_min = df['liquid_0'].item()  # also retrieve final melt mass fraction
        self.T_min = T_min
        self.mass_melt_min = mass_melt_min
        self.n_pressures = p_count

    def read_fo2_results(self, T_of_interest=1373.15, verbose=True, **kwargs):

        try:
            self.data = pd.read_csv(self.output_path + self.name + '_results' + str(int(T_of_interest)) + '.csv', sep='\t')
            if verbose:
                print('loaded df\n', self.data.head())
        except FileNotFoundError:
            print('results csv file not found for', int(T_of_interest), 'K, skipping', self.name)
            return None

        self.logfo2 = self.data['logfo2'].to_numpy()
        try:
            self.delta_qfm = self.data['delta_qfm'].to_numpy()
        except KeyError:
            pass

        pressure = self.data['P(bar)'].to_numpy() * 1e-4
        temperature = self.data['T(K)'].to_numpy()
        if (not hasattr(self, 'pressure')) or self.pressure is None:
            self.pressure = pressure
        elif (self.pressure != pressure).any():
            print('pressure', pressure)
            print('self.pressure', self.pressure)
            raise Exception('attempting to read in mismatching pressure data!')
        if (not hasattr(self, 'temperature')) or self.temperature is None:
            self.temperature = temperature
        elif (self.temperature != temperature).any():
            print('temperature', temperature)
            print('self.temperature', self.temperature)
            raise Exception('attempting to read in mismatching temperature data!')

        # save 1 GPa i.e. spinel stability field
        try:
            idx = find_nearest_idx(pressure, 1, ignorenan=True)
            self.logfo2_1GPa = self.logfo2[idx]
            self.delta_qfm_1GPa = self.delta_qfm[idx]
        except ValueError:
            # all nans
            self.logfo2_1GPa = np.nan
            idx = None
        except (KeyError, AttributeError):
            pass
        if np.isnan(self.logfo2_1GPa) and verbose:
            print('\nwarning:', self.name, 'logfo2_1GPa', self.logfo2_1GPa)
            print('idx', idx)
            print('pressure', pressure)
            print('logfo2', self.logfo2)
            print('self.data\n', self.data.head(10))

        # save 4 GPa i.e. garnet stability field
        try:
            idx = find_nearest_idx(pressure, 4, ignorenan=True)
            self.logfo2_4GPa = self.logfo2[idx]
            self.delta_qfm_4GPa = self.delta_qfm[idx]
            idx = None
        except ValueError:
            # all nans
            self.logfo2_4GPa = np.nan
        except (KeyError, AttributeError):
            pass

        if np.isnan(self.logfo2_4GPa) and verbose:
            print('\nwarning:', self.name, 'logfo2_4GPa', self.logfo2_4GPa)
            print('idx', idx)
            print('pressure', pressure)
            print('logfo2', self.logfo2)
            print('self.data\n', self.data.tail(10))



    def read_phase_main(self, phase, p_of_interest, T_of_interest, verbose=False):
        filename = self.output_path + str(int(p_of_interest)) + 'bar/Phase_main_tbl.txt'
        # print('filename l. 496', filename)
        # print('    p_of_interest', p_of_interest)
        tmp = []
        with open(filename) as file:
            start = False
            for line in file:
                line = line.rstrip()  # remove trailing whitespace
                if line.startswith(phase):
                    start = True
                if start:
                    if line.strip() == "":
                        break  # finished once you get to blank line
                    # find T_of_interest
                    bits = line.split()
                    if (bits[1] == 'Temperature') or (bits[1] == str(T_of_interest)):
                        tmp.append(line)
        # print('tmp', tmp)
        if not start:
            # never found this phase
            if verbose:
                print(phase, 'not found')
            return None
        if len(tmp) < 2:
            if verbose:
                # already print out missing T in read_melts_TP()
                print(T_of_interest, 'K not found at', p_of_interest, 'bar in Phase_main_tbl.txt: ', self.name)
            return None

        # write temp file (not sure how to make df otherwise)
        with open('tmp.txt', 'w') as f:
            for line in tmp:
                f.write(f"{line}\n")

        df = pd.read_csv('tmp.txt', skiprows=0, index_col=None, sep=r"\s+")
        # print(df.head())

        os.remove('tmp.txt')
        return df

    # def read_phase_comp(self, p_of_interest, T_of_interest, component='Fe2O3', phases=map_to_px_phase.keys(),
    #                     absolute_abundance=True, verbose=False, p_index=None):
    #     """
    #     reads component composition from Phase_main_table (wt%), stores in results df
    #
    #     Parameters
    #     ----------
    #     p_index :
    #     absolute_abundance :
    #     phases :
    #     T_of_interest :
    #     component :
    #     p_of_interest : bar
    #     verbose :
    #
    #     Returns
    #     -------
    #
    #     """
    def read_phase_main_components(self, p_of_interest, T_of_interest, component='Fe2O3', phases=map_to_px_phase.keys(),
                        absolute_abundance=True, verbose=False, p_index=None):
        """
        reads phase comp from Phase_main>table, returns dict of each phases' compositon of that component

        Parameters
        ----------
        p_index :
        absolute_abundance :
        phases :
        T_of_interest :
        component :
        p_of_interest : bar
        verbose :

        Returns
        -------

        """
        if not hasattr(self, 'data'):
            self.data = pd.read_csv(self.output_path + self.name + '_results' + str(int(T_of_interest)) + '.csv', sep='\t')

        if p_index:
            idx = p_index
        else:
            try:
                idx = self.pressures_of_interest.index(p_of_interest)
            except ValueError as e:
                # pressure not found
                print(p_of_interest, 'bar not found')
                raise e
            except AttributeError:
                # array?
                idx = np.abs(self.pressures_of_interest - p_of_interest).argmin()

        wt_pt_dict = {map_to_px_phase[ph]: None for ph in phases}
        try:
            for phase in phases:
                df = self.read_phase_main(phase, p_of_interest, T_of_interest, verbose=verbose)

                if df is None:  # phase not found
                    wt_pt_dict[map_to_px_phase[phase]] = np.nan
                else:
                    if verbose:
                        print(phase, 'loaded df\n', df.head())

                    if absolute_abundance:
                        # normalise to total quantity of phase
                        try:
                            mass_ph = self.data.loc[idx, 'X_' + map_to_px_phase[phase]] / 100  # these are wt%
                            print('mass', phase, mass_ph)
                        except KeyError as e:
                            # this is probably because T_of_interest not found - won't fill in _results.csv
                            # for 1M_88Ceff_HIP84856_999K_3,0fer, other weird issue
                            print(self.name, e)
                            print('old', '\n', self.data.head())

                            self.read_melts_phases(T_of_interest=T_of_interest, which='mass', verbose=True, reload_TP=True)
                            print('new', '\n', self.data.head())
                            return None
                    else:
                        mass_ph = 1

                    try:
                        wt_pt_dict[map_to_px_phase[phase]] = df.loc[0, component] * mass_ph
                    except KeyError:
                        if verbose:
                            print(component, 'not found in', phase)
                        wt_pt_dict[map_to_px_phase[phase]] = np.nan
        except FileNotFoundError as e:
            print(e)
            print('...file not found at', p_of_interest, 'bar, ', int(T_of_interest), 'K! skipping:', self.name, '\n')
            return None
        return wt_pt_dict
    
    def get_phase_composition_dict(self, p_of_interest=None, T_of_interest=None, component='Fe2O3',
                                   phases=map_to_px_phase.keys(),
                                   to_absolute_abundance=True, verbose=False):
        """ 
        read results csv to create phase composition dict (think only for use with ternay diagrams so far)
        phase component masses as stored in results.csv are relative to phase, not system!

        
        Parameters
        ----------
        p_of_interest :
        T_of_interest :
        phases :
        component :
        to_absolute_abundance : 

        Returns
        -------

        """

        # load results.csv dataframe
        if not hasattr(self, 'data'):
            self.data = pd.read_csv(self.output_path + self.name + '_results' + str(int(T_of_interest)) + '.csv', sep='\t')

        # find index in df of pressure of interest
        try:
            idx = self.pressures_of_interest.index(p_of_interest)
        except ValueError as e:
            # pressure not found
            # print(p_of_interest, 'bar not found in', self.pressures_of_interest, ':', self.name, '(probably melts crashed before cooling to', T_of_interest, 'K)')
            return None

        # initialise d
        try:
            wt_pt_dict = {map_to_px_phase[ph]: None for ph in phases}
        except KeyError:
            wt_pt_dict = {ph: None for ph in phases}
        for phase in phases:
            try:
                phase0 = map_to_px_phase[phase]
            except KeyError:
                phase0 = phase

            if to_absolute_abundance:
                # get phase mass as well so you can normalise to total quantity of phase
                try:
                    mass_ph = self.data.loc[idx, 'X_' + map_to_px_phase[phase]] / 100  # these are wt%
                    if verbose:
                        print('mass', phase, mass_ph)
                except KeyError as e:
                    # this is probably because T_of_interest not found - won't fill in _results.csv
                    # for 1M_88Ceff_HIP84856_999K_3,0fer, other weird issue
                    # HIP 80680: no spinel for some reason (but garnet)

                    if verbose:
                        print('\n', self.name, ': column', e, 'not found (p =', p_of_interest, ')')
                        print('           ', list(self.data.columns)[3:])
                    mass_ph = 0
            else:
                mass_ph = 1

            # extract column of phase component mass
            if component == 'Fe2O3':
                col = 'X_Fe3_' + phase0
            else:
                raise NotImplementedError('Not Implemented: component', component, 'in phase', phase0, 'not processed in results csv')
            try:
                wt_pt_dict[phase0] = self.data.loc[idx, col] * mass_ph
            except KeyError as e:
                if verbose and not to_absolute_abundance:
                    print('\n', self.name, ': column', e, 'not found (p =', p_of_interest, ')')
                    print('           ', list(self.data.columns)[3:])
                wt_pt_dict[phase0] = 0

        return wt_pt_dict


def init_from_results(name, output_parent_path=output_parent_default, alphamelts_path=alphamelts_path_default,
                      T_final=1373, load_results_csv=True, verbose=False, X_ferric=None, core_eff=None, **kwargs):
    import re
    parts = name.split('_')
    try:
        star = parts[2]
        X_ferric = float(re.findall('\d+\.\d+', parts[4].replace(',', '.'))[0] ) * 1e-2
        core_eff = float(re.findall('\d+\.\d+', parts[1].replace(',', '.'))[0] ) * 1e-2
        print('core eff', core_eff)
    except IndexError:
        # non conventional name
        star = None
    wt_oxides = {}

    subfolders = [f.name for f in os.scandir(output_parent_path + name + '/') if f.is_dir()]

    # check if solution
    if (len(subfolders) > 0) and ((os.path.isfile(output_parent_path + name + '/' + subfolders[0] + '/pl.melts')) or ((os.path.isfile(output_parent_path + name + '/' + subfolders[0] + '/' + name + '.melts')))):

        # load pressures
        pressures_of_interest = [float(s.replace(',', '.').replace('bar', '')) for s in subfolders]

        # sort ascending (alphabetical order might be different in os.scandir
        pressures_of_interest.sort()

        # print('initiating pressures of interest', pressures_of_interest)

        # parse melts file
        try:
            melts_file_contents = open(output_parent_path + name + '/' + subfolders[0] + '/pl.melts').readlines()
        except FileNotFoundError as e:
            try:
                melts_file_contents = open(output_parent_path + name + '/' + subfolders[0] + '/' + name + '.melts').readlines()
            except FileNotFoundError:
                print(e, 'pl.melts not found in', subfolders[0], 'skipping', name)
                return None

        for line in melts_file_contents:
            if 'Initial Composition:' in line:
                parts = line.split()
                wt_oxides[parts[2]] = float(parts[3])

        dat = MeltsFugacityData(name=name, star=star, X_ferric=X_ferric, core_eff=core_eff, wt_oxides=wt_oxides, output_parent_path=output_parent_path,
                                alphamelts_path=alphamelts_path, pressures_of_interest=pressures_of_interest, T_final=T_final, verbose=verbose)

        # load results csv
        if load_results_csv:
            try:
                # print('T_final in init_from_results()', T_final)
                dat.data = pd.read_csv(dat.output_path + name + '_results' + str(int(T_final)) + '.csv', sep='\t')
                dat.read_fo2_results(verbose=verbose, T_of_interest=T_final, **kwargs)

                # make sure results actually contains data (file may have been created with nans)
                if dat.data['P(bar)'].isnull().all() or dat.data['T(K)'].isnull().all():
                    if verbose:
                        print(name + '_results' + str(int(T_final)) + '.csv file is empty/nan (calculation failed?), skipped...')
                    return None

                # update (e.g. if local folders with only 10000 bar copied)
                dat.pressures_of_interest = list(dat.data['P(bar)'].to_numpy())

                # if verbose:
                #     print('loaded df\n', dat.data.head())

            except FileNotFoundError:
                if verbose:
                    print(name + '_results' + str(int(T_final)) + '.csv file not found, skipped...')
                return None

        return dat
    else:
        print('no melts files at', output_parent_path + name)
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


def common_Tmin(output_parent_path, store=True, **kwargs):
    names = [f.name for f in os.scandir(output_parent_path) if f.is_dir()]
    df = pd.DataFrame(columns=['name', 'T_min'
                                       , 'mass_melt_min', 'n_pressures'], index=range(len(names)))
    for row, sub in enumerate(names):
        dat = init_from_results(sub, output_parent_path=output_parent_path, verbose=False, **kwargs)
        if dat:  # existing runs
            dat.find_common_T_final_from_results(**kwargs)
            df.at[row, 'name'] = dat.name
            df.at[row, 'T_min'] = dat.T_min
            df.at[row, 'mass_melt_min'] = dat.mass_melt_min
            df.at[row, 'n_pressures'] = dat.n_pressures

    print('\ncompletion_analysis.csv\n', df)
    if store:
        df.to_csv(output_parent_path + 'completion_analysis.csv', sep="\t")


def fo2_from_local(output_parent_path, num=-1, restart=None, names=None, T_of_interest=None, **kwargs):
    """for existing melts runs (e.g. done remotely), do the same fo2 analysis"""

    # get all runs in directory
    subfolders = [f.name for f in os.scandir(output_parent_path) if f.is_dir()]
    bad = []

    if restart:
        # print('subfolders', subfolders)
        # star = restart.split('_')[2]
        # print('star', star)
        # print('[x.split("_")[2] for x in subfolders]', [x.split('_')[2] for x in subfolders])
        idx0 = [x.split('_')[2] for x in subfolders].index(restart)  # find directory matching star name only
        # idx0 = subfolders.index(restart)
    else:
        idx0 = 0
    for name in subfolders[idx0:num]:
        if ((names is not None) and (name in names)) or (names is None):
            print('-------------------------------\nStarting fo2 processing:', name)
            dat = init_from_results(name, output_parent_path=output_parent_path, load_results_csv=False,
                                    T_final=T_of_interest, **kwargs)
            if dat is not None:
                okay = dat.fo2_calc(run_alphamelts=False, T_of_interest=T_of_interest, **kwargs)
            if (dat is None) or (not okay):
                bad.append(name)
            print('\n\n\n\n\n\n')
    return bad


def create_isothermal_csv(output_parent_path, T, P, Xfer, coreeff, verbose=True, exclude_silica=True,
                          dropnan=True, **kwargs):
    """ create summary of this case , T in K, P in GPa """
    pd.set_option("display.max_columns", 99)

    # get all runs in directory
    subfolders = [f.name for f in os.scandir(output_parent_path) if f.is_dir()]

    n = len(subfolders)
    dfout = pd.DataFrame(index=range(n), columns=['star', 'T(K)', 'P(GPa)',
                                                  'logfO2', 'deltaQFM', 'Fe3/Fe', 'Fe_c/Fe_T',
                                                  'X_Ol', 'X_Opx', 'X_Cpx', 'X_Sp', 'X_Gt',
                                                  'X_Fe3_Opx', 'X_Fe3_Cpx', 'X_Fe3_Sp', 'X_Fe3_Gt',
                                                  'Mg/Si', 'Fe/Si', 'Al/Si', 'Ca/Si'])

    # add constant cols
    dfout['T(K)'] = T
    dfout['P(GPa)'] = P
    dfout['Fe3/Fe'] = Xfer
    dfout['Fe_c/Fe_T'] = coreeff

    for row, name in enumerate(subfolders):
        # ensure name is there
        parts = name.split('_')
        star = parts[2]
        dfout.at[row, 'star'] = star

        dat = init_from_results(name, output_parent_path=output_parent_path, load_results_csv=True, **kwargs)
        if dat is not None:
            if exclude_silica:
                if 'X_q' in dat.data.columns:  # don't print quartz saturation cases
                    print('dropping case with quartz:', dat.name)
                    continue

            # get run data at T, p
            df = dat.data.loc[(dat.data['T(K)'] == T) & (dat.data['P(bar)'] == P * 1e4)]
            if len(df.index) > 0:  # maybe melts crashed or didn't run for some reason
                # if verbose:
                #     print('retrieved\n', df)

                # add to output csv
                dfout.at[row, 'logfO2'] = df['logfo2'].values[0]
                if not np.isnan(df['logfo2'].values[0]):
                    dfout.at[row, 'deltaQFM'] = df['delta_qfm'].values[0]
                else:
                    dfout.at[row, 'deltaQFM'] = np.nan

                for s in ['X_Ol', 'X_Opx', 'X_Cpx', 'X_Sp', 'X_Gt', 'X_Fe3_Opx', 'X_Fe3_Cpx', 'X_Fe3_Sp', 'X_Fe3_Gt']:
                    try:
                        dfout.at[row, s] = df[s].values[0]
                    except KeyError as e:
                        # if verbose:
                        #     print(e)
                        pass # this phase not stable at T,p

                for s in ['Mg/Si', 'Fe/Si', 'Al/Si', 'Ca/Si']:
                    dfout.at[row, s] = bulk.get_element_ratio(s, dat.wt_oxides)

    if dropnan:
        dfout.dropna(subset=['logfO2'], inplace=True)

    dfout.to_csv(output_parent_path + 'summary_melts_' + str(T).replace('.', ',') + 'K_' + str(P).replace('.', ',') + 'GPa_' + str(Xfer*100).replace('.', ',') + 'fer_' + str(coreeff*100).replace('.', ',') + 'core' + '.csv', sep="\t", na_rep='NaN')
    if verbose:
        print(dfout.head(10))
    pd.reset_option("max_columns")
    print('Done!')
    return None









