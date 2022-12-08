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

output_parent_default = '/home/claire/Works/min-fo2/alphamelts_output/'
alphamelts_path_default = '/home/claire/Works/alphamelts/'

map_to_px_phase = {'olivine': 'Ol', 'orthopyroxene': 'Opx', 'clinopyroxene': 'Cpx', 'spinel': 'Sp', 'garnet': 'Gt',
                   'feldspar': 'Plag', 'quartz': 'q', 'coesite': 'coe'}


# solution_phases_default = ['olivine', 'spinel', 'garnet', 'orthopyroxene', 'clinopyroxene']  # opx and cpx automatic?

class MeltsFugacityData:

    def __init__(self, name='test', star=None,  # solution_phases=solution_phases_default,
                 output_parent_path=output_parent_default, wt_oxides=pf.wt_oxides_DMM, X_ferric=0.03,
                 alphamelts_path=alphamelts_path_default,
                 T_final=None, pressures_of_interest=None, verbose=False,
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

        # set main output path for this composition and subpaths for pressure runs
        self.output_path = output_parent_path + self.name + '/'
        self.output_p_paths = [self.output_path + str(pp).replace('.', ',') + 'bar/' for pp in pressures_of_interest]
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
            output_p_path = self.output_path + str(p_of_interest).replace('.', ',') + 'GPa/'
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
                                'ERROR: alphamelts did not complete, try running again with suppress_output=False')

            # return to original dir
            os.chdir(os.path.dirname(os.path.abspath(__file__)))
        else:
            print('Run', self.name, 'at p =', p_of_interest,
                  'bar already exists! To execute, delete files or set overwrite=True\n---------------------------------------------')

    def read_melts_TP(self, T_of_interest=1373.15, reload_TP=False, **kwargs):

        # check if already loaded
        if (not self.data['P(bar)'].isnull().values.any()) or (not self.data['T(K)'].isnull().values.any()):
            if not reload_TP:
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
            except IndexError as e:
                print('             T_of_interest not found, min:', df.Temperature.iloc[-1], '--- in System_main', self.name)
                return False

            # append P and T
            if np.isnan(self.data.loc[row, 'P(bar)']):
                self.data.loc[row, 'P(bar)'] = df['Pressure'].loc[idx]
            elif self.data.loc[row, 'P(bar)'] != df['Pressure'].loc[idx]:
                error_str = 'Error: pressure mismatch!\nLoaded {} from {}\nHave {} in self.data'.format(
                    df['Pressure'].loc[idx], output_file, self.data['P(bar)'].iloc[row])
                raise RuntimeError(error_str)
            if np.isnan(self.data['T(K)'].iloc[row]):
                self.data.loc[row, 'T(K)'] = df['Temperature'].iloc[idx]  # K
            elif self.data.loc[row, 'T(K)'] != df['Temperature'].iloc[idx]:
                error_str = 'Error: pressure mismatch!\nLoaded {} from {}\nHave {} in self.data'.format(
                    df['Temperature'].loc[idx], output_file, self.data['T(K)'].iloc[row])
                raise RuntimeError(error_str)
        return True

    def read_melts_phases(self, T_of_interest=1373.15, which='mass', verbose=False, **kwargs):
        import re
        """ make csv of phase proportions in same format as perple_x - isothermal x-section for several pressures """

        if which == 'mass':
            fname = 'Phase_mass_tbl.txt'
        elif which == 'volume':
            fname = 'Phase_vol_tbl.txt'

        okay = self.read_melts_TP(T_of_interest=T_of_interest, **kwargs)
        if okay:
            for row, path in enumerate(self.output_p_paths):
                # load output csv from pMELTS
                output_file = path + fname
                df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                                 dtype=np.float64)
                # print('loaded\n', df.head())

                # find idx of T of interest
                idx = df.loc[df['Temperature'] == T_of_interest].index[0]

                # append phases to self df
                m_tot = df['mass'].loc[idx]
                for ph in [col for col in df.columns if col not in ['Pressure', 'Temperature', 'mass', 'volume']]:
                    if ph != 'liquid_0':
                        ph2 = re.sub('_0', '', ph)
                        try:
                            label = 'X_' + map_to_px_phase[ph2]
                        except KeyError as e:
                            if ph2.endswith('_1'):
                                if df[ph].loc[idx] != 0:
                                    raise Exception(self.name, 'repeated phase', ph2, '(nonzero)')
                                else:
                                    continue
                            print('missing', ph2, 'in map_to_px_phase dictionary', self.name)
                            raise e
                        try:
                            self.data[label].iloc[row] = df[ph].loc[
                                                           idx] / m_tot * 100  # renormalise to 100 g total mass
                        except KeyError:
                            self.data[label] = np.nan  # add column
                            print('cols', self.data.columns)
                            self.data[label].iloc[row] = df[ph].loc[
                                                           idx] / m_tot * 100  # renormalise to 100 g total mass

            if verbose:
                print('             ...done loading phases!')
            # print(self.data.head())

        # also load Fe2O3 content

    def read_melts_fo2(self, T_of_interest=1373.15, verbose=False, **kwargs):
        """ make csv of logfo2 - isothermal x-section for several pressures """

        okay = self.read_melts_TP(T_of_interest=T_of_interest, **kwargs)
        if okay:
            fname = 'System_main_tbl.txt'
            for row, path in enumerate(self.output_p_paths):
                # load output csv from pMELTS
                output_file = path + fname
                df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                                 dtype=np.float64)
                # print('loaded\n', df.head())

                # find idx of T of interest
                idx = df.loc[df['Temperature'] == T_of_interest].index[0]

                # append fo2
                self.data['logfo2'].iloc[row] = df['logfO2(absolute)'].loc[idx]
                # print('idx', idx, 'T', df.Temperature.loc[idx], self.data['T(K)'].iloc[row])

            if verbose:
                print('             ...done loading fO2!')

    def fo2_calc(self, compare_buffer=None, save=True, perplex_path=px.perplex_path_default, run_alphamelts=True,
                 verbose=False, **kwargs):

        try:
            # run alphamelts -- note will check if run exists at run_alphamelts_at_p()

            if run_alphamelts:
                self.run_alphamelts_all_p(**kwargs)

            # retrieve fo2
            self.read_melts_fo2(**kwargs)

            # retrieve phases
            self.read_melts_phases(which='mass', **kwargs)

            if compare_buffer == 'qfm':
                try:
                    logfo2_buffer = pf.read_qfm_os(self.data['T(K)'].to_numpy(), self.data['P(bar)'].to_numpy(),
                                                   verbose=False, perplex_path=perplex_path)
                    # print('logfo2_qfm', logfo2_buffer)
                    # print('logfo2 = QFM + ', logfo2 - logfo2_buffer, 'at', P, 'bar,', T, 'K')
                    self.data['logfo2_qfm'] = logfo2_buffer
                    self.data['delta_qfm'] = self.data.logfo2 - logfo2_buffer
                except NotImplementedError:
                    # probably didn't get to 1100 C
                    pass
            okay = True
        except FileNotFoundError as e:
            print(e)
            # alphamelts did not work at all for whatever reason
            okay = False

        # store mega df
        df_save = self.data.loc[:, ~self.data.columns.duplicated()].copy()
        # print('df saved\n', df_save.head())
        if save:
            df_save.to_csv(self.output_path + self.name + '_results.csv', sep="\t")
            if verbose:
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

            try:
                df = pd.read_csv(output_file, skiprows=3, index_col=None, sep=r"\s+",
                                 dtype=np.float64).tail(1)
            except FileNotFoundError:
                # this run didn't complete
                print(self.name, 'did not complete, no alphamelts results found')
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

    def read_fo2_results(self, verbose=True):

        if np.isnan(self.data['P(bar)'].to_numpy().all()):
            try:
                self.data = pd.read_csv(self.output_path + self.name + '_results.csv', sep='\t')
                if verbose:
                    print('loaded df\n', self.data.head())
            except FileNotFoundError:
                print('...results.csv file not found! skipping')
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
        idx = find_nearest_idx(pressure, 1)
        self.logfo2_1GPa = self.logfo2[idx]
        try:
            self.delta_qfm_1GPa = self.delta_qfm[idx]
        except (KeyError, AttributeError):
            pass

        # save 4 GPa i.e. garnet stability field
        idx = find_nearest_idx(pressure, 4)
        self.logfo2_4GPa = self.logfo2[idx]
        try:
            self.delta_qfm_4GPa = self.delta_qfm[idx]
        except (KeyError, AttributeError):
            pass

    def read_phase_main(self, phase, p_of_interest, T_of_interest, verbose=False):
        filename = self.output_path + str(float(p_of_interest * 1e4)).replace('.', ',') + 'bar/Phase_main_tbl.txt'
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
        if not start:
            # never found this phase
            if verbose:
                print(phase, 'not found')
            return None
        if len(tmp) < 2:
            print(T_of_interest, 'K not found')
            return None

        # write temp file (not sure how to make df otherwise)
        with open('tmp.txt', 'w') as f:
            for line in tmp:
                f.write(f"{line}\n")

        df = pd.read_csv('tmp.txt', skiprows=0, index_col=None, sep=r"\s+")
        # print(df.head())

        os.remove('tmp.txt')
        return df

    def read_phase_comp(self, p_of_interest, T_of_interest, component='Fe2O3', phases=map_to_px_phase.keys(),
                        absolute_abundance=True, verbose=False):
        """
        Parameters
        ----------
        phases :
        T_of_interest :
        component :
        p_of_interest : GPa
        verbose :

        Returns
        -------

        """
        if absolute_abundance and (not hasattr(self, 'data')):
            self.data = pd.read_csv(self.output_path + self.name + '_results.csv', sep='\t')

        try:
            idx = self.pressures_of_interest.index(p_of_interest*1e4)
        except ValueError as e:
            # pressure not found
            print(p_of_interest, 'GPa not found')
            raise e

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
                            mass_ph = self.data['X_' + map_to_px_phase[phase]].iloc[idx] / 100  # these are wt%
                            print('mass', phase, mass_ph)
                        except KeyError as e:
                            # this is probably because T_of_interest not found - won't fill in _results.csv
                            # for 1M_88Ceff_HIP84856_999K_3,0fer, other weird issue
                            self.read_melts_phases(T_of_interest=T_of_interest, which='mass', verbose=True, reload_TP=True)
                            print(self.data.head())
                            print(self.name, e)
                            return None
                    else:
                        mass_ph = 1

                    try:
                        wt_pt_dict[map_to_px_phase[phase]] = df[component].iloc[0] * mass_ph
                    except KeyError:
                        print(component, 'not found in', phase)
                        wt_pt_dict[map_to_px_phase[phase]] = np.nan
        except FileNotFoundError:
            print('...results.csv file not found! skipping')
            return None
        return wt_pt_dict


def init_from_results(name, output_parent_path=output_parent_default, alphamelts_path=alphamelts_path_default,
                      T_final=1373, load_results_csv=False, verbose=False, X_ferric=None, **kwargs):
    import re
    parts = name.split('_')
    try:
        star = parts[2]
        X_ferric = float(re.findall('\d+\.\d+', parts[4].replace(',', '.'))[0] ) * 1e-2
    except IndexError:
        # non conventional name
        star = None
    wt_oxides = {}

    subfolders = [f.name for f in os.scandir(output_parent_path + name + '/') if f.is_dir()]

    # check if solution
    if (len(subfolders) > 0) and ((os.path.isfile(output_parent_path + name + '/' + subfolders[0] + '/pl.melts')) or ((os.path.isfile(output_parent_path + name + '/' + subfolders[0] + '/' + name + '.melts')))):

        # load pressures
        pressures_of_interest = [float(s.replace(',', '.').replace('bar', '')) for s in subfolders]  # TODO parse directories

        # sort ascending (alphabetical order might be different in os.scandir
        pressures_of_interest.sort()

        # parse melts file
        try:
            melts_file_contents = open(output_parent_path + name + '/' + subfolders[0] + '/pl.melts').readlines()
        except FileNotFoundError as e:
            try:
                melts_file_contents = open(output_parent_path + name + '/' + subfolders[0] + '/' + name + '.melts').readlines()
            except FileNotFoundError:
                print(e, 'skipping...')
                return None

        for line in melts_file_contents:
            if 'Initial Composition:' in line:
                parts = line.split()
                wt_oxides[parts[2]] = float(parts[3])

        dat = MeltsFugacityData(name=name, star=star, X_ferric=X_ferric, wt_oxides=wt_oxides, output_parent_path=output_parent_path,
                                alphamelts_path=alphamelts_path, pressures_of_interest=pressures_of_interest, T_final=T_final, verbose=verbose)

        # load results csv
        if load_results_csv:
            try:
                dat.data = pd.read_csv(dat.output_path + name + '_results.csv', sep='\t')
                dat.read_fo2_results(verbose=verbose)

                if verbose:
                    print('loaded df\n', dat.data.head())
            except FileNotFoundError:
                print('...results.csv file not found! skipping')

        return dat
    else:
        if verbose:
            print('no melts files at', name)
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
        dat = init_from_results(sub, output_parent_path=output_parent_path, verbose=False, **kwargs)
        if dat:  # existing runs
            dat.find_common_T_final_from_results()
            df.at[row, 'name'] = dat.name
            df.at[row, 'T_min'] = dat.T_min
            df.at[row, 'mass_melt_min'] = dat.mass_melt_min
            df.at[row, 'n_pressures'] = dat.n_pressures

    print(df)
    if store:
        df.to_csv(output_parent_path + 'completion_analysis.csv', sep="\t")


def fo2_from_local(output_parent_path, num=-1, **kwargs):
    """for existing melts runs (e.g. done remotely), do the same fo2 analysis"""

    # get all runs in directory
    subfolders = [f.name for f in os.scandir(output_parent_path) if f.is_dir()]
    bad = []
    for name in subfolders[:num]:
        dat = init_from_results(name, output_parent_path=output_parent_path, **kwargs)
        if dat is not None:
            okay = dat.fo2_calc(save=True, run_alphamelts=False, **kwargs)
        if (dat is None) or (not okay):
            bad.append(name)
    return bad








