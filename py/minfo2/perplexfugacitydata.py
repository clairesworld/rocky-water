import numpy as np
import pandas as pd
import os
import pathlib
import subprocess
from py import perplexdata as px
from py import bulk_composition as bulk
from py import main as rw
from scipy import interpolate
from py import parameters as p
from py.useful_and_bespoke import find_nearest_idx

# import ask_hypatia as hyp

Rb = 8.31446261815324
output_parent_default = '/home/claire/Works/min-fo2/perplex_output/'
output_parent_apollo = '/raid1/cmg76/perple_x/output/rocky-fo2/'
wt_oxides_DMM = {'SiO2': 44.71, 'MgO': 38.73, 'CaO': 3.17, 'Al2O3': 3.98, 'FeO': 8.18}  # Workman & Hart depleted mantle
solution_phases_default = ['Cpx(JH)', 'O(JH)', 'Sp(JH)', 'Grt(JH)', 'Opx(JH)', 'Pl(JH)']

wt_oxides_DMM_ext = {'SiO2': 44.71, 'MgO': 38.73, 'CaO': 3.17, 'Al2O3': 3.98, 'FeO': 8.17,
                     'Na2O': 0.28, 'TiO2': 0.13}  # Workman & Hart depleted mantle, extended

map_from_JH_phase = {'O(JH)': 'Ol', 'Opx(JH)': 'Opx', 'Cpx(JH)': 'Cpx', 'Sp(JH)': 'Sp', 'Grt(JH)': 'Gt',
                     'Pl(JH)': 'Plag'}
map_to_JH_phase = {v: k for k, v in map_from_JH_phase.items()}


class PerplexFugacityData(px.PerplexData):

    def __init__(self, name='test', star=None, solution_phases=solution_phases_default,
                 output_parent_path=output_parent_default, wt_oxides=wt_oxides_DMM, X_ferric=0.03,
                 perplex_path=px.perplex_path_default, core_efficiency=None,
                 **kwargs):
        # super().__init__(output_parent_path=output_parent_path, **kwargs)

        self.name = name
        self.star = star
        self.X_ferric = X_ferric
        self.core_eff = core_efficiency

        if 'O2' not in wt_oxides:
            # get O2 mantle fraction
            self.wt_oxides = bulk.insert_o2_from_wt_comp(wt_oxides, X_ferric=X_ferric)
        else:
            self.wt_oxides = wt_oxides
        self.oxide_list = [k for k in self.wt_oxides.keys()]
        self.get_mgsi()
        self.get_mg_number()
        self.solution_phases = solution_phases
        self.output_path = output_parent_path + self.name + '/'
        self.perplex_path = perplex_path

        self.pressure = None  # GPa
        self.temperature = None  # K
        self.logfo2 = None
        self.delta_qfm = None
        self.logfo2_1GPa = None  # spinel field
        self.logfo2_4GPa = None  # garnet field
        self.delta_qfm_1GPa = None
        self.delta_qfm_4GPa = None

        self.data = pd.DataFrame(columns=['P(bar)', 'T(K)', 'logfo2'])

        # remove TiO2 from bulk composition - no phase that incorporates it
        if 'TiO2' in self.oxide_list:
            self.oxide_list.remove('TiO2')
            self.wt_oxides.pop('TiO2')

    def print_comp(self):
        print('wt.% oxides\n-----------')
        for k in self.wt_oxides:
            print("{0:<7}".format(k), "{:5.2f}%".format(self.wt_oxides[k]))
        print('Mg/Si', self.mgsi)

        X_fer_back = bulk.test_ferric_ratio_from_O2(self.wt_oxides['FeO'], self.wt_oxides['O2'])
        print('Fe3+/Fe', X_fer_back)
        print('-----------')

        # insert Fe2O3 like MELTS - FeO already represents all mantle Fe - already gets printed
        wt_dict_o2 = bulk.insert_fe2o3({k: v for k, v in self.wt_oxides.items() if k != 'O2'}, X_fer_back)

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
        """ string for werami command file to get phase proportions """

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

    def command_werami_phase_feiii_composition_grid(self,
                                                    points_file=px.perplex_path_default + 'Earth_1600K_adiabat.dat',
                                                    build_file_end='', **kwargs):
        """ string for werami command file to get chemical potentials """

        s = self.name + build_file_end + '\n'  # Enter the project name (the name assigned in BUILD)
        s = s + '4\n'  # Select operational mode: 4 - as in 3 [properties along a 1d path], but input from file
        s = s + '2\n'  # Path will be described by:  2 - a file containing a list of x-y points
        s = s + points_file + '\n'  # Enter the file name:
        s = s + '1\n'  # every nth plot will be plotted, enter n:
        for ph in self.solution_phases:
            for i_ox in [self.oxide_list.index('FeO') + 1, self.oxide_list.index('O2') + 1]:
                # print('i_ox', i_ox, 'oxide_list', self.oxide_list)
                s = s + '8\n'  # Select a property: Composition (Mol, Mass, or Wt%) of a solution phase
                s = s + ph + '\n'  # Enter solution (left justified):
                s = s + 'n\n'  # Define the composition in terms of the species/endmembers of Opx(JH)    (y/n)?
                s = s + '1\n'  # How many components in the numerator of the composition (<15)?
                s = s + str(i_ox) + '\n1\n'  # Enter component indices and weighting factors for the numerator:
                s = s + str(len(self.wt_oxides)) + '\n'  # How many components in the denominator of the composition
                for i in range(1, len(self.wt_oxides) + 1):
                    s = s + str(i) + '\n1\n'  # Enter component indices and weighting factors for the denominator:
                s = s + 'n\n'  # Change it (y/n)? This composition will be designated: C[Opx(JH)1]
        s = s + '0\n'  # Select an additional property or enter 0 to finish:
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

    def ferric_composition_calc(self, T_iso=None, points_file=None, p_min=1000, p_max=40000, verbose=True,
                                build_file_end='', **kwargs):
        # decide which T, p points to calculate
        if T_iso is not None:
            if points_file is not None:
                raise Exception('ERROR: cannot take both points_file and isotherm, set one of these to None')
            try:
                # run werami over an isothermal section
                points_file = self.write_isothermal_points(T0=T_iso, p_min=p_min, p_max=p_max, delta_p=1000)
            except TypeError as e:
                print('ERROR: trying to write an isotherm section but isotherm is not a valid numeric')
                raise e
        elif points_file is None:
            raise NotImplementedError('ERROR: need either input points_file or input isotherm')
        self.run_perplex(werami_command_end='_werami_command_ferric_comp.txt',
                         werami_command_text_fn=self.command_werami_phase_feiii_composition_grid,
                         vertex_command_text_fn=self.command_vertex_grid,
                         output_file_end='_ferric_comp.tab', build_file_end=build_file_end, verbose=verbose,
                         clean=True, werami_kwargs={'points_file': points_file}, run_vertex=False,
                         # suppress_output=False,
                         **kwargs)

    def read_ferric_phase_comp(self, phase=None, T_of_interest=None, p_of_interest=None, scale=100):
        import re
        # load data
        df_w = self.read_werami(fend='_ferric_comp.tab')

        # find idx of T of interest
        try:
            irow = df_w.loc[(df_w['T(K)'] == T_of_interest) & (df_w['P(bar)'] == p_of_interest * 1e4)].index[0]
        except IndexError:
            print(self.name, T_of_interest, p_of_interest * 1e4, 'not found')
            print(df_w['P(bar)'].to_numpy())
            return None
        # convert FeO + O2 into Fe2O3
        for icol in range(2, len(df_w.columns), 2):  # skip first 2 (T, P), every 1st is FeO
            col = df_w.columns[icol]
            ph = map_from_JH_phase[re.sub(r'[0-9]', '', col[2:-1])]  # remove C[*] and trailing digit
            if ph == phase:
                ph_FeOstar = df_w.iloc[irow, icol]  # raw results has two columns per phase (FeO* and O2)
                ph_O2 = df_w.iloc[irow, icol + 1]
                # print('\n\n\n\nFeO* in', ph, ph_FeOstar)
                # print('O2 in', ph, ph_O2)

                # convert to Fe2O3
                if (ph_O2 == 0) or (ph_FeOstar == 0):
                    m_Fe2O3 = 0
                    m_FeO = 0  # not used
                else:
                    Xf = bulk.test_ferric_ratio_from_O2(ph_FeOstar, ph_O2)
                    m_FeO, m_Fe2O3 = bulk.partition_FeOstar(ph_FeOstar, Xf)

                # print('Xf', Xf, 'm_Fe2O3', m_Fe2O3)

                # print('returns n_o2', bulk.o2_molar_ratio(m_FeO / p.M_FeO, X_ferric=Xf))

                return pd.DataFrame({'P(bar)': [df_w.iloc[irow, 0]], 'T(K)':
                    [df_w.iloc[irow, 1]], 'Fe2O3': [m_Fe2O3 * scale]})

            # # process in results df
            # if len(self.data) == 0:
            #     try:
            #         self.data = pd.read_csv(self.output_path + self.name + '_results.csv', sep='\t', index_col=0)
            #     except FileNotFoundError:
            #         raise Exception(self.name, 'results.csv not found')

    def read_phase_comp(self, p_of_interest, T_of_interest, component='Fe2O3', phases=solution_phases_default,
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
            self.data = pd.read_csv(self.output_path + self.name + '_results' + str(int(T_of_interest)) + '.csv',
                                    sep='\t')

        try:
            wt_pt_dict = {map_from_JH_phase[ph]: None for ph in phases}
        except KeyError:
            # maybe already converted from JH phase to stx
            wt_pt_dict = {ph: None for ph in phases}
        try:
            for ph in phases:
                try:
                    phase = map_from_JH_phase[ph]
                except KeyError:
                    phase = ph
                df = self.read_ferric_phase_comp(phase=phase, T_of_interest=T_of_interest, p_of_interest=p_of_interest)
                if df is None:  # phase not found
                    wt_pt_dict[phase] = np.nan
                else:
                    if verbose:
                        print(phase, 'loaded df\n', df.head())

                    if absolute_abundance:
                        # normalise to total quantity of phase
                        try:
                            idx = self.data.loc[(self.data['P(bar)'] == p_of_interest * 1e4)].index[0]
                        except ValueError as e:
                            # pressure not found
                            print(p_of_interest, 'GPa not found')
                            raise e

                        try:
                            mass_ph = self.data['X_' + phase].iloc[idx]
                            # print('mass', phase, mass_ph)
                        except KeyError as e:
                            # print(self.data.head())
                            if verbose:
                                print('Phase', e, 'not found in', self.name, '--> assigning 0 wt%')
                            mass_ph = 0
                    else:
                        mass_ph = 1

                    try:
                        wt_pt_dict[phase] = df[component].iloc[0] * mass_ph
                    except KeyError:
                        print(component, 'not found in', phase)
                        wt_pt_dict[phase] = 0
        except FileNotFoundError:
            print('...results.csv file not found! skipping')
            return None
        return wt_pt_dict

    def get_phase_composition_dict(self, p_of_interest, T_of_interest, component='Fe2O3',
                                   phases=solution_phases_default,
                                   to_absolute_abundance=True, verbose=False):
        """
        load phase compositon from results.csv file
        todo maybe put these in results csv? for now can just copy over raw output files lol

        Parameters
        ----------
        p_of_interest :
        T_of_interest :
        component :
        phases :
        absolute_abundance :
        verbose :

        Returns
        -------

        """
        return self.read_phase_comp(p_of_interest, T_of_interest, component=component, phases=phases,
                                    absolute_abundance=to_absolute_abundance, verbose=verbose)

    def fo2_calc(self, p_min=1000, p_max=40000, T_min=1600, T_max=1603, T_iso=None, verbose=True, build_file_end='',
                 vertex_data='hp622ver', run=True, compare_buffer=None, do_system_comp=False, points_file=None,
                 mu0_file=None, save=True, excluded_phases=[], run_vertex='auto', **kwargs):
        """
        Parameters
        ----------
        p_min :
        p_max :
        T_min :
        T_max :
        T_iso : Numeric or None
            Temperature value for doing calculations for an isothermal section, None if not applicable
        verbose :
        build_file_end :
        vertex_data :
        run :
        compare_buffer :
        do_system_comp :
        points_file :
        mu0_file : Path or string
            Given relative to perplex_path
        save : bool
            If True, auto-save important results to csv

        Returns
        -------

        """
        # if verbose:
        #     print(self.__dict__)

        # decide which T, p points to calculate
        if T_iso is not None:
            if points_file is not None:
                raise Exception('ERROR: cannot take both points_file and isotherm, set one of these to None')
            try:
                # run werami over an isothermal section
                points_file = self.write_isothermal_points(T0=T_iso, p_min=p_min, p_max=p_max, delta_p=1000)
            except TypeError as e:
                print('ERROR: trying to write an isotherm section but isotherm is not a valid numeric')
                raise e
        elif points_file is None:
            raise NotImplementedError('ERROR: need either input points_file or input isotherm')

        if run:
            # create build file
            self.write_build(title='Planet', p_min=p_min, p_max=p_max,  # adiabat_file=points_file,
                             verbose=verbose, overwrite=True, vertex_data=vertex_data, option_file='perplex_option_fo2',
                             excluded_phases=excluded_phases,  # melt phases
                             use_solutions=True, build_file_end=build_file_end,
                             calculation_type='5',  # '5' for grid
                             T_min=T_min, T_max=T_max, **kwargs)
            # todo: re-writing build file unnecessarily in perplex wd if it already exists in output folder

            # run vertex and werami to get mu at T, P of interest
            # strategy is to run each case once in vertex over grid and keep all the files
            print('\nCalculating mu...')
            self.run_perplex(werami_command_end='_werami_command_mu.txt',
                             werami_command_text_fn=self.command_werami_mu,
                             vertex_command_text_fn=self.command_vertex_grid,
                             output_file_end='_mu.tab', build_file_end=build_file_end, verbose=verbose,
                             clean=True,  ## store_vertex_output means files don't get deleited
                             werami_kwargs={'points_file': points_file}, run_vertex=run_vertex,
                             store_vertex_output=True,
                             **kwargs)

            if do_system_comp:
                print('\nCalculating composition')
                self.run_perplex(werami_command_end='_werami_command_comp.txt',
                                 werami_command_text_fn=self.command_werami_composition_grid,
                                 vertex_command_text_fn=self.command_vertex_grid,
                                 output_file_end='_comp.tab', build_file_end=build_file_end, verbose=verbose,
                                 clean=True, werami_kwargs={'points_file': points_file}, run_vertex=False,
                                 store_vertex_output=True,
                                 **kwargs)
                self.ferric_composition_calc(points_file=points_file, T_iso=None, p_min=p_min, p_max=p_max,
                                             verbose=True)
                print('done Fe3+ composition for', self.name)

        # extract mu and T of interest from output data
        df_w = self.read_chem_potential_werami(verbose=False, return_df=True)
        mu = df_w['mu_O2(J/mol)'].to_numpy()
        T = df_w['T(K)'].to_numpy()
        P = df_w['P(bar)'].to_numpy()

        # prepare all data
        self.data = df_w.copy()

        if do_system_comp:
            df_c = self.read_werami(fend='_comp.tab')

            # concat with full results dataframe
            px_phases = [col for col in df_c.columns if col not in ['P(bar)', 'T(K)']]
            df_c.rename(columns={phase: 'X_' + phase for phase in px_phases}, inplace=True)
            self.data = pd.concat([self.data, df_c], axis=1)

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
        self.data['logfo2'] = logfo2

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

            logfo2_buffer = read_qfm_os(T, P, verbose=False, perplex_path=self.perplex_path)
            # print('logfo2_qfm', logfo2_buffer)
            # print('logfo2 = QFM + ', logfo2 - logfo2_buffer, 'at', P, 'bar,', T, 'K')
            self.data['logfo2_qfm'] = logfo2_buffer
            self.data['delta_qfm'] = logfo2 - logfo2_buffer

        # store mega df
        if save:
            df_save = self.data.loc[:, ~self.data.columns.duplicated()].copy()
            df_save.to_csv(self.output_path + self.name + '_results' + str(int(T_iso)) + '.csv', sep="\t")

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

    def read_fo2_results(self, verbose=False):
        if np.isnan(self.data['P(bar)'].to_numpy().all()):
            raise Exception('ERROR: PerplexFugacityData not initialised properly - not loaded dataframe')

        self.logfo2 = self.data['logfo2'].to_numpy()
        try:
            self.delta_qfm = self.data['delta_qfm'].to_numpy()
        except KeyError:
            pass

        pressure = self.data['P(bar)'].to_numpy() * 1e-4
        temperature = self.data['T(K)'].to_numpy()
        if self.pressure is None:
            self.pressure = pressure
        elif (self.pressure != pressure).any():
            print('pressure', pressure)
            print('self.pressure', self.pressure)
            raise Exception('attempting to read in mismatching pressure data!')
        if self.temperature is None:
            self.temperature = temperature
        elif (self.temperature != temperature).any():
            print('temperature', temperature)
            print('self.temperature', self.temperature)
            raise Exception('attempting to read in mismatching temperature data!')

        # save 1 GPa i.e. spinel stability field
        try:
            idx = find_nearest_idx(pressure, 1)
        except ValueError as e:
            print(self.name)
            print('pressure', pressure)
            print('self.data\n', self.data.head())
            raise e
        self.logfo2_1GPa = self.logfo2[idx]
        try:
            self.delta_qfm_1GPa = self.delta_qfm[idx]
        except KeyError:
            pass

        # save 4 GPa i.e. garnet stability field
        idx = find_nearest_idx(pressure, 4)
        self.logfo2_4GPa = self.logfo2[idx]
        try:
            self.delta_qfm_4GPa = self.delta_qfm[idx]
        except KeyError:
            pass


def init_from_results(name, X_ferric=None, T_iso=1373, load_results_csv=False, verbose=False, **kwargs):
    """ currently not saving X_ferric in directory name so must enter known value manually
    output_parent_path in kwargs"""
    # TODO parse star from name using possible prefixes

    # from py.parameters import M_E

    spl = name.split('_')
    # M_p = float(''.join(filter(str.isdigit, spl[0]))) * M_E  # get planet mass - not good if decimal
    try:
        core_efficiency = int(''.join(filter(str.isdigit, spl[1]))) / 100.0
        star = spl[2]
    except (IndexError, ValueError):
        star = None
        core_efficiency = None
    kwargs.update({'name': name, 'X_ferric': X_ferric, 'star': star, 'core_efficiency': core_efficiency})

    try:
        kwargs.update(**read_dict_from_build(**kwargs))  # get other important input stuff from build file
        dat = PerplexFugacityData(**kwargs)
    except FileNotFoundError:
        # if verbose:
        #     print('perplex input data not found:', name)
        return None

    # load saved results
    if load_results_csv:
        dat.data = pd.read_csv(dat.output_path + dat.name + '_results' + str(int(T_iso)) + '.csv', sep='\t',
                               index_col=0)
        # print('self.data\n', dat.data.head())
        dat.read_fo2_results(verbose=verbose)

    return dat


def read_dict_from_build(name=None, output_parent_path=output_parent_default, verbose=False, **kwargs):
    import re
    """ parse the parameters file into a python dictionary """

    fin = output_parent_path + name + '/' + name + '.dat'
    if verbose:
        print("Reading build parameters from", fin)

    build_file = open(fin).readlines()

    # read thermodynamic data (always 1st line?)
    vertex_data = build_file[0].split()[0].replace(".dat", "")

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
         'solution_phases': solution_phases, 'p_min': float(p_min), 'p_max': float(p_max), 'T_min': float(T_min),
         'T_max': float(T_max)}

    # if verbose:
    #     print('d', d)
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


def read_qfm_os(T, P, perplex_path=px.perplex_path_default, fin='data_tables/fmqNl_fo2_oli.dat', verbose=False,
                **kwargs):
    """ read in oli's qfm calculations, T in K, p in bar

    Parameters
    ----------
    T :
    verbose :
    fin :
    perplex_path :
    P :
    **kwargs :
    """

    def do(df, T0, p0):
        # first check if exact values already there
        exists = None
        if p0 in df['P(bar)'].unique():
            exists = 'P'
            df = df[df['P(bar)'] == p0]

        if T0 in df['T(K)'].unique():
            if exists:  # both T and P already there
                # print('p0', p0, 'T0', T0)
                df = df.loc[(df['P(bar)'] == p0)]
                df = df.loc[(df['T(K)'] == T0)]
                return np.squeeze(df['logfo2'])
                # return np.squeeze(df.loc[(df['P(bar)'] == p0) and (df['T(K)'] == T0)]['logfo2'])
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

    df = pd.read_csv(perplex_path + fin, delimiter=r"\s+", index_col=False, header=None,
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
            print('P', P)
            print('T', T)
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


def fo2_from_hypatia(p_min, p_max, n_sample=5, core_efficiency=0.88, T_iso=None, planet_kwargs={},
                     output_parent_path=output_parent_default, skip_existing=True, **kwargs):
    planet_kwargs.update({'core_efficiency': core_efficiency, 'solve_interior': False})
    pl_list = rw.planets_from_hypatia(n_sample=n_sample, plot_all=False,
                                      get_saturation=False,
                                      stopafter=None, output_parent_path=output_parent_path,
                                      **planet_kwargs, **kwargs)
    print('\nfinished generating compositions\n')
    bad = []
    for pl in pl_list:
        if skip_existing and os.path.exists(pl.output_path + pl.name + '_results' + str(int(T_iso)) + '.csv'):
            print('skipping', pl.name, ': results.csv file exists')
        else:
            # print('not found', pl.output_path + pl.name + '_results.csv')
            # print('running\n', pl.wt_oxides)
            okay = fo2_from_oxides(pl.name, p_min, p_max, pl=pl, output_parent_path=output_parent_path, T_iso=T_iso,
                                   **kwargs)
            if not okay:
                bad.append(pl.name)
    print('bad cases:', bad)
    return pl_list


def fo2_from_hypatia_1D(p_min, p_max, n_sample=5, core_efficiency=0.88, T_iso=None, planet_kwargs={},
                        output_parent_path=output_parent_default, skip_existing=True, **kwargs):
    """
    Go through cases which haven't finished, start in 1D - TODO

    Parameters
    ----------
    p_min :
    p_max :
    n_sample :
    core_efficiency :
    planet_kwargs :
    output_parent_path :
    skip_existing :
    kwargs :

    Returns
    -------

    """
    planet_kwargs.update({'core_efficiency': core_efficiency, 'solve_interior': False})
    pl_list = rw.planets_from_hypatia(n_sample=n_sample, plot_all=False,
                                      get_saturation=False,
                                      stopafter=None, output_parent_path=output_parent_path,
                                      **planet_kwargs, **kwargs)
    print('\nfinished generating compositions\n')
    bad = []
    for pl in pl_list:
        if skip_existing and os.path.exists(pl.output_path + pl.name + '_results' + str(int(T_iso)) + '.csv'):
            print('skipping', pl.name, ': results.csv file exists')
        else:
            # print('not found', pl.output_path + pl.name + '_results.csv')
            print('running\n', pl.wt_oxides)
            okay = fo2_from_oxides(pl.name, p_min, p_max, pl=pl, T_iso=T_iso, output_parent_path=output_parent_path,
                                   **kwargs)
            if not okay:
                bad.append(pl.name)
    print('bad cases:', bad)
    return pl_list


def fo2_from_local(output_parent_path=output_parent_default, T_iso=1373, run_werami=True,
                   do_ferric_comp=True, do_mu_comp=False, do_system_comp=False,
                   p_min=10000, p_max=40000, skip_names=[], start_after=None,
                   rewrite_options=True, cases=None, **kwargs):
    # perform fo2 calculations on local (existing) vertex data in entire directory
    if run_werami:
        run = True  # note this also re-writes build file (e.g. if options updated)
    else:
        run = False  # don't need to rerun any perplex stuff, just re-calc fo2.
    points_file = 'points_files/' + str(int(T_iso)) + 'K_isotherm.xy'

    subfolders = rw.get_run_dirs(output_path=output_parent_path)
    if not subfolders:
        print('no local output found in', output_parent_path)
        return None

    if start_after is None:
        idx0 = 0
    else:
        idx0 = [os.path.basename(sub) for sub in subfolders].index(start_after) + 1

    for sub in subfolders[idx0:]:
        name = os.path.basename(sub)
        # star = name.split('_')[2]

        if name in skip_names:
            continue
        elif (cases is not None) and name not in cases:
            continue

        try:
            d = read_dict_from_build(name=name, output_parent_path=output_parent_path, verbose=False)
        except FileNotFoundError as e:
            print('FileNotFound:', e)
            continue

        dat = PerplexFugacityData(name=name, output_parent_path=output_parent_path, **d, **kwargs)
        # print('pressures', dat.pressure)
        # print('d', d)
        # logfo2 = dat.fo2_calc(T_iso=None, run=run, points_file=points_file, run_vertex=False, **d, **kwargs)
        # print(dat.name, ': log fo2 of system:', logfo2)

        if do_ferric_comp:
            dat.ferric_composition_calc(points_file=points_file, T_iso=None, p_min=p_min, p_max=p_max,
                                        **kwargs)
            print('done Fe3+ composition for', dat.name)


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

    try:
        tmp = pl.wt_oxides.items()
    except AttributeError as e:
        print(e)
        print(name, 'has no wt_oxides composition')
        return False
    print('\n\n\n\n--------------------------------------------------------------------------------')
    print('Starting fO2 calculation for planet of', pl.star)
    dat = PerplexFugacityData(name=name, wt_oxides=pl.wt_oxides, verbose=verbose, output_parent_path=output_parent_path,
                              core_efficiency=core_efficiency, **kwargs)

    logfo2 = dat.fo2_calc(p_min=p_min, p_max=p_max, T_min=T_min, T_max=T_max, verbose=verbose, **kwargs)
    if verbose:
        print('log fo2 of system:', logfo2)
    return True


# def get_name(M_p=None, star=None, core_efficiency=None, Tp=None, test_CMF=None, test_oxides=None, suffix=None,
#              **kwargs):
#     mass_str = str(M_p / p.M_E) + 'M'
#     if test_CMF is not None:
#         cmf_str = str(int(test_CMF * 100)) + 'CMF'
#     else:
#         cmf_str = str(int(core_efficiency * 100)) + 'Ceff'
#     if test_oxides is not None:
#         comp_str = str(int(np.round(test_oxides['MgO']))) + 'Mg'
#     elif star is not None:
#         comp_str = star.replace(' ', '')
#     else:
#         raise NotImplementedError('no bulk composition scenario input (need star or test_oxides)')
#     if 'profile' in kwargs:
#         if kwargs['profile'] != 'adiabat':
#             temp_str = kwargs['profile']
#     elif Tp is not None:
#         temp_str = str(int(Tp)) + 'K'
#     else:
#         raise NotImplementedError('no temperature scenario input (need Tp or profile name)')
#     name = mass_str + '_' + cmf_str + '_' + comp_str + '_' + temp_str
#     if suffix is not None:
#         name = name + '_' + suffix
#     name = name.replace('.', ',')  # dots will crash file i/o
#     return name


def create_isothermal_csv(output_parent_path, T, P, Xfer, coreeff, verbose=True, **kwargs):
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
    dfout['Fe3/Fe(%)'] = Xfer
    dfout['Fe_c/Fe_T(%)'] = coreeff

    for row, name in enumerate(subfolders):
        # ensure name is there
        parts = name.split('_')
        star = parts[2]
        dfout.at[row, 'star'] = star

        dat = init_from_results(name, X_ferric=Xfer, load_results_csv=True, verbose=False, **kwargs)
        if dat is not None:
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
                        pass  # this phase not stable at T,p

                for s in ['Mg/Si', 'Fe/Si', 'Al/Si', 'Ca/Si']:
                    dfout.at[row, s] = bulk.get_element_ratio(s, dat.wt_oxides)

    dfout.to_csv(output_parent_path + 'summary_perplex_' + str(T).replace('.', ',') + 'K_' + str(P).replace('.',
                                                                                                            ',') + 'GPa_' + str(
        Xfer).replace('.', ',') + 'fer_' + str(coreeff).replace('.', ',') + 'core' + '.csv', sep="\t", na_rep='NaN')
    if verbose:
        print(dfout.head(10))
    pd.reset_option("max_columns")
    print('Done!')
    return None
