import numpy as np
# import matplotlib.pyplot as plt
import pandas as pd
from py.parameters import M_E, M_Fe, M_FeO, M_MgO, M_SiO2, M_Si, M_Mg, M_Ca, M_CaO, M_Al, M_Al2O3, G, R_E, rho_E
import os
import pathlib
import subprocess
from py.bulk_composition import stellar_mantle
import py.ask_hypatia as hyp

perplex_path_default = '/home/claire/Works/perple_x/'  # path to perple_x installation (https://www.perplex.ethz.ch/)
output_parent_default = perplex_path_default + 'output/'

# set defaults and benchmarks
wt_oxides_Mars = {'SiO2': 42.71, 'MgO': 31.5, 'CaO': 2.49, 'Al2O3': 2.31,
                  'FeO': 18.71}  # nominal Perple_x SNC meteorite composition, but outdated - FeO is ~14% (Khan+2022)
wt_oxides_Earth = {'SiO2': 44.48, 'MgO': 39.22, 'CaO': 3.44, 'Al2O3': 3.59, 'FeO': 8.10}  # weight percent, 'Na2O': 0.36
wt_oxides_MD95 = {'SiO2': 45.0, 'MgO': 37.8, 'CaO': 3.55, 'Al2O3': 4.45, 'FeO': 8.05}  # McDonough & Sun 1995 BSE
wt_oxides_pyrolite = {'SiO2': 38.71, 'MgO': 49.85, 'FeO': 6.17, 'CaO': 2.94, 'Al2O3': 2.22}  # 'Na2O': 0.11

oxide_list_default = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO']  # 'NaO
solution_phases_default = ['O', 'Sp', 'Cpx', 'Wad', 'Ring', 'Pv', 'Wus', 'C2/c', 'Opx', 'Aki', 'Ppv', 'Gt', 'CF',
                           'Pl']  # 'NAl'


class PerplexData:
    def __init__(self, name='test', core_efficiency=0.88, M_p=M_E, R_p=None,
                 oxides=None, solution_phases=None,
                 star='sun', perplex_path=perplex_path_default, output_parent_path=output_parent_default, verbose=False,
                 **kwargs):
        self.c_h2o_obm = None
        if solution_phases is None:
            solution_phases = solution_phases_default
        if oxides is None:
            oxides = oxide_list_default
        self.name = name
        self.oxide_list = oxides
        self.solution_phases = solution_phases
        self.star = star
        self.M_p = M_p
        self.core_eff = core_efficiency
        self.perplex_path = perplex_path
        self.output_path = output_parent_path + self.name + '/'

        # ensure output_parent_default is a valid directory, create directory if it doesn't exist yet
        # if not os.path.isdir(self.output_path):
        #     raise Exception(self.output_path, 'is not a valid directory with input parameter output_parent_path',
        #                     output_parent_path)
        if not os.path.exists(self.output_path):
            # create directory if it doesn't exist yet
            print('creating directory', self.output_path)
            os.makedirs(self.output_path)

        # initiate other variables
        self.nH_star = None
        self.wt_oxides = None
        self.mgsi = None
        self.mg_number = None
        self.R_p = R_p
        self.R_c = None
        self.CMF = None
        self.radius = None
        self.density = None
        self.gravity = None
        self.temperature = None
        self.pressure = None
        self.alpha = None
        self.cp = None
        self.cum_mass = None
        self.mass = None
        self.i_cmb = None
        self.p_mtz = None
        self.p_lm = None
        # self.pressure_m = None
        # self.temperature_m = None
        self.phases_px = None
        self.df_comp = None
        self.Fe_mass_fraction = None
        self.phase_mass = None
        self.mass_um = None
        self.mass_obm = None
        self.mass_h2o_total = None
        self.mass_h2o_um = None
        self.mass_h2o_obm = None
        self.c_h2o_mantle = None

        # if verbose:
        #     print('----------------------\ninitialising PerplexData object with M_p = {:.4e}%'.format(M_p), 'kg,',
        #           'core_eff = {:5.2f}%'.format(core_efficiency))

    def load_composition_px(self, build_file_end='', fillnan=True, check_consistent=True, verbose=True,
                            save=True, interp_missing=True, **kwargs):
        file = self.perplex_path + self.name + build_file_end + '_comp.tab'

        try:
            df = pd.read_csv(file, skiprows=8, index_col=None, sep=r"\s+",
                             # dtype=np.float64
                             )
        except FileNotFoundError:
            raise Exception('trying to load perple_x composition before running werami')

        if fillnan:
            df = df.fillna(0)

        P = df['P(bar)'].to_numpy() * 1e5  # in Pa
        T = df['T(K)'].to_numpy()

        if interp_missing:
            # occasionally get a layer with no compositon found in perplex - interpolate from bounding layers
            # weird because still get a density...
            sum_comp = df.iloc[:, 3:].sum(axis=1)
            bad_rows = sum_comp.eq(0)
            for idx in df.loc[bad_rows].index:
                try:
                    hi = df.iloc[idx - 1, 3:]
                    lo = df.iloc[idx + 1, 3:]
                except IndexError:
                    raise Exception('probably trying to interpolate 0 compositon at top or bottom layer')
                df.iloc[idx, 3:] = (hi + lo) / 2
                print('...interpolating missing layer at node', idx + 1)

        # if check_consistent and hasattr(self, 'pressure_m') and (self.pressure_m != P).any():
        #     raise Exception('current pressure inconsistent with loaded composition data - probably because you were '
        #                     'using thermodynamic output from an older iteration')
        # elif save:
        #     self.pressure_m = P
        # if check_consistent and hasattr(self, 'temperature_m') and (self.temperature_m != T).any():
        #     raise Exception('current temperature inconsistent with loaded composition data')
        # elif save:
        #     self.temperature_m = T

        # print('original phases list:', df.columns.values.tolist())

        # find repeating (e.g. Ppv in some cases)
        dup_cols = [col for col in df.columns.values.tolist() if '.' in col]  # dot put there from pandas reading duplicate behaviour
        for col in dup_cols:
            print('dealing with duplicate', col)
            dup_sum = df[col] + df[col[:-2]]
            df.drop([col, col[:-2]], axis=1, inplace=True)
            df[col[:-2]] = dup_sum
        phases = df.columns.values.tolist()

        # remove df columns that aren't a mineral phase
        phases.remove('node#')
        phases.remove('P(bar)')
        phases.remove('T(K)')
        if save:
            self.phases_px = phases
            self.df_comp = df
            if verbose:
                print('-----------------\nphases used in Perple_x:', [ph for ph in self.phases_px])
                print('-----------------')

        return df, phases

    def load_adiabat(self, build_file_end='', fillnan=True, check_consistent=True, fix_nan=True, run_bad=False,
                     store=False, head=False, rho_m_av=None, **kwargs):
        file = self.perplex_path + self.name + build_file_end + '_thermo.tab'

        df = pd.read_csv(file, skiprows=8, index_col=None, sep=r"\s+",
                         # dtype=np.float64
                         )
        if fillnan:
            df = df.fillna(0)
        if head:
            print(df.head())
        try:
            P = df['P(bar)'].to_numpy() * 1e5  # in Pa
            T = df['T(K)'].to_numpy()
            alpha_m = df['alpha,1/K'].to_numpy()
            cp_m = df['cp,J/K/kg'].to_numpy()
            rho_m = df['rho,kg/m3'].to_numpy()
        except AttributeError:
            # old version of pandas
            P = np.array(df['P(bar)']) * 1e5  # in Pa
            T = np.array(df['T(K)'])
            alpha_m = np.array(df['alpha,1/K'])
            cp_m = np.array(df['cp,J/K/kg'])
            rho_m = np.array(df['rho,kg/m3'])
        # check if any nans or zeros - flag and for now just extrapolate constant value
        # will want to see how many cases later (limits of Stixrude 2021 dataset e.g.)
        # need to do this for rho, alpha, cp as they could be different idk
        ierr_rho = np.argmax(rho_m == 0)
        ierr_alpha = np.argmax(alpha_m == 0)
        ierr_cp = np.argmax(cp_m == 0)
        if (ierr_cp > 0) or (ierr_rho > 0) or (ierr_alpha > 0):
            arri = np.array([ierr_cp, ierr_rho, ierr_alpha])
            maxp = P[np.min(arri[np.nonzero(arri)])]
            print('warning: NaN in thermodynamic output from p =', maxp * 1e-9, 'GPa')
            if run_bad:
                # run again with all output so you can see error - note this causes problems if overwriting
                vertex_command_file = self.name + build_file_end + '_vertex_command.txt'
                with open(self.perplex_path + vertex_command_file, 'w') as file:
                    s = self.name + build_file_end + '\n0'
                    file.write(s)
                os.system('./vertex < ' + vertex_command_file)
                werami_command_file = self.name + build_file_end + '_werami_command_thermo.txt'
                with open(self.perplex_path + werami_command_file, 'w') as file:
                    s = self.command_werami_thermo(build_file_end=build_file_end)
                    file.write(s)
                os.system('./werami < ' + werami_command_file)
                os.remove(werami_command_file)
                os.remove(vertex_command_file)
                for fend in ['_seismic_data.txt', '_auto_refine.txt', '_1.plt', '.tim', '.plt', '.blk',
                             '.arf', '.tof']:
                    if os.path.isfile(self.perplex_path + self.name + build_file_end + fend):
                        os.remove(self.perplex_path + self.name + build_file_end + fend)
            # raise Exception('Perple_x EoS error: cp = 0 at idx', i)
            if fix_nan:
                if ierr_cp > 0:
                    idx = np.where(cp_m == 0)[0]
                    for ii in idx:
                        cp_m[ii] = cp_m[ii + 1]
                    # cp_m[ierr_cp:] = cp_m[ierr_cp - 1]
                if ierr_alpha > 0:
                    idx = np.where(alpha_m == 0)[0]
                    for ii in idx:
                        alpha_m[ii] = alpha_m[ii + 1]
                    # alpha_m[ierr_alpha:] = alpha_m[ierr_alpha - 1]
                if ierr_rho > 0:
                    idx = np.where(rho_m == 0)[0]
                    for ii in idx:
                        rho_m[ii] = rho_m[ii + 1]
                    # rho_m[ierr_rho:] = rho_m[ierr_rho - 1]

        # if store:
        #     # update PerplexData object
        #     if check_consistent and hasattr(self, 'pressure_m') and (self.pressure_m != P).any():
        #         raise Exception('current pressure inconsistent with loaded thermodynamic data')
        #     else:
        #         self.pressure_m = P
        #     if check_consistent and hasattr(self, 'temperature_m') and (self.temperature_m != T).any():
        #         raise Exception('current temperature inconsistent with loaded thermodynamic data')
        #     else:
        #         self.temperature_m = T

            # self.alpha_m = alpha_m  # note upside down compared to full planet profile
            # self.cp_m = cp_m
            # self.rho_m = rho_m
        return rho_m, alpha_m, cp_m, P, T

    def core_mass_fraction(self):
        """
        Calculates core mass fraction given bulk FeO in mantle and the wt fraction of total FeO that ends up in core.
        Doesn't require knowing the planet mass
        core_eff = n_Fe_core / (n_FeO_mantle + n_Fe_core) = m_Fe_core / (m_Fe_mantle + m_Fe_core) as n_Fe = n_FeO in mtl
        """
        # todo: still might want to triple check this is consistent with stellar composition constraints
        # try:
        #     idx = self.oxide_list.index('FEO')
        # except ValueError:
        #     idx = self.oxide_list.index('FeO')
        # try:
        #     x_Fe_mm = self.wt_oxides[idx] * (M_Fe / M_FeO)  # mass fraction Fe in mantle wrt mtl
        # except KeyError:
        #     x_Fe_mm = self.wt_oxides['FeO'] * (M_Fe / M_FeO)  # mass fraction Fe in mantle wrt mtl

        try:
            x_Fe_mm = self.wt_oxides['FeO'] * (M_Fe / M_FeO)  # mass fraction Fe in mantle wrt total mantle mass
        except KeyError:
            raise Exception(
                'ERROR: cmf is undefined for a given molar core efficiency if no FeO in bulk mantle. try fixing input cmf instead.')

        if self.core_eff == 1:
            # TODO
            x_Fe_cm = 1 - x_Fe_mm  # mass fraction Fe in core wrt total mantle mass
        else:
            x_Fe_cm = -x_Fe_mm * self.core_eff / (self.core_eff - 1)  # mass fraction Fe in core wrt total mantle mass
        mantle_mass_fraction = 100 / (100 + x_Fe_cm)

        # scale by mass_mtl/M_p --> should be smaller than x_Fe_cm
        x_Fe_c = x_Fe_cm * mantle_mass_fraction
        # print('x_Fe_c wrt planet', x_Fe_c)
        self.CMF = x_Fe_c / 100

        self.Fe_mass_fraction = x_Fe_c + x_Fe_mm * mantle_mass_fraction  # store bulk planet Fe fraction for posterity
        # print('x_Fe_pl', x_Fe_pl)

        print('\ncore mass fraction = {:5.3f}\n'.format(self.CMF))
        return self.CMF

    def core_eff_from_cmf(self):
        M_c = self.CMF * self.M_p  # core mass in kg
        M_m = self.M_p - M_c  # mantle mass in kg
        n_Fe_mtl = self.wt_oxides['FeO'] / 100 * M_m / M_FeO  # n_Fe = n_FeO
        n_Fe_core = M_c / M_Fe  # pure Fe core
        self.core_eff = n_Fe_core / (n_Fe_core + n_Fe_mtl)
        return self.core_eff

    def get_star_compositon(self, **kwargs):
        if self.star == 'sun':
            from parameters import ca_sol, fe_sol, al_sol, mg_sol, si_sol, na_sol
            solar = {'ca_sol': ca_sol, 'fe_sol': fe_sol, 'al_sol': al_sol, 'mg_sol': mg_sol, 'si_sol': si_sol, 'na_sol': na_sol}
            self.nH_star = [solar[ox[:2].lower() + '_sol'] for ox in self.oxide_list if ox != 'O2']
        else:
            if 'oxide_list' in kwargs:
                pass
                # print('oxide list in kwargs', kwargs['oxide_list'])
            else:
                kwargs['oxide_list'] = [ox for ox in self.oxide_list if ox != 'O2']
            self.nH_star = hyp.star_composition(**kwargs)

    def write_star_composition(self, fname='nH_star.txt', path=None):
        if path is None:
            path = self.output_path
        try:
            np.savetxt(path + fname, self.nH_star, delimiter=',')  # X is an array
        except (IndexError, ValueError) as e:  # different errors for different python versions
            print(e)
            print('self.nH_star =', self.nH_star)
            print('ERROR: star composition possibly incomplete for', self.star, '; no file written')

    def get_mgsi(self, **kwargs):
        # if self.nH_star:
        #     self.mgsi = 10 ** self.nH_star[self.oxide_list.index('MgO')] / 10 ** self.nH_star[self.oxide_list.index('SiO2')]  # molar ratio
        # else:
        n_MgO = self.wt_oxides['MgO'] / M_MgO  # convert to moles
        n_SiO2 = self.wt_oxides['SiO2'] / M_SiO2
        self.mgsi = n_MgO / n_SiO2  # n_Mg = n_MgO etc
        return self.mgsi

    def get_mg_number(self):
        # Mg number of mantle
        n_MgO = self.wt_oxides['MgO'] / M_MgO  # convert to moles
        n_FeO = self.wt_oxides['FeO'] / M_FeO
        self.mg_number = n_MgO / (n_MgO + n_FeO) * 100
        return self.mg_number

    def get_femg_star(self):
        X_ratio_mol = 10 ** self.nH_star[self.oxide_list.index('FeO')] / 10 ** self.nH_star[self.oxide_list.index('MgO')]
        self.femg_star = X_ratio_mol

    def get_fesi_star(self):
        X_ratio_mol = 10 ** self.nH_star[self.oxide_list.index('FeO')] / 10 ** self.nH_star[self.oxide_list.index('SiO2')]
        self.fesi_star = X_ratio_mol

    def get_phase_masses(self):
        self.phase_mass = {ph: np.sum(self.df_comp[ph]*1e-2 * self.df_all['mass(kg)']) for ph in self.phases_px}
        return self.phase_mass

    def get_um_mass(self):
        from saturation import total_water_mass
        i_um_base = self.find_lower_mantle() - 1  # base of upper mantle
        # print('pressure at um base', self.df_all['P(bar)'][i_um_base], 'bar')
        self.mass_um = np.sum(self.df_all['mass(kg)'][:i_um_base + 1])

        self.mass_h2o_um = total_water_mass(self.df_all, i_min=0, i_max=i_um_base)
        return self.mass_um

    def star_to_oxide(self, **kwargs):
        """ get the bulk oxide compositon of the mantle (core Fe will have been extracted from returned dict) """
        wt_oxides = stellar_mantle(self.oxide_list, self.nH_star, self.core_eff)
        self.wt_oxides = wt_oxides
        return wt_oxides

    def write_build(self, build_file_end='', title='Planet', p_min=10000, p_max=245000, adiabat_file='aerotherm.dat',
                    verbose=False, overwrite=True, vertex_data='stx21ver', option_file='perplex_option_claire',
                    excluded_phases=None, use_solutions=True, calculation_type='10', T_min=0.0, T_max=0.0, **kwargs):
        """ write perple_x build file, p in bar but have fucked it before and np """
        if excluded_phases is None:
            excluded_phases = []  # no excluded phases / solution end members
        # if verbose:
        #     print('using thermodynamic dataset:', vertex_data)
        if 'stx21' in vertex_data:
            solution_model = 'stx21_solution_model'
        elif 'stx11' in vertex_data:
            solution_model = 'stx11_solution_model'
        else:
            solution_model = 'solution_model'

        build_file = self.perplex_path + self.name + build_file_end + '.dat'
        if os.path.isfile(build_file) and not overwrite:
            raise Exception('WARNING: build file', build_file, 'already exists, set overwrite=True')
        elif os.path.isfile(build_file):
            print('  overwriting', build_file)

        s = ''
        with open(build_file, 'w') as file:
            s = s + vertex_data + '.dat     thermodynamic data file\n'
            s = s + 'no_print | print generates print output\n'
            s = s + 'plot     | obsolete 6.8.4+\n'
            s = s + solution_model + '.dat     | solution model file, blank = none\n'
            s = s + title + '\n'
            s = s + option_file + '.dat | Perple_X option file\n'
            s = s + '   ' + calculation_type + ' calculation type: 0- composition, 1- Schreinemakers, 3- Mixed, 4- swash, 5- gridded min, 7- 1d fract, 8- gwash, 9- 2d fract, 10- 7 w/file input, 11- 9 w/file input, 12- 0d infiltration\n'
            if calculation_type == '10':
                s = s + adiabat_file + '     | coordinate file \n'
            s = s + '    0 unused place holder, post 06\n' * 9
            s = s + '    0 number component transformations\n'
            s = s + '    ' + str(len(self.wt_oxides)) + ' number of components in the data base\n'
            s = s + '    1 component amounts, 0 - mole, 1 mass\n'
            s = s + '    0 unused place holder, post 06\n' * 2
            s = s + '    0 unused place holder, post 05\n'
            s = s + '    5 ifug EoS for saturated phase\n'
            s = s + '    2 gridded minimization dimension (1 or 2)\n'
            s = s + '    0 special dependencies: 0 - P and T independent, 1 - P(T), 2 - T(P)\n'
            s = s + ' 0.00000      0.00000      0.00000      0.00000      0.00000     Geothermal gradient polynomial coeffs.\n\n'

            s = s + 'begin thermodynamic component list\n'
            for el, wt in self.wt_oxides.items():
                dig = len(str(int(wt)))
                if vertex_data == 'hp622ver':
                    el_database = el
                else:
                    el_database = el.upper()
                s = s + el_database.ljust(6) + '1'.ljust(3) + "{value:{width}.{precision}f}".format(value=float(wt),
                                                                                                   width=7,
                                                                                                   precision=6 - dig)
                s = s + '      0.00000      0.00000     mass  amount\n'
            s = s + 'end thermodynamic component list\n\n\n'

            s = s + 'begin saturated component list\nend saturated component list\n\n\n'
            s = s + 'begin saturated phase component list\nend saturated phase component list\n\n\n'
            s = s + 'begin independent potential/fugacity/activity list\nend independent potential list\n\n\n'

            s = s + 'begin excluded phase list\n'
            for ph in excluded_phases:
                s = s + ph + '\n'
            s = s + 'end excluded phase list\n\n\n'

            s = s + 'begin solution phase list'
            if use_solutions:
                for sol in self.solution_phases:
                    s = s + '\n' + sol
            s = s + '\nend solution phase list\n\n'

            try:
                s = s + "   {value:{width}.{precision}f}".format(value=float(p_max), width=7,
                                                                 precision=7 - len(str(int(p_max))))[
                        :-1] + "        {value:{width}.{precision}f}".format(value=float(T_max), width=7,
                                                                 precision=7 - len(str(int(T_max))))[
                        :-1] + '         0.00000        0.00000        0.00000     max p, t, xco2, mu_1, mu_2\n'
            except ValueError as e:  # not enough digits
                s = s + "   {value:{width}.{precision}f}".format(value=float(p_max), width=10,
                                                                 precision=10 - len(str(int(p_max))))[
                        :-1] + "        {value:{width}.{precision}f}".format(value=float(T_max), width=10,
                                                                 precision=10 - len(str(int(T_max))))[
                        :-1] + '        0.00000        0.00000        0.00000     max p, t, xco2, mu_1, mu_2\n'
            s = s + "   {value:{width}.{precision}f}".format(value=float(p_min), width=7,
                                                             precision=7 - len(str(int(p_min))))[
                    :-1] + "        {value:{width}.{precision}f}".format(value=float(T_min), width=7,
                                                                 precision=7 - len(str(int(T_min))))[
                        :-1] + '        0.00000        0.00000        0.00000     min p, t, xco2, mu_1, mu_2\n'
            s = s + '   0.00000        0.00000        0.00000        0.00000        0.00000     unused place holder post 06\n\n'
            s = s + ' 1  2  4  5  3   indices of 1st & 2nd independent & sectioning variables'

            file.write(s)
        if verbose:
            print('   wrote to', build_file)

    def write_adiabat(self, P, T, file_end='_adiabat', fout=None, verbose=False, overwrite=True, **kwargs):
        """ write perple_x build file, p in bar """
        if fout is None:
            fout = self.perplex_path + self.name + file_end + '.dat'
        if os.path.isfile(fout) and not overwrite:
            raise Exception('WARNING: file', fout, 'already exists, set overwrite=True')
        elif os.path.isfile(fout):
            print('  overwriting', fout)

        s = ''
        with open(fout, 'w') as file:
            for p, t in zip(P, T):
                s = s + "{:.6e}".format(p) + '	' + "{:.6e}".format(t) + '\n'
            file.write(s)
        if verbose:
            print('  wrote to', fout)

    def command_vertex(self, build_file_end='', **kwargs):
        """ string for vertex command file - entries are project name and 0 in operational mode 10 """
        s = self.name + build_file_end + '\n0'
        return s

    def command_werami_composition(self, build_file_end='', **kwargs):
        """string for werami command file to get compositional data"""
        s = self.name + build_file_end + '\n'  # Enter the project name (the name assigned in BUILD)
        s = s + '3\n'  # Select operational mode: 3 - properties along a 1d path
        s = s + '25\n'  # Select a property: 25 - Modes of all phases
        s = s + 'n\n'  # Output cumulative modes (y/n)?
        s = s + '0\n'  # Select operational mode: 0 - EXIT
        return s

    def command_werami_thermo(self, build_file_end='', **kwargs):
        """string for werami command file to get thermodynamic data"""
        s = self.name + build_file_end + '\n'  # Enter the project name (the name assigned in BUILD)
        s = s + '3\n'  # Select operational mode: 3 - properties along a 1d path
        s = s + '2\n'  # Select a property: 2 - Density (kg/m3)
        s = s + 'n\n'  # Calculate individual phase properties (y/n)?
        s = s + '19\n'  # Select an additional property or enter 0 to finish: 19 - Heat Capacity (J/K/kg)
        s = s + 'n\n'  # Calculate individual phase properties (y/n)?
        s = s + '4\n'  # Select an additional property or enter 0 to finish: 4 - Expansivity (1/K, for volume)
        s = s + 'n\n'  # Calculate individual phase properties (y/n)?
        s = s + '0\n'  # Select an additional property or enter 0 to finish:
        s = s + '0\n'  # Select operational mode: 0 - EXIT
        return s

    def run_perplex(self, werami_command_end='_werami_command.txt', build_file_end='', suppress_output=True, clean=True,
                    verbose=True, vertex_command_text_fn=None, werami_command_text_fn=None, run_vertex=True,
                    store_vertex_output=False, output_file_end='.tab', werami_kwargs=None, **kwargs):
        """ run_vertex = 'auto' will check for vertex output files in run output directory, False will require output
        files already in perple_x working directory"""
        if werami_kwargs is None:
            werami_kwargs = {}
        cwd = os.getcwd()
        os.chdir(self.perplex_path)
        if suppress_output:
            stderr, stdout = subprocess.DEVNULL, subprocess.DEVNULL
        else:
            stderr, stdout = None, None

        vertex_copy_flag = False  # track if you moved vertex files
        if run_vertex == 'auto':
            run_vertex = False
            for fend in ['_seismic_data.txt', '_auto_refine.txt', '.tim', '.plt', '.blk', '.arf', '.tof']:
                if not os.path.isfile(self.output_path + self.name + build_file_end + fend):
                    run_vertex = True  # if any of these files are missing, need to run
                    break
                else:
                    # temporarily move vertex output files to perple_x working directoy
                    vertex_copy_flag = True
                    # print('copying to', self.perplex_path + self.name + build_file_end + fend)
                    os.rename(self.output_path + self.name + build_file_end + fend,
                              self.perplex_path + self.name + build_file_end + fend)
            print('   vertex output files already exist for', self.name, ', skipping to werami')

        if run_vertex:  # this takes longer so if you kept files (clean=False) might not need to run again
            # create vertex command file
            vertex_command_file = self.name + build_file_end + '_vertex_command.txt'
            with open(self.perplex_path + vertex_command_file, 'w') as file:
                s = vertex_command_text_fn(build_file_end=build_file_end)
                file.write(s)

            # run vertex
            output = subprocess.run('./vertex < ' + vertex_command_file, shell=True,
                                    stdout=stdout, stderr=stderr)
            if verbose:
                print('  ', output)

        # create werami command file
        werami_command_file = self.name + build_file_end + werami_command_end
        with open(self.perplex_path + werami_command_file, 'w') as file:
            s = werami_command_text_fn(build_file_end=build_file_end, **werami_kwargs)
            file.write(s)

        print('   created werami command file', self.name + build_file_end + werami_command_end)

        # run werami
        output = subprocess.run('./werami < ' + werami_command_file, shell=True,
                                stdout=stdout, stderr=stderr)
        if verbose:
            print('  ', output)

        # delete extra files and rename .tab file to something meaningful
        try:
            os.rename(self.perplex_path + self.name + build_file_end + '_1.tab',
                      self.perplex_path + self.name + build_file_end + output_file_end)
        except FileNotFoundError as e:
            print('ERROR: vertex did not complete, try running again with suppress_output=False')
            # something probably went wrong with vertex or werami, run again with full output
            # os.system('./vertex < ' + vertex_command_file)
            # os.system('./werami < ' + werami_command_file)
            raise e

        if store_vertex_output or vertex_copy_flag:
            # move vertex files to this run's output directory
            for fend in ['_seismic_data.txt', '_auto_refine.txt', '.tim', '.plt', '.blk', '.arf', '.tof',
                         '_vertex_command.txt']:
                if os.path.isfile(self.perplex_path + self.name + build_file_end + fend):
                    os.rename(self.perplex_path + self.name + build_file_end + fend,
                              self.output_path + self.name + build_file_end + fend)
        if clean:
            # get rid of any files left in
            for fend in ['_seismic_data.txt', '_auto_refine.txt', '_1.plt', '.tim', '.plt', '.blk', '.arf', '.tof',
                         '_vertex_command.txt', werami_command_end]:
                if os.path.isfile(self.perplex_path + self.name + build_file_end + fend):
                    os.remove(self.perplex_path + self.name + build_file_end + fend)

        # return to original dir
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        # os.chdir(cwd)

    def werami_adiabat(self, **kwargs):
        self.run_perplex(werami_command_end='_werami_command_thermo.txt',
                         werami_command_text_fn=self.command_werami_thermo,
                         vertex_command_text_fn=self.command_vertex,
                         output_file_end='_thermo.tab', **kwargs)

    def werami_composition(self, **kwargs):
        self.run_perplex(werami_command_end='_werami_command_comp.txt',
                         werami_command_text_fn=self.command_werami_composition,
                         vertex_command_text_fn=self.command_vertex,
                         output_file_end='_comp.tab', **kwargs)

    def get_last_composition(self, build_file_end='', **kwargs):
        self.werami_composition(build_file_end=build_file_end, **kwargs)
        df, _ = self.load_composition_px(build_file_end=build_file_end, save=False, **kwargs)
        df = df.iloc[-1, 3:]  # get last row, drop P and T columns
        df = df.loc[df != 0]  # drop columns with 0
        phases = [s.replace('-', '') for s in df.index]  # this is just to work with burnman later
        proportions = df.to_numpy()
        # print('phases', phases)
        # print('proportions', proportions)
        return phases, proportions

    def extend_lm_composition(self, build_file_end='', save=True, **kwargs):
        """ get deepest composition ran in perple_x and hold constant to CMB"""
        self.werami_composition(build_file_end=build_file_end, **kwargs)
        df_um, phases = self.load_composition_px(build_file_end=build_file_end, **kwargs)

        n = len(self.pressure)
        n_um = len(df_um)
        n_lm = n - n_um - len(self.pressure[:self.i_cmb + 1])

        # print('UM')
        # print(df_um.tail())
        # print(df_um.index)

        idx_lm = range(len(df_um), len(df_um) + n_lm)

        df_lm = pd.DataFrame(index=idx_lm)
        for phase in phases:
            # print('UM phase', phase)
            df_lm[phase] = df_um[phase].iloc[-1]  # get last row
            # print('value', df_um[phase].iloc[-1])

        st, en = self.i_cmb + 1, n - n_um
        pressure_lm = self.pressure[st: en]
        temperature_lm = self.temperature[st: en]

        # print('st, en', st, en)

        df_lm['P(bar)'] = pressure_lm[::-1] * 1e-5
        df_lm['T(K)'] = temperature_lm[::-1]

        # print('LM')
        # print(df_lm.head())
        # print(df_lm.index)

        df_comp = pd.concat((df_um, df_lm), axis=0, join='outer')

        if save:
            self.df_comp = df_comp
            # self.pressure_m = df_comp['P(bar)'] * 1e5
            # self.temperature_m = df_comp['T(K)']
        return df_comp

    def werami_garnet_composition(self, build_file_end='', **kwargs):
        def command_werami_gt_composition(**kwargs):
            # string for werami command file to get compositional data
            s = self.name + build_file_end + '\n'  # Enter the project name (the name assigned in BUILD)
            s = s + '3\n'  # Select operational mode: 3 - properties along a 1d path
            for endmember in (1, 2, 3, 4):  # (maj, gr, alm, py)
                s = s + '8\n'  # Select a property: 8 - Composition (Mol, Mass, or Wt%) of a solution phase
                s = s + 'Gt\n'  # Enter solution (left justified):
                s = s + 'y\n'  # Define the composition in terms of the species/endmembers of Gt         (y/n)?
                s = s + '1\n'  # How many species in the numerator of the composition (<15)?
                s = s + str(endmember) + ' 1\n'  # Enter species indices and weighting factors for the numerator:
                s = s + '0\n'  # How many species in the denominator of the composition? Enter zero to use the numerator
                s = s + 'n\n'  # Change it (y/n)?
            s = s + '0\n'  # Select an additional property or enter 0 to finish:
            s = s + '0\n'  # Select operational mode: 0 - EXIT
            return s

        # need to copy build file and adiabat back to perplex directory if missing
        flag = False
        for fend in ('.dat', '_adiabat.dat'):
            if not os.path.isfile(self.perplex_path + self.name + build_file_end + fend):
                flag = True
                os.rename(self.output_path + self.name + build_file_end + fend,
                          self.perplex_path + self.name + build_file_end + fend)
        self.run_perplex(werami_command_end='_werami_command_gtcomp.txt',
                         werami_command_text_fn=command_werami_gt_composition,
                         vertex_command_text_fn=self.command_vertex,
                         output_file_end='_gtcomp.tab', **kwargs)
        # move back
        if flag:
            for fend in ('.dat', '_adiabat.dat'):
                os.rename(self.perplex_path + self.name + build_file_end + fend,
                          self.output_path + self.name + build_file_end + fend)
        os.rename(self.perplex_path + self.name + build_file_end + '_gtcomp.tab',
                  self.output_path + self.name + build_file_end + '_gtcomp.tab')
        # todo: replace Gt[1] with maj, etc...

    def iterate_structure(self, Psurf=1000, Tp=1600, n='auto', maxIter=100, tol=0.0001, clean=True,
                          deltaT_cmb=0, rho_m0=None, profile='adiabat', parameterise_lm=True,
                          p_max_perplex=200e9 * 1e-5,  # maximum pressure in bar for perplex if parameterise_lm=True
                          core_density_fudge=None,
                          **kwargs):
        """ tweaked from Noack+
        Tsurf is potential surface temperature in K (was 1600), Psurf is surface pressure in bar
        n is radial resolution from center of planet to surface, must be <~ 2000 else perple_x error
        tol is fractional convergence criterion for iteration solver
        rho_m0 is initial guess for mantle density (optional)
        parameterise_lm extrapolates a constant composition deeper than 200 GPa"""
        import eos
        import math

        M = self.M_p  # planet mass in kg
        Mc = M * self.CMF
        try:
            x_Fe = self.CMF  # this is just for starting guess on Rp from parameterisation
        except AttributeError:
            raise Exception('must run star_to_oxide() and core_mass_fraction() to get planet Fe fraction')

        if n == 'auto':
            if M / M_E <= 1:
                n = 1200 #300  # n = 200 saves about 7 sec per run, misses a little Aki phase
            elif M / M_E <= 2:
                n = 1200 #500
            elif M / M_E <= 2.5:
                n = 1200
            elif M / M_E <= 5:
                n = 1600

        if profile == 'warm':
            Tp = self.melting_temperature(**kwargs)  # use compositionally-dependent solidus temperature
            print('using Tp', Tp, 'K')
        # elif profile == 'adiabat': Tp = Tp

        # Initialization - guesses
        # rho_c_av = 11000  # guess for core density in kg/m^3
        # rho_m_av = 4000  # guess for mantle density in kg/m^3
        cp_c = 800  # guess for core heat capacity in J/kg K
        cp_m = 1300  # guess for mantle heat capacity in J/kg K
        alpha_c = 0.00001  # guess for core thermal expansion ceoff. in 1/K
        alpha_m = 0.000025  # guess for mantle thermal expansion ceoff. in 1/K
        Rp = 1e3 * (7030 - 1840 * x_Fe) * (
                    M / M_E) ** 0.282  # initial guesss, Noack & Lasbleis 2020 (5) ignoring mantle Fe
        if self.CMF > 0:
            Rc = 1e3 * 4850 * x_Fe ** 0.328 * (M / M_E) ** 0.266  # initial guess, hot case, ibid. (9)
            rho_c_av = x_Fe * M / (4 / 3 * np.pi * Rc ** 3)
        else:
            Rc = 0
            rho_c_av = 0
        if rho_m0 is None:
            rho_m_av = (1 - x_Fe) * M / (4 / 3 * np.pi * (Rp ** 3 - Rc ** 3))  # Noack & Lasbleis parameterisation
        else:
            rho_m_av = rho_m0  # use initial guess as given
            Rc, Rp = eos.core_planet_radii(x_Fe, M, rho_c_av, rho_m_av)  # get consistent radius

        # Arrays
        radius = np.zeros(n)  # corresponds to height at top of layer
        density = np.zeros(n)
        alpha = np.zeros(n)
        cp = np.zeros(n)
        mass = np.zeros(n)  # cumulative mass, not differential

        # Initialization of arrays: surface values
        radius[-1] = Rp  # radius of top of layer, from centre
        density[-1] = rho_m_av
        alpha[-1] = alpha_m
        cp[-1] = cp_m
        mass[-1] = M
        if x_Fe == 1:  # pure iron shell
            density[-1] = rho_c_av
            alpha[-1] = alpha_c
            cp[-1] = cp_c
            Rc = Rp

        # Initialization of arrays: interpolation over depth from top
        for i in range(2, n + 1):  # goes to i = n-1
            mass[n - i] = mass[n - i + 1] - M / n  # such that it never goes to 0 (would cause numerical errors)
        if self.CMF > 0:
            i_cmb = np.argmax(mass > Mc)  # index of cmb in profiles
        else:
            i_cmb = 0

        density[i_cmb + 1:] = rho_m_av
        alpha[i_cmb + 1:] = alpha_m
        cp[i_cmb + 1:] = cp_m
        density[:i_cmb + 1] = rho_c_av
        alpha[:i_cmb + 1] = alpha_c
        cp[:i_cmb + 1] = cp_c

        # get radius from volumes
        for i in range(2, n + 1):  # goes to i = n-1
            dmass = mass[n - i + 1] - mass[n - i]
            radius[n - i] = np.cbrt(-(dmass / density[n - i] / (4 * np.pi / 3)) + radius[n - i + 1] ** 3)

        gravity = eos.g_profile(n, radius, density)
        pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, Psurf, Tp, i_cmb, deltaT_cmb)
        p_cmb = pressure[i_cmb]  # the pressure (at cell top edge) of the cell that Rc passes through
        p_cmb_guess = np.mean(gravity[i_cmb + 1:]) * rho_m_av * (
                    Rp - Rc)  # not used but neat that it can be not far off
        p_mantle_bar = pressure[i_cmb + 1:][::-1] * 1e-5  # convert Pa to bar and invert
        T_mantle = temperature[i_cmb + 1:][::-1]
        print('initial guesses | p_cen =', pressure[0] * 1e-9, 'GPa | p_cmb =', p_cmb * 1e-9, 'GPa | Rp =', Rp / R_E,
              'R_E')

        # Iteration
        print('>>>>>>>>>\nIterating interior structure...')
        it = 1
        iter_param_old = 1e-5
        iter_param = p_cmb  # Rp
        while (abs((iter_param - iter_param_old) / iter_param_old) > tol) and (it < maxIter):
            #     print(it, 'criterion:', abs((Rp - Rp_old)/Rp_old), '>', tol)
            # store old value to determine convergence
            iter_param_old = iter_param

            run_flag = True
            for i in range(n):  # M: 1:n, P: 0:n-1  index from centre to surface
                if i <= i_cmb:  # get local thermodynamic properties - core - layer by layer
                    if pressure[i] > 10e3 * 1e9:
                        print('pressure error, i', i, pressure[i] * 1e-9, 'GPa', 'rho cen', density[0], 'rho cmb',
                              density[i_cmb])
                        # plt.plot(radius * 1e-3, pressure * 1e-9)
                        # plt.axvline(radius[i_cmb] * 1e-3)
                        # plt.show()
                    _, density[i], alpha[i], cp[i] = eos.EOS_all(pressure[i] * 1e-9, temperature[i], 4)

                    if core_density_fudge:
                        density[i] = density[i] * core_density_fudge

                    if cp[i] == 0:
                        print('i', i, 'cp[i]', cp[i], 'problem with core EoS')
                        raise ZeroDivisionError

                else:  # run perple_x to get mantle EOS
                    if run_flag:

                        if parameterise_lm and pressure[
                            i_cmb + 1] * 1e-5 > p_max_perplex:  # only run perple_x to 200 GPa and parameterise below
                            i_lm = np.argmax(
                                p_mantle_bar > p_max_perplex)  # index of 200 GPa in profiles (this layer is inclusive)
                            p_mantle_compute = p_mantle_bar[:i_lm + 1]
                            T_mantle_compute = T_mantle[:i_lm + 1]
                            # print('running perple_x over pressures', p_mantle_compute, 'bar')
                        else:
                            p_mantle_compute = p_mantle_bar
                            T_mantle_compute = T_mantle

                        # run perple_x over mantle p, T to calculate thermodynamic properties (inverse adiabat idx)
                        self.write_adiabat(p_mantle_compute, T_mantle_compute, file_end='_temp' + str(it) + '_adiabat',
                                           **kwargs)
                        self.write_build(build_file_end='_temp' + str(it),
                                         p_min=p_mantle_compute[0], p_max=p_mantle_compute[-1],
                                         adiabat_file=self.name + '_temp' + str(it) + '_adiabat.dat',
                                         use_solutions=False, **kwargs)  # no solutions
                        self.werami_adiabat(build_file_end='_temp' + str(it), clean=clean, **kwargs)

                        # then after running vertex & werami, extract density, alpha, cp
                        density_wer, alpha_wer, cp_wer, p_wer, T_wer = self.load_adiabat(
                            build_file_end='_temp' + str(it),
                            head=False, check_consistent=False, store=False, **kwargs)

                        if parameterise_lm and pressure[i_cmb + 1] * 1e-5 > p_max_perplex:
                            # store perplex stuff for UM: invert again because perple_x idx in opposite directions
                            n_um = len(density_wer)
                            # print('index cmb', i_cmb)
                            # print('index n - n_um', n - n_um)
                            # print('length of um, density[n - n_um:]', len(density[n - n_um:]))
                            alpha[n - n_um:] = alpha_wer[::-1]  # TODO: check index
                            cp[n - n_um:] = cp_wer[::-1]
                            density[n - n_um:] = density_wer[::-1]

                            # if parameterising, now need to fill in extra bits. start from composition of deepest layer
                            # run werami to get composition
                            phases, proportions = self.get_last_composition(build_file_end='_temp' + str(it),
                                                                            clean=clean, **kwargs)
                            # get properties
                            st, en = i_cmb + 1, n - n_um
                            # print('length of lm', len(range(st, en)))
                            _, density[st: en], alpha[st: en], cp[st: en] = \
                                eos.EOS_lm(pressure[st: en] * 1e-9, temperature[st: en], phases, proportions)

                        else:
                            # store perplex, invert again because perple_x idx in opposite directions
                            alpha[i_cmb + 1:] = alpha_wer[::-1]
                            cp[i_cmb + 1:] = cp_wer[::-1]
                            density[i_cmb + 1:] = density_wer[::-1]

                        run_flag = False  # only run perple_x once per profile

            # update mass and radius - 1st mass entry never changes
            for i in range(2, n + 1):  # goes to i = n-1
                mass[n - i] = mass[n - i + 1] - M / n  # such that it never goes to 0 (would cause numerical errors)
            radius[0] = 0  #np.cbrt(mass[0] / density[0] / (4 * np.pi / 3))
            print('radius[0]', radius)
            for i in range(1, n):
                dmass = mass[i] - mass[i - 1]
                radius[i] = np.cbrt((dmass / density[i] / (4 * np.pi / 3)) + radius[i - 1] ** 3)

            gravity = eos.g_profile(n, radius, density)

            # pressure and temperature are interpolated from surface downwards
            pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, Psurf, Tp, i_cmb,
                                                   deltaT_cmb)

            # i_cmb = np.argmax(radius > Rc)   # update index of top of core in profiles
            i_cmb = np.argmax(mass > Mc)
            Rp = radius[-1]
            Rc = radius[i_cmb]

            if profile == 'solidus':
                # test with melting temperature instead
                temperature[i_cmb + 1:] = self.melting_temperature(p_GPa=pressure[i_cmb + 1:] * 1e-9)

            p_mantle_bar = pressure[i_cmb + 1:][::-1] * 1e-5  # convert Pa to bar
            T_mantle = temperature[i_cmb + 1:][::-1]

            # update parameter that should be converging
            iter_param = pressure[i_cmb]  # Rp
            print(it, "R_p = {:.6e}".format(Rp / R_E), 'R_E', "R_c = {:.4e}".format(Rc / R_E), 'R_E',
                  "p_cmb = {:.4e}".format(pressure[i_cmb] * 1e-9), 'GPa', "p_cen = {:.4e}".format(pressure[0] * 1e-9),
                  'GPa',
                  # "rho_c_av = {:.2e}".format(rho_c_av), 'kg/m3', "rho_m_av = {:.2e}".format(rho_m_av), 'kg/m3',
                  "CMF*M_p = {:.3e}".format(Mc / M_E), 'M_E',
                  'mass[i_cmb] = {:.5e}'.format(mass[i_cmb] / M_E), 'M_E', 'i_cmb', i_cmb)

            # if parameterise_lm:
            #     fig, axes = plt.subplots(1, 5)
            #     axes[0].plot(p_mantle_bar)
            #     axes[1].plot(density[i_cmb + 1:][::-1])
            #     axes[2].plot(alpha[i_cmb + 1:][::-1])
            #     axes[3].plot(cp[i_cmb + 1:][::-1])
            #     axes[4].plot(radius[i_cmb + 1:][::-1])
            #     axes[0].set_ylabel('p bar')
            #     axes[1].set_ylabel('density')
            #     axes[2].set_ylabel('alpha')
            #     axes[3].set_ylabel('cp')
            #     axes[4].set_ylabel('radius')
            #     for ax in axes:
            #         ax.axvline(i_cmb)
            #         ax.axvline(n - n_um)

            if clean:
                for fend in ['.dat', '_adiabat.dat', '_thermo.tab', '_comp.tab', '_WERAMI_options.txt',
                             '_VERTEX_options.txt']:
                    if os.path.isfile(self.perplex_path + self.name + '_temp' + str(it) + fend):
                        os.remove(self.perplex_path + self.name + '_temp' + str(it) + fend)
            it = it + 1
            if it == maxIter:
                print('WARNING: reached maximum iterations for interior solver')
            # end while

        # ultimately only need adiabat for mantle - write this to file
        self.write_adiabat(p_mantle_compute, T_mantle_compute, file_end='_adiabat', **kwargs)

        # update profiles in object attributes - these are centre-of-planet to surface
        self.radius = radius
        self.density = density
        self.gravity = gravity
        print('gravity stored', self.gravity, 'len', len(self.gravity))
        print('radius stored', self.radius, 'len', len(self.radius))
        self.temperature = temperature
        self.pressure = pressure
        self.alpha = alpha
        self.cp = cp
        self.cum_mass = mass
        self.mass = np.diff(np.insert(mass, 0, 0))
        self.R_p = Rp
        self.R_c = Rc
        self.i_cmb = i_cmb
        # print('final length of mantle', np.shape(T_mantle), np.shape(p_mantle_bar))

    def find_lower_mantle(self):
        """ return the index denoting the top of the lower mantle where 0 is the surface
         although not all planets have MTZ minerals or olivine polymorphs, can generalise to where perovskite kicks in,
         even if most of the upper mantle is just SiO2. index will change slightly depending on resolution, but deciding
         to take the first layer because e.g. not all planets have over 50% perovskite in lower mantle. Pv reaches max
         modality within a handful of layers.
         EDIT: when X_pv > X_ring - because if transition is gradual will be biased to lower p - but makes viol plot look weird because can include lm phases
         but at 0 FeO content can get no ring sometimes, so here take first pv?
         """
        try:
            # series = self.df_comp['Pv'].divide(self.df_comp['Ring'])
            # i_lm = series[series.gt(1)].index[0]
            # print('i_lm', i_lm)
            i_lm = self.df_comp['Pv'].ne(0).idxmax()
        except KeyError:
            # if 'X_ring' not in self.df_comp.columns:  # happens at 0 FeO, get aki instead
            #     i_lm = self.df_comp['Pv'].ne(0).idxmax()
            # else:
            print('      no Pv phase found! (planet too small?)')
            return len(
                self.df_comp)  # so retrieving UM will retrieve whole mantle, and retrieving LM will retrieve empty
        except IndexError as e:  # this happens if there is pv but never reaches larger fraction than ring
            print('      not enough of a Pv phase found! (planet too small?)')
            # print(self.name, self.output_path)
            return len(
                self.df_comp)  # so retrieving UM will retrieve whole mantle, and retrieving LM will retrieve empty
        self.p_lm = self.df_comp['P(bar)'].iloc[i_lm] * 1e5  # Pa
        self.z_lm = self.df_all['z(m)'].iloc[i_lm]  # m
        return i_lm

    def find_transition_zone(self):
        """ return the index denoting the top of the transition zone where 0 is the surface
         for planets with no Ol or Wad, still get Ring so use this as MTZ definition...
         index will change slightly depending on resolution """
        try:
            i_mtz = self.df_comp['Wad'].ne(0).idxmax()
        except KeyError:
            print('      no Wad phase found! Defining TZ at base of opx')
            return len(self.df_comp)
        self.p_mtz = self.df_comp['P(bar)'].iloc[i_mtz] * 1e5   # Pa
        return i_mtz

    def get_obm_water(self):
        from saturation import total_water_mass
        i_mtz_base = self.find_transition_zone() - 1  # base of upper mantle
        self.mass_h2o_obm = total_water_mass(self.df_all, i_min=0, i_max=i_mtz_base)
        self.mass_obm = np.sum(self.df_all['mass(kg)'][:i_mtz_base])  # also get mass of layer
        self.c_h2o_obm = self.mass_h2o_obm / self.mass_obm  # also get concentration (avg)

    def get_lm_water(self):
        from saturation import total_water_mass
        i_lm = self.find_lower_mantle()
        self.mass_h2o_lm = total_water_mass(self.df_all, i_min=i_lm)

    def setup_interior(self, test_CMF=None, test_oxides=None, oxides=None, x_Si_core=None, test_nH_star=None,
                     parameterise_lm=True, p_max_perplex=200e9 * 1e-5, solve_interior=True, **kwargs):
        """ procedure for setting up interior composition and structure of planet """

        # first, get oxide bulk composition and CMF (skip if input a testing value)
        if oxides is None:
            oxides = oxide_list_default
        if test_oxides is None:
            # print('kwargs get_interior', kwargs)
            if test_nH_star is None:
                self.get_star_compositon(**kwargs)
                if self.nH_star is None:
                    return None  # e.g. missing element in hypatia catalogue
            else:  # use fixed for testing
                self.nH_star = test_nH_star
            self.star_to_oxide(**kwargs)
        else:
            if self.star == 'sun':
                self.get_star_compositon()
            elif self.star not in ('sun', None):
                print('Did you mean to input both star name and fixed oxide wt%s?')

            # ensure bulk oxide composition sums to 100% (necessary for cmf calculations
            sum_wt_oxides = sum(test_oxides.values())
            self.wt_oxides = {k: test_oxides[k] / sum_wt_oxides * 100 for k in oxides}

        self.get_mgsi()  # molar mg/si ratio
        self.get_mg_number()
        if self.nH_star is not None:
            self.get_femg_star()

        # second, get core size for either a fixed CMF or a core partitioning
        if test_CMF is None:
            self.core_mass_fraction()
        else:
            if test_CMF == 0:
                test_CMF = 1e-10  # can't be 0 for analytical reasons???
            self.CMF = test_CMF
            self.core_eff_from_cmf()

        # pull some Si out of mantle into core - not self-consistently doing core EoS with Si though
        if x_Si_core is not None:
            self.partition_core_Si(x_Si_core)

        if solve_interior:
            # iterate perplex to get interior structure and geotherm
            self.iterate_structure(parameterise_lm=parameterise_lm, p_max_perplex=p_max_perplex, **kwargs)

            # finally, run perple_x again to get mantle compositional profile
            # N. B. this is different than final iteration because including solution phases
            self.write_build(adiabat_file=self.name + '_adiabat.dat', p_min=self.pressure[-1] * 1e-5,
                             p_max=self.pressure[self.i_cmb + 1] * 1e-5, **kwargs)
            if parameterise_lm and self.pressure[self.i_cmb + 1] * 1e-5 > p_max_perplex:
                self.extend_lm_composition(build_file_end='', **kwargs)  # includes get_composition()
            else:
                self.werami_composition(**kwargs)
                self.load_composition_px(**kwargs)
            return True  # checks that it works
        else:
            return True  # e.g if you just want to check out bulk compositions

    def melting_temperature(self, T_sol=1927, p0=1, **kwargs):
        """ get compositionally-dependent top-of-mantle solidus temperature
         1927 K from Miyazaki & Korenaga (2022), p0 in GPa"""

        def dorn(p, X_Fe):
            # p in GPa, X_Fe in mass fraction
            dT = (102 + 64.1 * p - 3.62 * p ** 2) * (0.1 - X_Fe)
            return dT

        # convert wt% FeO to mass fraction Fe
        m_Fe = self.wt_oxides['FeO'] * M_Fe / M_FeO / 100
        return T_sol + dorn(p0, m_Fe)

    def partition_core_Si(self, x_Si_core):
        # x_Si_core is wt% Si in core

        print('old oxides', self.wt_oxides)
        print('old CMF', self.CMF)
        print('old mg/si', self.mgsi)

        # get total Si
        mass_mtl_old = self.M_p - (self.CMF * self.M_p)  # kg
        mass_Si_tot = self.wt_oxides['SiO2'] * M_Si / M_SiO2 / 100 * mass_mtl_old  # kg, old total

        # boost core mass
        mass_core_old = self.CMF * self.M_p  # kg
        mass_core_new = mass_core_old / (1 - x_Si_core * 1e-2)
        self.CMF = mass_core_new / self.M_p
        print('new CMF', self.CMF)

        mass_Si_core = mass_core_new * x_Si_core * 1e-2  # kg
        mass_mtl_new = self.M_p - (self.CMF * self.M_p)  # kg
        mass_Si_mtl = mass_Si_tot - mass_Si_core  # kg, new

        mass_SiO2_mtl_wtpt = mass_Si_mtl * M_SiO2 / M_Si / mass_mtl_new * 100
        # print('new SiO2 mantle in wt%', mass_SiO2_mtl_wtpt)
        self.wt_oxides['SiO2'] = mass_SiO2_mtl_wtpt

        # renormalise all wt percent oxides
        sum_wt_oxides = sum(self.wt_oxides.values())
        self.wt_oxides = {k: self.wt_oxides[k] / sum_wt_oxides * 100 for k in self.oxide_list}
        self.get_mgsi()
        print('new mg/si', self.mgsi)






def read_from_input(star_name, which='oxide_list', output_path=output_parent_default, oxide_list=None, **kwargs):
    # this isn't tested! trying to read information from perple_x input file
    oldpath = hyp.find_existing_directory(star_name)
    name = os.path.dirname(oldpath)
    oldfile = os.path.join(oldpath, name + '.dat')
    input = []
    with open(oldfile, 'r') as f:
        go = False
        for line in f:
            if which == 'oxide_list':
                if go:
                    s = line.split()
                    for ox in oxide_list.upper():
                        if s[0] == (ox):
                            input.append(s[2])
                if line == 'begin thermodynamic component list':
                    go = True
                if line == 'end thermodynamic component list':
                    go = False
            else:
                raise NotImplementedError('read_from_input(): only oxide list implemented')
    return input
