import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from useful_and_bespoke import colorize, iterable_not_string
from parameters import M_E, M_Fe, M_FeO, M_MgO, M_SiO2, G, R_E, rho_E
import os
import subprocess
import saturation as sat
from bulk_composition import bulk_composition
import time
from scipy.interpolate import splev

perplex_path = '/home/claire/Works/perple_x/'  # path to perple_x installation (https://www.perplex.ethz.ch/)
fig_path = '/home/claire/Works/rocky-water/figs_scratch/'  # path to save figures

# set defaults
wt_oxides_Mars = {'SiO2': 42.71, 'MgO': 31.5, 'CaO': 2.49, 'Al2O3': 2.31,
                  'FeO': 18.71, 'Na2O': 0.5}  # nominal Perple_x SNC meteorite composition, like example 24 but + Na
wt_oxides_Earth = {'SiO2': 44.48, 'MgO': 39.22, 'CaO': 3.44, 'Al2O3': 3.59, 'FeO': 8.10,
                   'Na2O': 0.36}  # weight percent
wt_oxides_pyrolite = {'SiO2': 38.71, 'MgO': 49.85, 'FeO': 6.17, 'CaO': 2.94, 'Al2O3': 2.22, 'Na2O': 0.11}
# solution_phases_default = ['Wus(fab)', 'Pv(fab)', 'O(stx)', 'Wad(stx)', 'Ring(stx)', 'C2/c(stx)', 'Opx(stx)',
#                            'Cpx(stx)', 'Sp(stx)', 'Gt(stx)', 'Aki(fab)']  # like perple_x example 24
oxide_list_default = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'Na2O', 'FeO']
solution_phases_default = ['O', 'Sp', 'Cpx', 'Wad', 'Ring', 'Pv', 'Wus', 'C2/c', 'Opx', 'Aki', 'Ppv', 'Gt_maj', 'CF',
                           'Pl', 'NAl']


def update_MgSi(MgSi=None, oxides=wt_oxides_Earth):
    # vary Mg/Si, use original otherwise
    # update m_MgO and m_SiO2, preserving total mass
    m_tot = oxides['MgO'] + oxides['SiO2']  # total mass (%) of MgO and SiO2 - conserve this
    n_MgO = oxides['MgO'] / M_MgO  # convert to moles
    n_SiO2 = oxides['SiO2'] / M_SiO2
    n_MgSi_old = n_MgO / n_SiO2  # n_Mg = n_MgO etc
    if MgSi is None:
        MgSi = n_MgSi_old  # no update

    m_ratio_new = MgSi * (M_MgO / M_SiO2)  # m = nM, in % | new mass ratio
    x = m_ratio_new
    y = m_tot
    m_MgO_new = x * y / (1 + x)
    m_SiO2_new = m_tot - m_MgO_new
    oxides['MgO'] = m_MgO_new
    oxides['SiO2'] = m_SiO2_new
    return oxides


class PerplexData:
    def __init__(self, name='default', core_efficiency=0.5, M_p=M_E, R_p=None,
                 oxides=oxide_list_default, solution_phases=solution_phases_default,
                 star='sun', perplex_path=perplex_path, verbose=False, **kwargs):
        self.name = name
        self.oxide_list = oxides
        self.solution_phases = solution_phases
        self.star = star
        self.M_p = M_p
        self.R_p = R_p
        self.core_eff = core_efficiency
        self.data_path = perplex_path
        # if verbose:
        #     print('----------------------\ninitialising PerplexData object with M_p = {:.4e}%'.format(M_p), 'kg,',
        #           'core_eff = {:5.2f}%'.format(core_efficiency))

    def load_composition(self, build_file_end='', fillnan=True, check_consistent=True, head=False, verbose=True,
                         **kwargs):
        file = self.data_path + self.name + build_file_end + '_comp.tab'

        df = pd.read_csv(file, skiprows=8, index_col=None, sep=r"\s+",
                         # dtype=np.float64
                         )
        if fillnan:
            df = df.fillna(0)
        if head:
            print(df.head())
        P = df['P(bar)'].to_numpy() * 1e5  # in Pa
        T = df['T(K)'].to_numpy()
        if check_consistent and hasattr(self, 'pressure_m') and (self.pressure_m != P).any():
            raise Exception('current pressure inconsistent with loaded composition data - probably because you were '
                            'using thermodynamic output from an older iteration')
        else:
            self.pressure_m = P
        if check_consistent and hasattr(self, 'temperature_m') and (self.temperature_m != T).any():
            raise Exception('current temperature inconsistent with loaded composition data')
        else:
            self.temperature_m = T
        self.df_comp = df

        phases = df.columns.values.tolist()

        # remove df columns that aren't a mineral phase
        phases.remove('node#')
        phases.remove('P(bar)')
        phases.remove('T(K)')
        self.phases_px = [ph.replace('.', '') for ph in phases]  # not sure why this notation happens sometimes

        if verbose:
            print('-----------------\nphases used in Perple_x:', [ph for ph in self.phases_px])
            print('-----------------')

    def load_adiabat(self, build_file_end='', fillnan=True, check_consistent=True, fix_nan=True, run_bad=False,
                     store=False, head=False, rho_m_av=None, **kwargs):
        file = self.data_path + self.name + build_file_end + '_thermo.tab'

        df = pd.read_csv(file, skiprows=8, index_col=None, sep=r"\s+",
                         # dtype=np.float64
                         )
        if fillnan:
            df = df.fillna(0)
        if head:
            print(df.head())
        P = df['P(bar)'].to_numpy() * 1e5  # in Pa
        T = df['T(K)'].to_numpy()
        alpha_m = df['alpha,1/K'].to_numpy()
        cp_m = df['cp,J/K/kg'].to_numpy()
        rho_m = df['rho,kg/m3'].to_numpy()

        # check if any nans or zeros - flag and for now just extrapolate constant value
        # will want to see how many cases later (limits of Stixrude 2021 dataset e.g.)
        # need to do this for rho, alpha, cp as they could be different idk
        ierr_rho = np.argmax(rho_m == 0)
        ierr_alpha = np.argmax(alpha_m == 0)
        ierr_cp = np.argmax(cp_m == 0)
        if (ierr_cp > 0) or (ierr_rho > 0) or (ierr_alpha > 0):
            arri = np.array([ierr_cp, ierr_rho, ierr_alpha])
            maxp = P[np.min(arri[np.nonzero(arri)])]
            print('warning: NaN in thermodynamic output from p =', maxp*1e-9, 'GPa')
            if run_bad:
                # run again with all output so you can see error - note this causes problems if overwriting
                vertex_command_file = self.name + build_file_end + '_vertex_command.txt'
                with open(self.data_path + vertex_command_file, 'w') as file:
                    s = self.name + build_file_end + '\n0'
                    file.write(s)
                os.system('./vertex < ' + vertex_command_file)
                werami_command_file = self.name + build_file_end + '_werami_command_thermo.txt'
                with open(self.data_path + werami_command_file, 'w') as file:
                    s = self.command_text_thermo(build_file_end=build_file_end)
                    file.write(s)
                os.system('./werami < ' + werami_command_file)
                os.remove(werami_command_file)
                os.remove(vertex_command_file)
                for fend in ['_seismic_data.txt', '_auto_refine.txt', '_1.plt', '.tim', '.plt', '.blk',
                             '.arf', '.tof']:
                    if os.path.isfile(self.data_path + self.name + build_file_end + fend):
                        os.remove(self.data_path + self.name + build_file_end + fend)
            # raise Exception('Perple_x EoS error: cp = 0 at idx', i)
            if fix_nan:
                if ierr_cp > 0:
                    cp_m[ierr_cp:] = cp_m[ierr_cp - 1]
                if ierr_alpha > 0:
                    alpha_m[ierr_alpha:] = alpha_m[ierr_alpha - 1]
                if ierr_rho > 0:
                    rho_m[ierr_rho:] = rho_m[ierr_rho - 1]
                    # use a constant value such that average is preserved from before?
                    # print('old rho_m_av', rho_m_av)
                    # rho_fill = (rho_m_av * len(P) - np.sum(rho_m[:ierr_rho])) / len(rho_m[ierr_rho:])
                    # rho_m[ierr_rho:] = [rho_fill] * len(rho_m[ierr_rho:])
                    # print('new rho_m_av', np.mean(rho_m))

        if store:
            # update PerplexData object
            if check_consistent and hasattr(self, 'pressure_m') and (self.pressure_m != P).any():
                raise Exception('current pressure inconsistent with loaded thermodynamic data')
            else:
                self.pressure_m = P
            if check_consistent and hasattr(self, 'temperature_m') and (self.temperature_m != T).any():
                raise Exception('current temperature inconsistent with loaded thermodynamic data')
            else:
                self.temperature_m = T

            self.alpha_m = alpha_m  # note upside down compared to full planet profile
            self.cp_m = cp_m
            self.rho_m = rho_m
        return rho_m, alpha_m, cp_m, P, T

    def core_mass_fraction(self):
        # doesn't strictly require knowing the planet mass
        # todo: still might want to triple check this is consistent with stellar composition constraints
        # core from leftover iron
        idx = self.oxide_list.index('FEO')
        x_Fe_mm = self.wt_oxides[idx] * (M_Fe / M_FeO)  # mass fraction Fe in mantle wrt mtl
        # print('M Fe/FeO', p.M_Fe / p.M_FeO)
        # print('x_Fe_m wrt mantle', x_Fe_mm)
        x_Fe_cm = -x_Fe_mm * self.core_eff / (self.core_eff - 1)  # mass fraction Fe in core wrt mantle mass
        # print('x_Fe_c wrt mantle', x_Fe_cm)

        # CMF = x_Fe_cm / (x_Fe_cm + 100)

        mantle_mass_fraction = 100 / (100 + x_Fe_cm)
        # print('mantle mass fraction', mantle_mass_fraction)
        # print('core mass fraction', 1 - mantle_mass_fraction)
        x_Fe_c = x_Fe_cm * mantle_mass_fraction  # scale by mass_mtl/M_p --> should be smaller than x_Fe_cm
        # print('x_Fe_c wrt planet', x_Fe_c)
        self.CMF = x_Fe_c / 100

        self.x_Fe_pl = x_Fe_c + x_Fe_mm * mantle_mass_fraction
        # print('x_Fe_pl', x_Fe_pl)

        print('\ncore mass fraction = {:4.2f}\n'.format(self.CMF))
        return self.CMF

    def get_hypatia(self, api_key='136611fd615f8a90aafb1da408f0e5b3',
                    ca_sol=6.33 - 12, al_sol=6.47 - 12, fe_sol=7.45 - 12, si_sol=7.52 - 12, mg_sol=7.54 - 12,
                    na_sol=6.30 - 12):
        """ star id is same as hypatia catalog with spaces e.g. "HIP 12345" """
        import requests

        els = []
        for ox in self.oxide_list:
            els.append(ox[:2].lower())

        if self.star != 'sun':  # not needed for solar values
            params = {"name": [self.star] * len(els), "element": els, "solarnorm": ["lod09"] * len(els)}
            entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", auth=(api_key, "api_token"),
                                 params=params)

            if np.size(entry.json()) == 0:
                raise Exception('No entry found in Hypatia Catalog:', self.star)
            # print('loaded json', entry.json())
        nH_star = []
        for ii, el in enumerate(els):
            # absolute not working for some reason so get difference from solar via lodders norm
            if self.star == 'sun':
                nH = 0
            else:
                try:
                    nH = entry.json()[ii]['mean']
                    print('found', el, entry.json()[ii]['element'], 'mean:', entry.json()[ii]['mean'])
                except IndexError:
                    raise Exception('Catalog does not contain enough elements, try changing list of oxides')

            sol_val = eval(el + '_sol')
            nH_star.append(nH + sol_val)
        self.nH_star = nH_star  # will always be in same order as oxides list
        return nH_star

    def star_to_oxide(self):
        wt_oxides = bulk_composition(self.oxide_list, self.nH_star, self.core_eff)
        self.wt_oxides = wt_oxides
        return wt_oxides

    def write_build(self, build_file_end='', title='Planet', p_min=10000, p_max=245000, adiabat_file='aerotherm.dat',
                    verbose=False, overwrite=True, vertex_data='stx21ver', option_file='perplex_option_new',
                    excluded_phases=None, **kwargs):
        """ write perple_x build file, p in bar but have fucked it before and np """
        if excluded_phases is None:
            excluded_phases = []  # no excluded phases / solution end members
        if verbose:
            print('using thermodynamic dataset:', vertex_data)
        if 'stx21' in vertex_data:
            solution_model = 'stx21_solution_model'
        elif 'stx11' in vertex_data:
            solution_model = 'stx11_solution_model'
        else:
            solution_model = 'solution_model'

        build_file = self.data_path + self.name + build_file_end + '.dat'
        if os.path.isfile(build_file) and not overwrite:
            raise Exception('WARNING: build file', build_file, 'already exists, set overwrite=True')
        elif os.path.isfile(build_file):
            print('  overwriting', build_file)

        s = ''
        with open(build_file, 'w') as file:
            # s = s + 'sfo05ver.dat     thermodynamic data file\n'
            s = s + vertex_data + '.dat     thermodynamic data file\n'
            s = s + 'no_print | print generates print output\n'
            s = s + 'plot     | obsolete 6.8.4+\n'
            s = s + solution_model + '.dat     | solution model file, blank = none\n'
            s = s + title + '\n'
            s = s + option_file + '.dat | Perple_X option file\n'
            s = s + '   10 calculation type: 0- composition, 1- Schreinemakers, 3- Mixed, 4- swash, 5- gridded min, 7- 1d fract, 8- gwash, 9- 2d fract, 10- 7 w/file input, 11- 9 w/file input, 12- 0d infiltration\n'
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
                s = s + el.upper().ljust(6) + '1'.ljust(3) + "{value:{width}.{precision}f}".format(value=float(wt),
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
            for sol in self.solution_phases:
                s = s + '\n' + sol
            s = s + '\nend solution phase list\n\n'

            try:
                s = s + "   {value:{width}.{precision}f}".format(value=float(p_max), width=7,
                                                                 precision=7 - len(str(int(p_max))))[
                        :-1] + '        0.00000        0.00000        0.00000        0.00000     max p, t, xco2, mu_1, mu_2\n'
            except ValueError as e:  # not enough digits
                s = s + "   {value:{width}.{precision}f}".format(value=float(p_max), width=10,
                                                                 precision=10 - len(str(int(p_max))))[
                        :-1] + '        0.00000        0.00000        0.00000        0.00000     max p, t, xco2, mu_1, mu_2\n'
            s = s + "   {value:{width}.{precision}f}".format(value=float(p_min), width=7,
                                                             precision=7 - len(str(int(p_min))))[
                    :-1] + '        0.00000        0.00000        0.00000        0.00000     min p, t, xco2, mu_1, mu_2\n'
            s = s + '   0.00000        0.00000        0.00000        0.00000        0.00000     unused place holder post 06\n\n'
            s = s + ' 1  2  4  5  3   indices of 1st & 2nd independent & sectioning variables'

            file.write(s)
        if verbose:
            print('  wrote to', build_file)

    def write_adiabat(self, P, T, file_end='_adiabat', verbose=False, overwrite=True, **kwargs):
        """ write perple_x build file, p in bar """
        adiabat_file = self.data_path + self.name + file_end + '.dat'
        if os.path.isfile(adiabat_file) and not overwrite:
            raise Exception('WARNING: build file', adiabat_file, 'already exists, set overwrite=True')
        elif os.path.isfile(adiabat_file):
            print('  overwriting', adiabat_file)

        s = ''
        with open(adiabat_file, 'w') as file:
            for p, t in zip(P, T):
                s = s + "{:.6e}".format(p) + '	' + "{:.6e}".format(t) + '\n'
            file.write(s)
        if verbose:
            print('  wrote to', adiabat_file)

    def command_text_composition(self, build_file_end=''):
        # string for werami command file to get compositional data
        s = self.name + build_file_end + '\n'  # Enter the project name (the name assigned in BUILD)
        s = s + '3\n'  # Select operational mode: 3 - properties along a 1d path
        s = s + '25\n'  # Select a property: 25 - Modes of all phases
        s = s + 'n\n'  # Output cumulative modes (y/n)?
        s = s + '0\n'  # Select operational mode: 0 - EXIT
        return s

    def command_text_thermo(self, build_file_end=''):
        # string for werami command file to get thermodynamic data
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
                    verbose=True, werami_command_text_fn=None, output_file_end='.tab', **kwargs):
        os.chdir(self.data_path)
        if suppress_output:
            stderr, stdout = subprocess.DEVNULL, subprocess.DEVNULL
        else:
            stderr, stdout = None, None

        # create vertex command file
        vertex_command_file = self.name + build_file_end + '_vertex_command.txt'
        with open(self.data_path + vertex_command_file, 'w') as file:
            s = self.name + build_file_end + '\n0'
            file.write(s)

        # run vertex
        output = subprocess.run('./vertex < ' + vertex_command_file, shell=True,
                                stdout=stdout, stderr=stderr)
        if verbose:
            print('  ', output)

        # create werami command file
        werami_command_file = self.name + build_file_end + werami_command_end
        with open(self.data_path + werami_command_file, 'w') as file:
            s = werami_command_text_fn(build_file_end=build_file_end)
            file.write(s)

        # run werami
        output = subprocess.run('./werami < ' + werami_command_file, shell=True,
                                stdout=stdout, stderr=stderr)
        if verbose:
            print('  ', output)

        # delete extra files and rename .tab file to something meaningful
        try:
            os.rename(self.data_path + self.name + build_file_end + '_1.tab',
                  self.data_path + self.name + build_file_end + output_file_end)
        except FileNotFoundError as e:
            # something probably went wrong with vertex or werami, run again with full output
            os.system('./vertex < ' + vertex_command_file)
            os.system('./werami < ' + werami_command_file)
            raise e

        if clean:
            for fend in ['_seismic_data.txt', '_auto_refine.txt', '_1.plt', '.tim', '.plt', '.blk', '.arf', '.tof']:
                if os.path.isfile(self.data_path + self.name + build_file_end + fend):
                    os.remove(self.data_path + self.name + build_file_end + fend)
            os.remove(werami_command_file)
            os.remove(vertex_command_file)

    def get_adiabat(self, **kwargs):
        self.run_perplex(werami_command_end='_werami_command_thermo.txt',
                         werami_command_text_fn=self.command_text_thermo,
                         output_file_end='_thermo.tab', **kwargs)

    def get_composition(self, **kwargs):
        self.run_perplex(werami_command_end='_werami_command_comp.txt',
                         werami_command_text_fn=self.command_text_composition,
                         output_file_end='_comp.tab', **kwargs)

    def iterate_structure(self, Psurf=1000, Tsurf=1600, n=500, maxIter=100, tol=0.0001, clean=True,
                          deltaT_cmb=0, rho_m0=None, profile='adiabat', **kwargs):
        """ tweaked from Noack+
        Tsurf is potential surface temperature in K (was 1600), Psurf is surface pressure in bar
        n is radial resolution from center of planet to surface, must be <~ 2000 else perple_x error
        tol is fractional convergence criterion for iteration solver
        rho_m0 is initial guess for mantle density (optional)"""
        import eos
        import math

        M = self.M_p  # planet mass in kg
        Mc = M * self.CMF
        try:
            # x_Fe = self.x_Fe_pl / 100  # wt% planet iron content (here for now without iron in mantle)
            x_Fe = self.CMF
        except AttributeError:
            raise Exception('must run star_to_oxide() and core_mass_fraction() to get planet Fe fraction')

        if n == 'auto':
            if M / M_E <= 2:
                n = 500
            elif M / M_E <= 4:
                n = 1000
            elif M / M_E <= 5:
                n = 1500
        if profile == 'warm':
            Tsurf = self.melting_temperature(p_GPa=Psurf*1e5*1e-9)  # Noack "warm" case

        # Initialization - guesses
        # rho_c_av = 11000  # guess for core density in kg/m^3
        # rho_m_av = 4000  # guess for mantle density in kg/m^3
        cp_c = 800  # guess for core heat capacity in J/kg K
        cp_m = 1300  # guess for mantle heat capacity in J/kg K
        alpha_c = 0.00001  # guess for core thermal expansion ceoff. in 1/K
        alpha_m = 0.000025  # guess for mantle thermal expansion ceoff. in 1/K
        Rp = 1e3 * (7030 - 1840*x_Fe) * (M / M_E)**0.282  # initial guesss, Noack & Lasbleis 2020 (5) ignoring mantle Fe
        if self.CMF > 0:
            Rc = 1e3 * 4850 * x_Fe ** 0.328 * (M / M_E) ** 0.266  # initial guess, hot case, ibid. (9)
            rho_c_av = x_Fe * M / (4/3 * np.pi * Rc**3)
        else:
            Rc = 0
            rho_c_av = 0
        if rho_m0 is None:
            rho_m_av = (1 - x_Fe) * M / (4/3 * np.pi * (Rp**3 - Rc**3))  # Noack & Lasbleis parameterisation
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
        if x_Fe == 1:  # pure iron shell
            density[-1] = rho_c_av
            alpha[-1] = alpha_c
            cp[-1] = cp_c
            Rc = Rp

        # Initialization of arrays: interpolation over depth from top
        for i in range(2, n + 1):  # goes to i = n-1
            radius[n - i] = radius[n - i + 1] - Rp / n  # such that it never goes to 0 (would cause numerical errors)
        if self.CMF > 0:
            i_cmb = np.argmax(radius > Rc)  # index of cmb in profiles
        else:
            i_cmb = 0
        density[i_cmb + 1:] = rho_m_av
        alpha[i_cmb + 1:] = alpha_m
        cp[i_cmb + 1:] = cp_m
        density[:i_cmb + 1] = rho_c_av
        alpha[:i_cmb + 1] = alpha_c
        cp[:i_cmb + 1] = cp_c
        gravity = eos.g_profile(n, radius, density)
        pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, Psurf, Tsurf, i_cmb, deltaT_cmb)
        p_cmb = pressure[i_cmb]  # the pressure (at cell top edge) of the cell that Rc passes through
        p_cmb_guess = np.mean(gravity[i_cmb + 1:]) * rho_m_av * (Rp - Rc)
        p_mantle_bar = pressure[i_cmb + 1:][::-1] * 1e-5  # convert Pa to bar and invert
        T_mantle = temperature[i_cmb + 1:][::-1]
        print('initial p_cen =', pressure[0]*1e-9, 'GPa, p_cmb =', p_cmb*1e-9, 'GPa, p_cmb guess', p_cmb_guess*1e-9,
              'GPa, Rp =', Rp/R_E, 'R_E')

        # Iteration
        print('>>>>>>>>>\nIterating interior structure...')
        it = 1
        iter_param_old = 1e-5
        iter_param = p_cmb  # Rp
        while (abs((iter_param - iter_param_old) / iter_param_old) > tol) and (it < maxIter):
            #     print(it, 'criterion:', abs((Rp - Rp_old)/Rp_old), '>', tol)
            # store old value to determine convergence
            iter_param_old = iter_param

            # average values will be re-determined from material properties
            rho_c = 0
            rho_m = 0
            vol_c = 0
            vol_m = 0

            run_flag = True
            for i in range(n):  # M: 1:n, P: 0:n-1  index from centre to surface
                if i <= i_cmb: # get local thermodynamic properties - core - layer by layer
                    _, density[i], alpha[i], cp[i] = eos.EOS_all(pressure[i] * 1e-9, temperature[i], 4)
                    if cp[i] == 0:
                        print('i', i, 'cp[i]', cp[i], 'problem with core EoS')
                        raise ZeroDivisionError
                else:  # run perple_x to get mantle EOS
                    if run_flag:
                        # run perple_x over mantle p, T to calculate thermodynamic properties (inverse adiabat idx)
                        self.write_adiabat(p_mantle_bar, T_mantle, file_end='_temp' + str(it) + '_adiabat',
                                           **kwargs)
                        self.write_build(build_file_end='_temp' + str(it),
                                         p_min=np.min(p_mantle_bar), p_max=np.max(p_mantle_bar),
                                         adiabat_file=self.name + '_temp' + str(it) + '_adiabat.dat',
                                         **kwargs)
                        self.get_adiabat(build_file_end='_temp' + str(it), clean=clean, **kwargs)

                        # then after running vertex & werami, extract density, alpha, cp
                        density_wer, alpha_wer, cp_wer, p_wer, T_wer = self.load_adiabat(
                            build_file_end='_temp' + str(it),
                            head=False, check_consistent=False, store=False, **kwargs)

                        alpha[i_cmb + 1:] = alpha_wer[::-1]
                        cp[i_cmb + 1:] = cp_wer[::-1]
                        density[i_cmb + 1:] = density_wer[
                                              ::-1]  # invert again because perple_x idx in opposite directions
                        run_flag = False  # only run perple_x once per profile
                        print('just loaded density from werami, top:', density[-1], 'mantle base', density[i_cmb + 1],
                              'cmb', density[i_cmb])
#
                # Get mass vector - not used at the moment though
                # TODO: confused by equation
                if i == 0:
                    mass[i] = 4 * math.pi * radius[i] ** 3 * density[i]  # originally no 3 in denominator
                else:
                    # mass[i] = mass[i - 1] + 4 / 3 * math.pi * (radius[i] ** 3 - radius[i - 1] ** 3) * density[i]
                    mass[i] = mass[i - 1] + (radius[i] - radius[i - 1]) * 4 * math.pi * radius[i] ** 2 * density[i]  # ??

                # Get new average densities for radius update
                if i == 0:
                    vol = 4 / 3 * math.pi * radius[i] ** 3
                else:
                    vol = 4 / 3 * math.pi * (radius[i] ** 3 - radius[i - 1] ** 3)

                if radius[i] <= Rc:
                    rho_c = rho_c + density[i] * vol
                    vol_c = vol_c + vol
                else:
                    rho_m = rho_m + density[i] * vol
                    vol_m = vol_m + vol

            # volume-averaged densities
            if self.CMF > 0:
                rho_c_av = rho_c / vol_c
            else:
                rho_c_av = 0
            rho_m_av = rho_m / vol_m

            # update radius for core and mantle using averages
            Rc, Rp = eos.core_planet_radii(x_Fe, M, rho_c_av, rho_m_av)

            # update vectors for radius, gravity, pressure and temperature
            radius[-1] = Rp
            for i in range(2, n + 1):  # M: 1:n-1; n-i from n-1...1; P: from n-2...0 -> i from 2...n
                radius[n - i] = radius[n - i + 1] - Rp / n  # always a constant divison of Rp

            gravity = eos.g_profile(n, radius, density)

            # pressure and temperature are interpolated from surface downwards
            pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, Psurf, Tsurf, i_cmb,
                                                   deltaT_cmb)

            i_cmb = np.argmax(radius > Rc)   # update index of top of core in profiles
            # i_cmb = np.argmax(mass > Mc)

            if profile == 'solidus':
                # test with melting temperature instead
                temperature[i_cmb + 1:] = self.melting_temperature(p_GPa=pressure[i_cmb + 1:] * 1e-9)

            p_mantle_bar = pressure[i_cmb + 1:][::-1] * 1e-5  # convert Pa to bar
            T_mantle = temperature[i_cmb + 1:][::-1]

            # update parameter that should be converging
            iter_param = pressure[i_cmb]  # Rp
            print(it, "R_p = {:.6e}".format(Rp / R_E), 'R_E', "R_c = {:.4e}".format(Rc / R_E), 'R_E',
                  "p_cmb = {:.5e}".format(pressure[i_cmb] * 1e-9), 'GPa',
                  # "rho_c_av = {:.2e}".format(rho_c_av), 'kg/m3', "rho_m_av = {:.2e}".format(rho_m_av), 'kg/m3',
                  "CMF*M_p = {:.3e}".format(Mc / M_E), 'M_E',
                  'mass[i_cmb] = {:.5e}'.format(mass[i_cmb] / M_E), 'M_E', 'i_cmb', i_cmb)

            if clean:
                for fend in ['.dat', '_adiabat.dat', '_thermo.tab']:
                    if os.path.isfile(self.data_path + self.name + '_temp' + str(it) + fend):
                        os.remove(self.data_path + self.name + '_temp' + str(it) + fend)
            it = it + 1
            # end while

        # ultimately only need adiabat for mantle - write this to file
        self.write_adiabat(p_mantle_bar, T_mantle, file_end='_adiabat', **kwargs)

        # update profiles in object attributes - these are centre-of-planet to surface
        self.radius = radius
        self.density = density
        self.gravity = gravity
        self.temperature = temperature
        self.pressure = pressure
        self.alpha = alpha
        self.cp = cp
        self.cum_mass = mass
        self.mass = np.diff(np.insert(mass, 0, 0))
        self.R_p = Rp
        self.R_c = Rc
        self.i_cmb = i_cmb


    # def iterate_structure_trying(self, Psurf=1000, Tsurf=1600, r_res=2e3, maxIter=100, tol=0.001, clean=True,
    #                       deltaT_cmb=0, **kwargs):
    #     """ problem with this is that if array shapes keep changing can't compute density with array multiplication
    #     so"""
    #     import eos
    #     import math
    #
    #     M = self.M_p  # planet mass in kg
    #     try:
    #         # x_Fe = self.x_Fe_pl / 100  # wt% planet iron content (here for now without iron in mantle)
    #         x_Fe = self.CMF
    #         Mc = self.CMF * M
    #     except AttributeError:
    #         raise Exception('must run star_to_oxide() and core_mass_fraction() to get planet Fe fraction')
    #
    #     # Initialization - guesses
    #     # rho_c = 11000  # guess for core density in kg/m^3
    #     # rho_m = 4000  # guess for mantle density in kg/m^3
    #     cp_c = 800  # guess for core heat capacity in J/kg K
    #     cp_m = 1200  # guess for mantle heat capacity in J/kg K
    #     alpha_c = 0.00001  # guess for core thermal expansion ceoff. in 1/K
    #     alpha_m = 0.000025  # guess for mantle thermal expansion ceoff. in 1/K
    #     Rp = 1e3 * (7030 - 1840*x_Fe) * (M/M_E)**0.282  # initial guesss, Noack & Lasbleis 2020 (5) ignoring mantle Fe
    #     Rc = 1e3 * 4850 * x_Fe**0.328 * (M/M_E)**0.266  # initial guess, hot case, ibid. (9)
    #     rho_c = x_Fe * M / (4 / 3 * np.pi * Rc ** 3)
    #     rho_m = (1 - x_Fe) * M / (4/3 * np.pi * (Rp**3 - Rc**3))
    #
    #     # Initialization of arrays: surface values
    #     radius = [Rp]
    #     density = [rho_m]
    #     alpha = [alpha_m]
    #     cp = [cp_m]
    #     mass = [M]
    #     pressure = [psurf*1e5]
    #     gravity = [G * mass[0] / radius[0] ** 2]
    #     temperature = [Tsurf]
    #
    #     # Keep on depleting mass from surface until you reach centre
    #     while radius[0] - r_res > 0:
    #         mass.insert(0, mass[0] - (r_res * 4 * math.pi * radius[0] ** 2 * density[0]))
    #         radius.insert(0, radius[0] - r_res)
    #         if mass[0] > Mc:
    #             # updated mass is in mantle
    #             density.insert(0, rho_m)
    #             alpha.insert(0, alpha_m)
    #             cp.insert(0, cp_m)
    #         else:  # layer is in core
    #             density.insert(0, rho_c)
    #             alpha.insert(0, alpha_c)
    #             cp.insert(0, cp_c)
    #         gravity.insert(0, G * mass[0] / radius[0] ** 2)
    #         pressure.insert(0, pressure[0] - r_res * gravity[0] * density[0])
    #         temperature.insert(0, temperature[0] - r_res * alpha[0] / cp[0] * gravity[0] * temperature[0])
    #
    #     # to arrays
    #     mass = np.array(mass)
    #     radius = np.array(radius)
    #     density = np.array(density)
    #     alpha = np.array(alpha)
    #     cp = np.array(cp)
    #     gravity = np.array(gravity)
    #     temperature = np.array(temperature)
    #     pressure = np.array(pressure)
    #     n = len(radius)
    #     print('radius', radius*1e-3)
    #     i_cmb = np.argmax(mass > Mc) - 1  # index of cmb in profiles
    #     p_cmb = pressure[i_cmb + 1]
    #
    #     # Iteration
    #     print('Iterating interior structure...')
    #     it = 1
    #     iter_param_old = 1e-5
    #     iter_param = p_cmb  # Rp
    #     while (abs((iter_param - iter_param_old) / iter_param_old) > tol) and (it < maxIter):
    #         #     print(it, 'criterion:', abs((Rp - Rp_old)/Rp_old), '>', tol)
    #         # store old value to determine convergence
    #         iter_param_old = iter_param
    #
    #         i_cmb = np.argmax(mass > Mc) - 1  # index of cmb in profiles - i_cmb is still in the core
    #         p_mantle_bar = pressure[i_cmb + 1:] * 1e-5  # convert Pa to bar
    #         T_mantle = temperature[i_cmb + 1:]
    #
    #         # run perple_x to get mantle materials
    #         self.write_adiabat(p_mantle_bar[::-1], T_mantle[::-1], file_end='_temp' + str(it) + '_adiabat',
    #                            **kwargs)
    #         self.write_build(build_file_end='_temp' + str(it),
    #                          p_min=np.min(p_mantle_bar), p_max=np.max(p_mantle_bar),
    #                          adiabat_file=self.name + '_temp' + str(it) + '_adiabat.dat',
    #                          **kwargs)
    #         self.get_adiabat(build_file_end='_temp' + str(it), clean=clean, **kwargs)
    #
    #         # then after running vertex & werami, extract density, alpha, cp - note these are upside down
    #         density_wer, alpha_wer, cp_wer, p_wer, T_wer = self.load_adiabat(
    #             build_file_end='_temp' + str(it),
    #             head=False, check_consistent=False, store=False)
    #
    #         # now start over planet - but stuck on radius??
    #         mass = [M]
    #         pressure = [psurf * 1e5]
    #         gravity = [G * mass[0] / radius[0] ** 2]
    #         temperature = [Tsurf]
    #
    #         # Keep on depleting mass from surface until you reach centre
    #         while radius[0] - r_res > 0:
    #             mass.insert(0, mass[0] - (r_res * 4 * math.pi * radius[0] ** 2 * density[0]))
    #             radius.insert(0, radius[0] - r_res)
    #             if mass[0] > Mc:
    #                 # updated mass is in mantle
    #                 density.insert(0, rho_m)
    #                 alpha.insert(0, alpha_m)
    #                 cp.insert(0, cp_m)
    #             else:  # layer is in core
    #                 density.insert(0, rho_c)
    #                 alpha.insert(0, alpha_c)
    #                 cp.insert(0, cp_c)
    #             gravity.insert(0, G * mass[0] / radius[0] ** 2)
    #             pressure.insert(0, pressure[0] - r_res * gravity[0] * density[0])
    #             temperature.insert(0, temperature[0] - r_res * alpha[0] / cp[0] * gravity[0] * temperature[0])
    #         for i in range(n - 1):  # from cmb to centre keep on adding mass
    #             # get layer
    #             if n - i > i_cmb:
    #                 if run_flag:
    #                     # run perple_x to get mantle EOS, consistent with given adiabat and M, R
    #
    #
    #                     # run perple_x over mantle p, T to calculate thermodynamic properties (inverse adiabat idx)
    #
    #
    #                     density[i_cmb + 1:] = density_wer[
    #                                           ::-1]  # invert again because perple_x idx in opposite directions
    #                     alpha[i_cmb + 1:] = alpha_wer[::-1]
    #                     cp[i_cmb + 1:] = cp_wer[::-1]
    #
    #
    #             else:
    #                 # get local thermodynamic properties - core
    #                 _, density[n - i], alpha[n - i], cp[n - i] = eos.EOS_all(pressure[n - i] * 1e-9, temperature[n - i], 4)
    #                 if cp[n - i] == 0:
    #                     print('n - i', n - i, 'cp[n - i]', cp[n - i])
    #                     raise ZeroDivisionError
    #
    #         # now re-build planet from bottom up given new material properties - may change array length
    #         mass = [0]
    #         radius = [0]
    #         gravity = [0]
    #         while mass[-1] < M:
    #             radius.append(radius[-1] + r_res)
    #             mass.append(mass[-1] + (r_res * 4 * math.pi * radius[-1] ** 2 * density[-1]))
    #             gravity.append(gravity[-1] + r_res * (
    #                     4 * math.pi * G * density[-1] - 2 * gravity[-1] / radius[-1]))
    #
    #         # to arrays
    #         mass = np.array(mass)
    #         radius = np.array(radius)
    #         gravity = np.array(gravity)
    #         n = len(radius)
    #         print('radius', n, radius * 1e-3)
    #         i_cmb = np.argmax(mass > Mc) - 1  # index of cmb in profiles
    #         pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, Psurf, Tsurf, i_cmb,
    #                                                deltaT_cmb)
    #         p_cmb = pressure[i_cmb + 1]
    #         p_mantle_bar = pressure[i_cmb + 1:] * 1e-5  # convert Pa to bar
    #         T_mantle = temperature[i_cmb + 1:]
    #
    #         # update parameter that should be converging
    #         iter_param = p_cmb  # Rp
    #         print(it, "R_p = {:.8e}".format(Rp / 1000), 'km', "R_c = {:.8e}".format(Rc / 1000), 'km',
    #               "p_cmb = {:.2e}".format(pressure[i_cmb + 1] * 1e-9), 'GPa',
    #               # "rho_c_av = {:.2e}".format(rho_c), 'kg/m3', "rho_m_av = {:.2e}".format(rho_m), 'kg/m3',
    #               "CMF*M_p = {:.5e}".format(self.CMF * self.M_p / M_E), 'M_E',
    #               'mass[i_cmb] = {:.5e}'.format(mass[i_cmb] / M_E), 'M_E')
    #
    #         if clean:
    #             for fend in ['.dat', '_adiabat.dat']:  #, '_thermo.tab']:
    #                 if os.path.isfile(self.data_path + self.name + '_temp' + str(it) + fend):
    #                     os.remove(self.data_path + self.name + '_temp' + str(it) + fend)
    #         it = it + 1
    #         # end while
    #
    #     # ultimately only need adiabat for mantle - write this to file
    #     self.write_adiabat(p_mantle_bar[::-1], T_mantle[::-1], file_end='_adiabat', **kwargs)
    #
    #     # update profiles in object attributes - these are centre-of-planet to surface
    #     self.radius = radius
    #     self.density = density
    #     self.gravity = gravity
    #     self.temperature = temperature
    #     self.pressure = pressure
    #     self.alpha = alpha
    #     self.cp = cp
    #     self.cum_mass = mass
    #     self.mass = np.diff(np.insert(mass, 0, 0))
    #     self.R_p = Rp
    #     self.R_c = Rc
    #     self.i_cmb = i_cmb

    def plot_structure(self, fig_path=fig_path, save=True):

        fig, ax = plt.subplots(4, 2, sharex=True, figsize=[12, 8])

        ax[0, 0].plot(self.radius / 1000, self.density, color='black', alpha=1)
        ax[0, 0].set_ylabel(r"$\rho$ (kg m$^{-3}$)")

        ax[1, 0].plot(self.radius / 1000, self.gravity, color='black', alpha=1)
        ax[1, 0].set_ylabel("$g$ (m s$^{-2}$)")

        ax[0, 1].plot(self.radius / 1000, self.pressure / 10 ** 9, color='black', alpha=1)
        ax[0, 1].set_ylabel("$P$ (GPa)")

        ax[1, 1].plot(self.radius / 1000, self.temperature, color='black', alpha=1)
        ax[1, 1].set_ylabel("$T$ (K)")

        ax[2, 0].plot(self.radius / 1000, self.alpha * 10 ** 5, color='black', alpha=1)
        ax[2, 0].set_ylabel(r"$\alpha$ ($10^{-5}$ K$^{-1}$)")

        ax[2, 1].plot(self.radius / 1000, self.cp, color='black', alpha=1)
        ax[2, 1].set_ylabel("$C_p$ (J kg$^{-1}$ K$^{-1}$)")

        ax[3, 0].plot(self.radius / 1000, self.mass, color='black', alpha=1)
        ax[3, 0].set_ylabel("$m$ (kg)")
        ax[3, 0].set_xlabel("radius (km)")

        plt.suptitle(self.name)
        plt.tight_layout()
        if save:
            fig.savefig(fig_path + self.name + '_structure.png', bbox_inches='tight')
        else:
            plt.show()

    def plot_composition(self, which='pressure', comp_stacked=False, fig_path=fig_path, save=True, cmap='tab20',
                         labelsize=16, **kwargs):
        # print('df_comp\n', self.df_comp)
        fig, ax = plt.subplots(1, 1)
        if which == 'pressure':
            x = self.pressure_m * 1e-9  # plot P in GPa
            ax.set_xlabel('Pressure (GPa)', fontsize=labelsize)
        else:
            raise NotImplementedError('independent variables other than pressure not implemented')
        ax.set_ylabel('Modal abundance (mol %)', fontsize=labelsize)
        colours = colorize(range(len(self.phases_px)), cmap=cmap)[0]

        y_stacked = np.zeros_like(x)
        for ii, phase in enumerate(self.phases_px):
            y = self.df_comp[phase]
            if comp_stacked:
                y_stacked = np.vstack((y_stacked, y))
            else:
                ax.plot(x, y, c=colours[ii], label=phase)
        if comp_stacked:
            ax.stackplot(x, y_stacked, labels=self.phases_px, colors=colours)
            ax.set_ylim(0, 100)
        ax.set_xlim(0, np.max(x))
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
        plt.title(self.name, fontsize=labelsize)

        plt.tight_layout()
        if save:
            fig.savefig(fig_path + self.name + '_composition.png', bbox_inches='tight')
        else:
            plt.show()

    def plot_gety(self, yvar, ylabel=None, yscale=None):
        if yvar == 'p':
            y = self.df_all['P(bar)'].to_numpy()
            if ylabel is None:
                ylabel = 'Pressure (GPa)'
            if yscale is None:
                yscale = 1e-4  # bar to GPa
        elif yvar == 'T':
            y = self.df_all['T(K)'].to_numpy()
            if ylabel is None:
                ylabel = 'T (K)'
            if yscale is None:
                yscale = 1
        elif yvar == 'r':
            y = self.df_all['r'].to_numpy()
            if ylabel is None:
                ylabel = 'Radius (km)'
            if yscale is None:
                yscale = 1e-3
        elif yvar == 'z':
            y = self.R_p - self.df_all['r'].to_numpy()
            if ylabel is None:
                ylabel = 'Depth (km)'
            if yscale is None:
                yscale = 1e-3
        return y * yscale, ylabel

    def plot_profile(self, col, yvar='p', xlabel=None, ylabel=None, xscale=1, yscale=None, c='k', lw=2,
                     fig=None, ax=None, logx=False, reverse_y=True, label=None, labelsize=14,
                     y2var=None, y2label=None, y2scale=None, save=True, **kwargs):
        if fig is None and ax is None:
            fig, ax = plt.subplots(1, 1)

        # get data to plot
        x = self.df_all[col].to_numpy() * xscale
        y, ylabel = self.plot_gety(yvar, ylabel, yscale)

        ax.plot(x, y, c=c, lw=lw, label=label, **kwargs)
        ax.set_xlabel(xlabel, fontsize=labelsize)
        ax.set_ylabel(ylabel, fontsize=labelsize)
        ax.set_ylim(min(y), max(y))

        if logx:
            ax.set_xscale('log')
        if reverse_y:
            ax.invert_yaxis()
        if y2var is not None:
            y2, y2label = self.plot_gety(y2var, y2label, y2scale)
            ax2 = ax.twinx()
            ax2.plot(x, y2, c=c, lw=lw, label=label, **kwargs)
            ax2.set_ylabel(y2label, fontsize=labelsize)
            ax2.set_ylim(min(y2), max(y2))
            # if reverse_y:
            #     ax2.invert_yaxis()

        if save:
            plt.tight_layout()
            fig.savefig(fig_path + self.name + '_' + col + '.png', bbox_inches='tight')
        return fig, ax

    def melting_temperature(self, p_GPa=None):
        # Dorn & Lichtenberg 2021 (eq 2-3)
        #  following Belonoshko et al.(2005) and Stixrude (2014)
        if p_GPa is None:
            p_GPa = self.pressure_m * 1e-9  # GPa

        # use your bs spline to peridotite melting, Fiquet+ 2010
        tck = (np.array([31.11111111, 31.11111111, 34.44444444, 38.66666667,
                41.33333333, 44., 47.77777778, 53.77777778,
                58.88888889, 62.88888889, 67.33333333, 71.77777778,
                76.88888889, 82.88888889, 84.66666667, 89.11111111,
                95.77777778, 99.77777778, 104.44444444, 108.44444444,
                115.77777778, 126., 132.88888889, 139.77777778,
                139.77777778]), np.array([2650., 2731.25, 2837.5, 2900., 2962.5, 3043.75, 3168.75,
                                       3268.75, 3350., 3425., 3506.25, 3581.25, 3675., 3700.,
                                       3762.5, 3843.75, 3887.5, 3943.75, 3981.25, 4043.75, 4118.75,
                                       4156.25, 4187.5, 0., 0.]), 1)

        T_melt = splev(p_GPa, tck)

        self.T_melt = T_melt
        return T_melt

    def plot_melting(self, labelsize=14):
        from useful_and_bespoke import dark_background
        p = self.pressure_m*1e-9  # GPa
        T_a = self.temperature_m  # K
        # Fe_num = self.wt_oxides['FeO'] / self.wt_oxides['MgO']
        # T_sol = 1409.15 + 134.2 * p - 6.581 * p ** 2 + 0.1054 * p ** 3 + (102.0 + 64.1 * p - 3.62 * p ** 2) * (0.1 - Fe_num)
        # T_m_stx = 5400 * (p / 140) ** 0.48 / (1 - np.log(1 - Fe_num))
        T_melt = self.melting_temperature()

        # compare to Dorn parameterisation
        i_189 = np.argmax(p >= 189.75)
        x_FeO = self.wt_oxides['FeO'] * 1e-2  # wt fraction FeO in mantle, lowers melting temp
        a1 = 1831
        a2 = 4.6
        a3 = 0.33
        b1 = 5400
        b2 = 140
        b3 = 0.48
        c1 = 360
        c2 = 0.0818
        c3 = 102
        c4 = 64.1
        c5 = 3.62
        if i_189 == 0:  # no pressures above this
            T_melt_Dorn = a1 * (1 + p / a2) ** a3 + c1 * (c2 - x_FeO)
        else:
            T_melt_Dorn[:i_189] = a1 * (1 + p / a2) ** a3 + c1 * (c2 - x_FeO)
            T_melt_Dorn[i_189:] = b1 * (p / b2) ** b3 + (c3 + c4 * p - c5 * p ** 2) * (c2 - x_FeO)

        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        ax.plot(T_a, p, label='Mantle adiabat', c='w', ls='-')
        ax.plot(T_melt, p, label='Peridotite solidus (Fiquet+ 2010)', c='w', ls='--')
        # ax.plot(T_melt_Dorn, p, label='solidus Dorn eq (dry)')
        ax.invert_yaxis()
        ax.set_xlabel('T (K)',  fontsize=labelsize)
        ax.set_ylabel('p (GPa)',  fontsize=labelsize)
        leg = ax.legend(frameon=False, fontsize=12, bbox_to_anchor=(-0.1, 1.02, 1., .102), loc='lower left',)
        fig, *ax = dark_background(fig, ax, )
        fig.savefig(fig_path + 'melt.png', bbox_inches='tight',bbox_extra_artists=(leg,),
                    facecolor=fig.get_facecolor())
        plt.show()



def build_planet(name=None, M_p=M_E, star='sun', core_efficiency=0.8, oxides=oxide_list_default, Tp=1600,
                 get_saturation=True, test_CMF=None, test_oxides=None,
                 plot=False, store_all=True, clean=True, suffix=None, plot_kwargs={}, **kwargs):
    """ kwargs include
     overwrite (bool) : auto overwrite build etc files,
     head (bool) : to print DataFrame headers when reading them
     n (int) : interior structure radial resolution,
     tol (float) : iteration tolerance for R_p in fraction
     maxIter (int): max iterations
     verbose (bool) : lots of places """
    time_start = time.time()

    # name object systematically
    if name is None:
        mass_str = str(int(np.round(M_p/M_E))) + 'M'
        if test_CMF is not None:
            cmf_str = str(int(test_CMF*100)) + 'CMF'
        else:
            cmf_str = str(int(core_efficiency * 100)) + 'Ceff'
        if test_oxides is not None:
            comp_str = str(int(np.round(test_oxides['MgO']))) + 'Mg'
        elif star is not None:
            comp_str = star
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

    # create class object
    dat = PerplexData(name=name, M_p=M_p, star=star, core_efficiency=core_efficiency, oxides=oxides, **kwargs)

    # get oxide bulk composition and CMF (skip if input a testing value)
    if test_oxides is None:
        dat.get_hypatia()
        dat.star_to_oxide()
    else:
        if star not in ('sun', None):
            print('Did you mean to input both star name and fixed oxide wt%s?')
        dat.wt_oxides = {k: test_oxides[k] for k in oxides}
    if test_CMF is None:
        dat.core_mass_fraction()
    else:
        dat.CMF = test_CMF

    # iterate perplex to get interior structure and geotherm
    dat.iterate_structure(Tsurf=Tp, clean=clean, **kwargs)

    # run perple_x again to get mantle compositional profile
    dat.write_build(adiabat_file=dat.name + '_adiabat.dat', p_min=dat.pressure[-1]*1e-5,
                    p_max=dat.pressure[dat.i_cmb + 1]*1e-5, **kwargs)
    dat.get_composition(clean=clean, **kwargs)
    dat.load_composition(**kwargs)

    # rename and scale columns (for use with saturation calcs but being consistent regardless)
    df_all = dat.df_comp.copy()
    for col in dat.phases_px:
        # print('found column:', col)
        try:
            df_all[col] = df_all[col] * 1e-2  # convert from % to frac
            new_col = 'X_' + col.lower()
            if col.endswith('(stx)') or col.endswith('(fab)'):  # edit as needed
                new_col = new_col[:-5]
            # need to rename some more
            if new_col == 'X_o':
                new_col = 'X_ol'
            elif new_col == 'X_c2/c':
                new_col = 'X_hpcpx'
            elif new_col == 'X_ca-pv':
                new_col = 'X_capv'
            df_all.rename(columns={col: new_col}, inplace=True)
        except KeyError as e:
            print('something very wrong:', e)
            print('df columns:', df_all.columns)

    # put important mantle data into single dataframe for convenience - note these are all SI units
    df_all['mass(kg)'] = dat.mass[dat.i_cmb + 1:][::-1]  # invert and just store mantle (surface down)
    df_all['z(m)'] = dat.R_p - dat.radius[dat.i_cmb + 1:][::-1]
    df_all['rho'] = dat.density[dat.i_cmb + 1:][::-1]
    df_all['alpha'] = dat.alpha[dat.i_cmb + 1:][::-1]
    df_all['cp'] = dat.cp[dat.i_cmb + 1:][::-1]

    # calculate saturation water content profile per mineral phase
    if get_saturation:
        df_all = sat.mineral_water_contents(dat.pressure_m, dat.temperature_m, df=df_all)
        # print(df_all.head(20))
        c_h2o_tot = sat.total_water(df_all)
        dat.c_h2o_tot = c_h2o_tot  # weight fraction
        print('\n\ntotal water in mantle:', c_h2o_tot * 1e2, 'wt% =', c_h2o_tot * dat.M_p / sat.TO, 'OM')

    dat.df_all = df_all

    if plot:
        # plot
        dat.plot_structure()
        dat.plot_composition(**kwargs)
        if get_saturation:
            dat.plot_profile(col='c_h2o', xlabel='H$_2$O storage capacity (ppm)', xscale=1e6,
                             yvar='p', y2var=None, logx=True, **plot_kwargs)
    if store_all:
        # make mega df to save all mantle data
        file = dat.data_path + dat.name + '_profiles.dat'
        # df_all = dat.df_all.drop('node#', axis=1, inplace=False)
        df_all.to_csv(path_or_buf=file, sep='\t', na_rep='NaN', index=True, index_label=None)
    if clean:
        # move perplex files to their own folder
        new_dir = os.path.dirname(dat.data_path + 'output/' + dat.name + '/')
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        for fend in ['_comp.tab', '_adiabat.dat', '_profiles.dat', '.dat']:
            file = dat.name + fend
            if os.path.exists(dat.data_path + file):
                os.rename(dat.data_path + file, new_dir + '/' + file)

    time_end = time.time()
    print('\n>>>>>> COMPLETED', dat.name, 'in', time_end - time_start, 's\n\n')
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
