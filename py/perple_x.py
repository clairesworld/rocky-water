import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from useful_and_bespoke import colorize, iterable_not_string
import parameters as p
import os
import subprocess

perplex_path = '/home/claire/Works/perple_x/'  # path to perple_x installation (https://www.perplex.ethz.ch/)
fig_path = '/home/claire/Works/rocky-water/figs_scratch/'  # path to save figures

wt_oxides_Mars = [42.71, 31.5, 2.49, 2.31, 18.71]  # nominal Perple_x SNC meteorite composition, for testing


class PerplexData:
    def __init__(self, name='default', core_efficiency=0.5, M_p=p.M_E, oxides=['SIO2', 'MGO', 'CAO', 'AL2O3', 'FEO'],
                 star='sun', perplex_path=perplex_path, verbose=False, **kwargs):
        self.name = name
        self.oxide_list = oxides
        self.star = star
        self.M_p = M_p  # Holzapfel EoS for hcp iron valid to 2 M_E
        self.core_eff = core_efficiency
        self.data_path = perplex_path
        if verbose:
            print('----------------------\ninitialising PerplexData object with M_p = {:.4e}%'.format(M_p), 'kg,',
                  'core_eff = {:5.2f}%'.format(core_efficiency))

    def load_composition(self, build_file_end='', fillnan=True, check_consistent=True, head=False, **kwargs):
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
        self.phases = phases

    def load_adiabat(self, build_file_end='', fillnan=True, check_consistent=True, store=False, head=False, **kwargs):
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

            self.alpha_m = alpha_m
            self.cp_m = cp_m
            self.rho_m = rho_m
        return rho_m, alpha_m, cp_m, P, T

    def core_mass_fraction(self):
        # todo: still might want to triple check this is consistent with stellar composition constraints
        # core from leftover iron
        idx = self.oxide_list.index('FEO')
        x_Fe_mm = self.wt_oxides[idx] * (p.M_Fe / p.M_FeO)  # mass fraction Fe in mantle wrt mtl
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

    def get_hypatia(self, api_key='c53fd88b7e9c7e0c2e719ea4ea3c5e46',
                    ca_sol=6.33 - 12, al_sol=6.47 - 12, fe_sol=7.45 - 12, si_sol=7.52 - 12, mg_sol=7.54 - 12):
        """ star id is same as hypatia catalog with spaces e.g. "HIP 12345" """
        import requests

        els = []
        for ox in self.oxide_list:
            els.append(ox[:2].lower())

        if self.star != 'sun':  # not needed for solar values
            params = {"name": [self.star] * len(els), "element": els, "solarnorm": ["lod09"] * len(els)}
            entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", auth=(api_key, "api_token"),
                                 params=params)

        nH_star = []
        for ii, el in enumerate(els):
            # absolute not working for some reason so get difference from solar via lodders norm
            if self.star == 'sun':
                nH = 0
            else:
                nH = entry.json()[ii]['mean']
                # print(el, entry.json()[ii]['element'], entry.json()[ii]['mean'])
            sol_val = eval(el + '_sol')
            nH_star.append(nH + sol_val)
        self.nH_star = nH_star
        return nH_star

    def star_to_oxide(self, M_oxides=[p.M_SiO2, p.M_MgO, p.M_CaO, p.M_Al2O3, p.M_FeO]):
        """ nH_star is list of the absolute log(N_X/N_H) for the star, in same order
         would need to convert from wrt solar if necessary
         by convention base element (denomenator for oxide wt ratios) is first in list
         core_efficiency is fraction of moles Fe that go into core instead of mantle FeO wrt total moles Fe """

        wt_oxides = [1]  # give denomenator 1 for now
        for ii, ox in enumerate(self.oxide_list):
            if ii > 0:
                if ox == 'AL2O3':
                    X_ratio_mol = 0.5 * 10 ** self.nH_star[ii] / 10 ** self.nH_star[0]  # 2 mols Al per 1 mol Al2O3
                else:
                    X_ratio_mol = 10 ** self.nH_star[ii] / 10 ** self.nH_star[0]  # cancel out H abundance

                if ox == 'FEO':
                    # some molar percentage of Fe goes to core instead of to mantle FeO
                    X_ratio_mol = X_ratio_mol * (1 - self.core_eff)

                # convert from mole ratio to mass ratio assuming all metals in oxides
                M_ratio = M_oxides[ii] / M_oxides[0]
                m_ratio = X_ratio_mol * M_ratio
                wt_oxides.append(m_ratio)

        # now have ratios in wt% 1:X1:X2:X3... normalise so total wt% is 100
        wt_oxides = np.array(wt_oxides)
        wt_oxides = wt_oxides / sum(wt_oxides) * 100

        print('wt.% oxides\n-----------')
        for chem, val in zip(self.oxide_list, wt_oxides):
            print("{0:<7}".format(chem), "{:5.2f}%".format(val))

        self.wt_oxides = wt_oxides
        return wt_oxides

    def write_build(self, build_file_end='', title='Planet', p_min=10000, p_max=245000, adiabat_file='aerotherm.dat',
                    verbose=True, overwrite=False, **kwargs):
        """ write perple_x build file, p in bar """
        build_file = self.data_path + self.name + build_file_end + '.dat'
        if os.path.isfile(build_file) and not overwrite:
            raise Exception('WARNING: build file', build_file, 'already exists, set overwrite=True')
        elif os.path.isfile(build_file):
            print('  overwriting', build_file)

        s = ''
        with open(build_file, 'w') as file:
            s = s + 'sfo05ver.dat     thermodynamic data file\n'
            s = s + 'no_print | print generates print output\n'
            s = s + 'plot     | obsolete 6.8.4+\n'
            s = s + 'solution_model.dat     | solution model file, blank = none\n'
            s = s + title + '\n'
            s = s + 'perplex_option.dat | Perple_X option file\n'
            s = s + '   10 calculation type: 0- composition, 1- Schreinemakers, 3- Mixed, 4- swash, 5- gridded min, 7- 1d fract, 8- gwash, 9- 2d fract, 10- 7 w/file input, 11- 9 w/file input, 12- 0d infiltration\n'
            s = s + adiabat_file + '     | coordinate file \n'
            s = s + '    0 unused place holder, post 06\n' * 9
            s = s + '    0 number component transformations\n'
            s = s + '    5 number of components in the data base\n'
            s = s + '    1 component amounts, 0 - mole, 1 mass\n'
            s = s + '    0 unused place holder, post 06\n' * 2
            s = s + '    0 unused place holder, post 05\n'
            s = s + '    5 ifug EoS for saturated phase\n'
            s = s + '    2 gridded minimization dimension (1 or 2)\n'
            s = s + '    0 special dependencies: 0 - P and T independent, 1 - P(T), 2 - T(P)\n'
            s = s + ' 0.00000      0.00000      0.00000      0.00000      0.00000     Geothermal gradient polynomial coeffs.\n\n'

            s = s + 'begin thermodynamic component list\n'
            for ii, el in enumerate(self.oxide_list):
                wt = self.wt_oxides[ii]
                dig = len(str(int(wt)))
                s = s + el.ljust(6) + '1'.ljust(3) + "{value:{width}.{precision}f}".format(value=float(wt), width=7,
                                                                                           precision=6 - dig)
                s = s + '      0.00000      0.00000     mass  amount\n'
            s = s + 'end thermodynamic component list\n\n\n'

            s = s + 'begin saturated component list\nend saturated component list\n\n\n'
            s = s + 'begin saturated phase component list\nend saturated phase component list\n\n\n'
            s = s + 'begin independent potential/fugacity/activity list\nend independent potential list\n\n\n'
            s = s + 'begin excluded phase list\nstv\nend excluded phase list\n\n\n'
            s = s + 'begin solution phase list\nWus(fab)\nPv(fab)\nO(stx)\nWad(stx)\nRing(stx)\nC2/c(stx)\nOpx(stx)\nCpx(stx)\nSp(stx)\nGt(stx)\nAki(fab)\nend solution phase list\n\n'

            s = s + "   {value:{width}.{precision}f}".format(value=float(p_max), width=7,
                                                             precision=7 - len(str(int(p_max))))[
                    :-1] + '        0.00000        0.00000        0.00000        0.00000     max p, t, xco2, mu_1, mu_2\n'
            s = s + "   {value:{width}.{precision}f}".format(value=float(p_min), width=7,
                                                             precision=7 - len(str(int(p_min))))[
                    :-1] + '        0.00000        0.00000        0.00000        0.00000     min p, t, xco2, mu_1, mu_2\n'
            s = s + '   0.00000        0.00000        0.00000        0.00000        0.00000     unused place holder post 06\n\n'
            s = s + ' 1  2  4  5  3   indices of 1st & 2nd independent & sectioning variables'

            file.write(s)
        if verbose:
            print('  wrote to', build_file)

    def write_adiabat(self, P, T, file_end='_adiabat', verbose=True, overwrite=False, **kwargs):
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
        # os.system('./vertex < ' + self.name + build_file_end + '_vertex_command.txt')

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
        # os.system('./werami < ' + self.name + build_file_end + '_werami_therm_command.txt')

        # delete extra files and rename .tab file to something meaningful
        os.rename(self.data_path + self.name + build_file_end + '_1.tab',
                  self.data_path + self.name + build_file_end + output_file_end)
        if clean:
            for fend in ['_seismic_data.txt', '_auto_refine.txt', '_1.plt', '.tof', '.tim', '.plt', '.blk', '.arf']:
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

    def iterate_structure(self, Psurf=1000, Tsurf=1300, n=200, maxIter=100, tol=10, clean=True, **kwargs):
        """ tweaked from Lena Noack - todo find citation
        Tsurf is potential surface temperature in K (was 1600), Psurf is surface pressure in bar
        n is radial resolution from center of planet to surface"""
        import eos
        import math

        M = self.M_p  # planet mass in kg
        try:
            # x_Fe = self.x_Fe_pl / 100  # wt% planet iron content (here for now without iron in mantle)
            x_Fe = self.CMF
        except AttributeError:
            raise Exception('must run star_to_oxide() and core_mass_fraction() to get planet Fe fraction')

        # Initialization - guesses
        rho_c = 11000  # guess for core density in kg/m^3
        rho_m = 4000  # guess for mantle density in kg/m^3
        cp_c = 800  # guess for core heat capacity in J/kg K
        cp_m = 1300  # guess for mantle heat capacity in J/kg K
        alpha_c = 0.00001  # guess for core thermal expansion ceoff. in 1/K
        alpha_m = 0.000025  # guess for mantle thermal expansion ceoff. in 1/K
        Rc, Rp = eos.core_planet_radii(x_Fe, M, rho_c, rho_m)

        # Arrays
        radius = np.zeros(n)  # array goes from 0 to 199
        density = np.zeros(n)
        alpha = np.zeros(n)
        cp = np.zeros(n)
        mass = np.zeros(n)

        # Initialization of arrays: surface values
        radius[-1] = Rp
        density[-1] = rho_m
        alpha[-1] = alpha_m
        cp[-1] = cp_m
        if x_Fe == 1:  # pure iron shell
            density[-1] = rho_c
            alpha[-1] = alpha_c
            cp[-1] = cp_c
            Rc = Rp

        # Initialization of arrays: interpolation over depth
        for i in range(2, n + 1):  # goes to i = n-1
            radius[n - i] = radius[n - i + 1] - Rp / n
            if radius[n - i] > Rc:
                density[n - i] = rho_m
                alpha[n - i] = alpha_m
                cp[n - i] = cp_m
            else:
                density[n - i] = rho_c
                alpha[n - i] = alpha_c
                cp[n - i] = cp_c  ## goes to idx = n - i = n - n

        gravity = eos.g_profile(n, radius, density)
        pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, Psurf, Tsurf)

        # Iteration
        print('Iterating interior structure...')
        it = 1
        Rp_old = 0
        # Rp_store = []
        # TODO: converge on p_cen or p_cmb instead? shouldn't change results
        while (abs(Rp - Rp_old) > tol) and (it < maxIter):
            # store old Rp value to determine convergence
            # Rp_store.append(Rp * 1e-3)
            Rp_old = Rp

            # average values will be re-determined from material properties
            rho_c = 0
            rho_m = 0
            vol_c = 0
            vol_m = 0

            i_cmb = np.argmax(radius > Rc) - 1  # index of cmb in profiles
            run_flag = True
            for i in range(n):  # M: 1:n, P: 0:n-1
                # get layer
                if i <= i_cmb:
                    # get local thermodynamic properties - core
                    _, density[i], alpha[i], cp[i] = eos.EOS_all(pressure[i] * 1e-9, temperature[i], 4)
                    if cp[i] == 0:
                        print('i', i, 'cp[i]', cp[i])
                        raise ZeroDivisionError
                else:
                    if run_flag:
                        p_mantle_bar = pressure[i_cmb + 1:] * 1e-5  # convert Pa to bar
                        T_mantle = temperature[i_cmb + 1:]

                        # run perple_x over mantle p, T to calculate thermodynamic properties
                        self.write_adiabat(p_mantle_bar[::-1], T_mantle[::-1], file_end='_temp' + str(it) + '_adiabat',
                                           **kwargs)
                        self.write_build(build_file_end='_temp' + str(it),
                                         p_min=np.min(p_mantle_bar), p_max=np.max(p_mantle_bar),
                                         adiabat_file=self.name + '_temp' + str(it) + '_adiabat.dat',
                                         **kwargs)
                        self.get_adiabat(build_file_end='_temp' + str(it), clean=clean, **kwargs)

                        # then after running vertex & werami, extract density, alpha, cp
                        density_wer, alpha_wer, cp_wer, p_wer, T_wer = self.load_adiabat(
                            build_file_end='_temp' + str(it),
                            head=False, check_consistent=False, store=False)

                        density[i:] = density_wer[::-1]
                        alpha[i:] = alpha_wer[::-1]
                        cp[i:] = cp_wer[::-1]
                        run_flag = False  # only run perple_x once per profile

                    if cp[i] == 0:
                        # Perple_x EoS error
                        cp[i] = cp[i+1]
                        print('warning: correcting cp = 0 at idx', i)

                # get mass vector
                if i == 0:
                    mass[i] = 4 * math.pi * radius[i] ** 3 * density[i]
                else:
                    mass[i] = mass[i - 1] + (radius[i] - radius[i - 1]) * 4 * math.pi * radius[i] ** 2 * density[i]

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

                # for ii, rr in zip(range(n), density):
                #     print(ii, rr)

            # volume-averaged densities
            rho_c = rho_c / vol_c
            rho_m = rho_m / vol_m

            # update radius for core and mantle
            Rc, Rp = eos.core_planet_radii(x_Fe, M, rho_c, rho_m)

            # update vectors for radius, gravity, pressure and temperature
            radius[-1] = Rp
            for i in range(2, n + 1):  # M: 1:n-1; n-i from n-1...1; P: from n-2...0 -> i from 2...n
                radius[n - i] = radius[n - i + 1] - Rp / n  # always a constant divison of Rp

            gravity = eos.g_profile(n, radius, density)

            # pressure and temperature are interpolated from surface downwards
            pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, Psurf, Tsurf)

            i_cmb = np.argmax(radius > Rc) - 1  # update index of top of core in profiles
            p_mantle_bar = pressure[i_cmb + 1:] * 1e-5  # convert Pa to bar
            T_mantle = temperature[i_cmb + 1:]

            print(it, "R_p = {:.8e}".format(Rp / 1000), 'km', "R_c = {:.8e}".format(Rc / 1000), 'km',
                  "p_cmb = {:.2e}".format(pressure[i_cmb] * 1e-9), 'GPa',
                  "rho_c_av = {:.2e}".format(rho_c), 'kg/m3', "rho_m_av = {:.2e}".format(rho_m), 'kg/m3',
                  "CMF*M_p = {:.5e}".format(self.CMF * self.M_p / p.M_E), 'M_E',
                  'mass[i_cmb] = {:.5e}'.format(mass[i_cmb] / p.M_E), 'M_E')

            if clean:
                for fend in ['.dat', '_adiabat.dat', '_thermo.tab']:
                    if os.path.isfile(self.data_path + self.name + '_temp' + str(it) + fend):
                        os.remove(self.data_path + self.name + '_temp' + str(it) + fend)

            it = it + 1

            # if it == 2:
            #     print('did one iter')
            #     return None

        # ultimately only need adiabat for mantle - write this to file
        self.write_adiabat(p_mantle_bar[::-1], T_mantle[::-1], file_end='_adiabat', **kwargs)

        # update profiles in object attributes - these are centre-of-planet to surface
        self.radius = radius
        self.density = density
        self.gravity = gravity
        self.temperature = temperature
        self.pressure = pressure
        self.alpha = alpha
        self.cp = cp
        self.mass = mass
        self.R_p = Rp
        self.R_c = Rc

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

    def plot_composition(self, which='pressure', fig_path=fig_path, save=True, **kwargs):

        fig, ax = plt.subplots(1, 1)
        if which == 'pressure':
            x = self.pressure_m * 1e-9  # plot P in GPa
            ax.set_xlabel('Pressure (GPa)')
        else:
            raise NotImplementedError('independent variables other than pressure not implemented')

        ax.set_ylabel('Modal abundance (mol %)')
        colors = colorize(range(len(self.phases)), cmap='rainbow')[0]
        for ii, phase in enumerate(self.phases):
            y = self.df_comp[phase]
            plt.plot(x, y, c=colors[ii], label=phase)
        ax.legend()
        plt.title(self.name)

        plt.tight_layout()
        if save:
            fig.savefig(fig_path + self.name + '_composition.png', bbox_inches='tight')
        else:
            plt.show()


def start(name=None, M_p=p.M_E, stars='sun', core_efficiency=0.8, clean=True, **kwargs):
    """ kwargs include
     overwrite (bool) : auto overwrite build etc files,
     head (bool) : to print DataFrame headers when reading them
     n (int) : interior structure radial resolution,
     tol (float) : iteration tolerance for R_p
     verbose (bool) : lots of places """
    if isinstance(stars, str):
        stars = [stars]
    if name is None:
        # name after star by default
        name = [s.replace(" ", "") for s in stars]
    elif isinstance(name, str):
        name = [name]  # just one star given

    for ii, star in enumerate(stars):
        dat = PerplexData(name=name[ii], M_p=M_p, star=star, core_efficiency=core_efficiency, **kwargs)

        # get interior structure and geotherm
        dat.get_hypatia()
        dat.star_to_oxide()
        dat.core_mass_fraction()
        dat.iterate_structure(clean=clean, **kwargs)

        # get composition
        dat.write_build(adiabat_file=dat.name + '_adiabat.dat', **kwargs)
        dat.get_composition(clean=clean, **kwargs)
        dat.load_composition(**kwargs)

        # plot
        dat.plot_structure()
        dat.plot_composition()

        if clean:
            # move perplex files to their own folder
            new_dir = os.path.dirname(dat.data_path + 'output/' + dat.name + '/')
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)
            for fend in ['_comp.tab', '_adiabat.dat', '.dat']:
                file = dat.name + fend
                os.rename(dat.data_path + file, new_dir + '/' + file)


