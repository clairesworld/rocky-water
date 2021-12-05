import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from useful_and_bespoke import colorize
import os

data_path = '/home/claire/Works/perple_x/'
fig_path = '/home/claire/Works/rccky-water/figs_scratch/'

M_O = 15.999
M_Si = 28.0855
M_Mg = 24.305
M_Ca = 40.078
M_Al = 26.981539
M_Fe = 55.845

# molar masses of oxides
M_SiO2 = M_Si + 2*M_O
M_CaO = M_Ca + M_O
M_MgO = M_Mg + M_O
M_Al2O3 = 2*M_Al + 3*M_O
M_FeO = M_Fe + M_O

class PerplexData:
    def __init__(self, name='default', oxides=['SIO2', 'MGO', 'CAO', 'AL2O3', 'FEO'], star='sun', **kwargs):
        self.name = name
        self.oxide_list = oxides
        self.star = star
        print('oxide list', oxides)

    def load_results(self, data_path=data_path, fig_path=fig_path, fillnan=True, **kwargs):
        file = data_path + self.name + '_1.tab'

        df = pd.read_csv(file, skiprows=8, index_col=None, sep=r"\s+",
                         # dtype=np.float64
                         )
        if fillnan:
            df = df.fillna(0)
        print(df.head())
        self.P = df['P(bar)'].to_numpy() * 1e5  # in Pa
        self.T = df['T(K)']
        self.df = df

        phases = df.columns.values.tolist()
        # remove df columns that aren't a mineral phase
        phases.remove('node#')
        phases.remove('P(bar)')
        phases.remove('T(K)')
        self.phases = phases

    # def get_adiabat(self, **kwargs):

    def get_hypatia(self, api_key='c53fd88b7e9c7e0c2e719ea4ea3c5e46',
                    ca_sol=6.33-12, al_sol=6.47-12, fe_sol=7.45-12, si_sol=7.52-12, mg_sol=7.54-12):
        """ star id is same as hypatia catalog with spaces e.g. "HIP 12345" """
        import requests

        els = []
        for ox in self.oxide_list:
            els.append(ox[:2].lower())

        if self.star != 'sun':  # not needed for solar values
            params = {"name": [self.star]*len(els), "element": els, "solarnorm": ["lod09"]*len(els)}
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

    def star_to_oxide(self, M_oxides=[M_SiO2, M_MgO, M_CaO, M_Al2O3, M_FeO],
                      core_efficiency=0.5):
        """ nH_star is list of the absolute log(N_X/N_H) for the star, in same order
         would need to convert from wrt solar if necessary
         by convention base element (denomenator for oxide wt ratios) is first in list
         core_efficiency is fraction of moles Fe that go into core instead of mantle FeO """

        wt_oxides = [1]  # give denomenator 1 for now
        for ii, el in enumerate(self.oxide_list):
            if ii > 0:
                if el == 'AL2O3':
                    X_ratio_mol = 0.5 * 10 ** self.nH_star[ii] / 10 ** self.nH_star[0]  # 2 mols Al per 1 mol Al2O3
                else:
                    X_ratio_mol = 10 ** self.nH_star[ii] / 10 ** self.nH_star[0]  # cancel out H abundance

                if el == 'FEO':
                    # some percentage of molar Fe goes to core instead of FeO
                    X_ratio_mol = X_ratio_mol * (1 - core_efficiency)

                # convert from mole ratio to mass ratio assuming all metals in oxides
                M_ratio = M_oxides[ii] / M_oxides[0]
                m_ratio = X_ratio_mol * M_ratio
                wt_oxides.append(m_ratio)

        # now have ratios in wt% 1:X1:X2:X3... normalise so total wt% is 100
        wt_oxides = np.array(wt_oxides)
        wt_oxides = wt_oxides / sum(wt_oxides) * 100
        print('wt % oxides', wt_oxides)
        self.wt_oxides = wt_oxides
        return wt_oxides

    def write_build(self, title='Planet', p_min=10000, p_max=245000, adiabat_file='aerotherm.dat', data_path=data_path,
                    overwrite=False):
        """ write perple_x build file, p in bar """
        build_file = data_path + self.name + '.dat'
        if os.path.isfile(build_file) and not overwrite:
            raise Exception('WARNING: build file', build_file, 'already exists, set overwrite=True')
        elif os.path.isfile(build_file):
            print('overwriting', build_file)

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
            s = s + '    0 unused place holder, post 06\n'*9
            s = s + '    0 number component transformations\n'
            s = s + '    5 number of components in the data base\n'
            s = s + '    1 component amounts, 0 - mole, 1 mass\n'
            s = s + '    0 unused place holder, post 06\n'*2
            s = s + '    0 unused place holder, post 05\n'
            s = s + '    5 ifug EoS for saturated phase\n'
            s = s + '    2 gridded minimization dimension (1 or 2)\n'
            s = s + '    0 special dependencies: 0 - P and T independent, 1 - P(T), 2 - T(P)\n'
            s = s + ' 0.00000      0.00000      0.00000      0.00000      0.00000     Geothermal gradient polynomial coeffs.\n\n'

            s = s + 'begin thermodynamic component list\n'
            for ii, el in enumerate(self.oxide_list):
                wt = self.wt_oxides[ii]
                dig = len(str(int(wt)))
                s = s + el.ljust(6) + '1'.ljust(3) + "{value:{width}.{precision}f}".format(value=float(wt), width=7, precision=6-dig)
                s = s + '      0.00000      0.00000     mass  amount\n'
            s = s + 'end thermodynamic component list\n\n\n'

            s = s + 'begin saturated component list\nend saturated component list\n\n\n'
            s = s + 'begin saturated phase component list\nend saturated phase component list\n\n\n'
            s = s + 'begin independent potential/fugacity/activity list\nend independent potential list\n\n\n'
            s = s + 'begin excluded phase list\nstv\nend excluded phase list\n\n\n'
            s = s + 'begin solution phase list\nWus(fab)\nPv(fab)\nO(stx)\nWad(stx)\nRing(stx)\nC2/c(stx)\nOpx(stx)\nCpx(stx)\nSp(stx)\nGt(stx)\nAki(fab)\nend solution phase list\n\n'

            s = s + "   {value:{width}.{precision}f}".format(value=float(p_max), width=7, precision=7-len(str(int(p_max))))[:-1] + '        0.00000        0.00000        0.00000        0.00000     max p, t, xco2, mu_1, mu_2\n'
            s = s + "   {value:{width}.{precision}f}".format(value=float(p_min), width=7,
                                                             precision=7 - len(str(int(p_min))))[
                    :-1] + '        0.00000        0.00000        0.00000        0.00000     min p, t, xco2, mu_1, mu_2\n'
            s = s + '   0.00000        0.00000        0.00000        0.00000        0.00000     unused place holder post 06\n\n'
            s = s + ' 1  2  4  5  3   indices of 1st & 2nd independent & sectioning variables'

            file.write(s)


    def get_composition(self, data_path=data_path):
        # import subprocess
        os.chdir(data_path)

        # create vertex command file
        with open(data_path + self.name + '_vertex_command.txt', 'w') as file:
            s = self.name + '\n0'
            file.write(s)

        # run vertex
        os.system('./vertex < ' + self.name + '_vertex_command.txt')

        # create werami command file
        with open(data_path + self.name + '_werami_command.txt', 'w') as file:
            s = self.name + '\n'  # Enter the project name (the name assigned in BUILD)
            s = s + '3\n'  # Select operational mode: 3 - properties along a 1d path
            s = s + '25\n'  # Select a property: 25 - Modes of all phases
            s = s + 'n\n'  # Output cumulative modes (y/n)?
            s = s + '0\n'  # Select operational mode: 0 - EXIT
            file.write(s)

        # run werami
        os.system('./werami < ' + self.name + '_werami_command.txt')

    def plot_composition(self, which='pressure', **kwargs):

        fig, ax = plt.subplots(1, 1)
        if which == 'pressure':
            x = self.P * 1e-9  # plot P in GPa
            ax.set_xlabel('Pressure (GPa)')
        else:
            raise NotImplementedError('independent variables other than pressure not implemented')

        ax.set_ylabel('Modal abundance (mol %)')
        colors = colorize(range(len(self.phases)), cmap='rainbow')[0]
        for ii, phase in enumerate(self.phases):
            y = self.df[phase]
            plt.plot(x, y, c=colors[ii], label=phase)
        ax.legend()
        plt.title(self.name)
        plt.show()


stars = ["HIP 98355", "sun"]
for star in stars:

    dat = PerplexData(name=star.replace(" ", ""), star=star)
    dat.get_hypatia()
    dat.star_to_oxide()
    dat.write_build(overwrite=True)
    dat.get_composition()
    dat.load_results()
    dat.plot_composition()
