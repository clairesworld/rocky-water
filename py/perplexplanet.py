import perplexdata as px
import numpy as np
# import matplotlib.pyplot as plt
import pandas as pd
from parameters import M_E, M_Fe, M_FeO, M_MgO, M_SiO2, M_Si, M_Mg, M_Ca, M_CaO, M_Al, M_Al2O3, G, R_E, rho_E
import os
import pathlib
import subprocess
from bulk_composition import stellar_mantle
import ask_hypatia as hyp

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

class PerplexPlanet(px.PerplexData):

    def __init__(self, name='test', core_efficiency=0.88, M_p=M_E, R_p=None,
                 oxides=None, solution_phases=None,
                 star='sun', perplex_path=perplex_path_default, output_parent_path=output_parent_default, verbose=False,
                 **kwargs):

        super().__init__(**kwargs)
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