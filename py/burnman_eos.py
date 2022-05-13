import os
import sys
import numpy as np

# sys.path.insert(1, '/home/claire/Works/burnman/')
# os.chdir('/home/claire/Works/burnman/')

import burnman
# from burnman.tools.chemistry import dictionarize_formula, formula_mass
from burnman import minerals

# define mineral properties from stixrude 2021 database


# update with stx21 dataset
params_ppv = {'name': 'Mg_Post_Perovskite', 'formula': {'Mg': 1.0, 'Si': 1.0, 'O': 3.0},
               'equation_of_state': 'slb3',
               'F_0': -1313625.76,  # -1348641.0,  Helmoltz energy, J/mol
               'V_0': 2.3525e-05,  # 2.4419e-05,  Volume at 300 K and 1e5 Pa, m3/mol
               'K_0': 2920000e5,  # 231200000000.0,  Isothermal bulk modulus, Pa
               'Kprime_0': 3.740,  # 4.0,  Pressure derivative of K_0
               'Debye_0': 941.49795,  # 855.8173,  theta0 - stx gives 941 but Sakai gives 971
               'grueneisen_0': 1.76958,  # 1.89155,  gamma0 - stx21 gives 1.77 but Sakai gives 1.31
               'q_0': 2.04631,  # 1.09081,  Mie-Gruneisen exponent - stx gives 2.04 but Sakai gives 2.5
               'G_0': 1711358.7e5,  # 150167000000.0,  Adiabatic shear modulus, Pa
               'Gprime_0': 1.85188,  # 1.97874,  Pressure derivative of G_0
               'eta_s_0': 1.28810,  # 1.16704,  Shear strain derivative of tensorian Gruneisen parameter
               'n': 5.0,  # number of atoms per formula unit
               'molar_mass': 0.1003887,
               'T_0': 300.0,
               'E_0': 0.0, 'P_0': 0.0}  # aka mppv

params_fppv = {'name': 'Fe_Post_Perovskite', 'formula': {'Fe': 1.0, 'Si': 1.0, 'O': 3.0},
               'equation_of_state': 'slb3',
               'F_0': -982035.50,  # -981806.9,
               'V_0': 2.4652e-5,  # 2.5459e-05,
               'K_0': 2920000e5,  # 231200000000.0,
               'Kprime_0': 3.740,  # 4.0,
               'Debye_0': 794.15823,  # 781.3465,
               'grueneisen_0': 1.76958,  # 1.89155,
               'q_0': 2.04631,  # 1.09081,
               'G_0': 1295000e5,  # 129500000000.0,
               'Gprime_0': 1.31526,  # 1.44675,
               'eta_s_0': 1.72601,  # 1.36382,
               'n': 5.0,
               'molar_mass': 0.1319287,
               'T_0': 300.0, 'E_0': 0.0, 'P_0': 0.0}

params_appv = {'name': 'Al_Post_Perovskite', 'formula': {'Al': 2.0, 'O': 3.0}, 'equation_of_state': 'slb3',
               'F_0': -1336464.73,  # -1377582.0,
               'V_0': 2.3847e-5,  # 2.3847e-05,
               'K_0': 2490000e5,  # 249000000000.0,
               'Kprime_0': 4.0,
               'Debye_0': 722.93835,  # 762.1951,
               'grueneisen_0': 1.88758,  # 1.64573,
               'q_0': 2.04631,  # 1.09081,
               'G_0': 919652.6e5,  # 91965310000.0,
               'Gprime_0': 1.81603,  # 1.81603,
               'eta_s_0': 2.52605,  # 2.83762,
               'n': 5.0, 'molar_mass': 0.1019612, 'T_0': 300.0, 'E_0': 0.0, 'P_0': 0.0}

params_capv = {'name': 'Ca_Perovskite', 'formula': {'Ca': 1.0, 'Si': 1.0, 'O': 3.0}, 'equation_of_state': 'slb3',
               'F_0': -1459910.18,  # -1463358.0,
               'V_0': 2.7450e-5,  # 2.745e-05,
               'K_0': 236000000000.0,
               'Kprime_0': 3.9,
               'Debye_0': 798.78581,  # 795.779,
               'grueneisen_0': 1.88943,  # 1.88839,
               'q_0': 0.89662,  # 0.89769,
               'G_0': 1552052.4e5,  # 156831500000.0,
               'Gprime_0': 2.22637,  # 2.22713,
               'eta_s_0': 1.23493,  # 1.28818,
               'n': 5.0, 'molar_mass': 0.1161617, 'T_0': 300.0, 'E_0': 0.0, 'P_0': 0.0}

params_fapv = {'name': 'FeAlO3-Perovskite', 'formula': {'Fe': 1.0, 'Al': 1.0, 'O': 3.0}, 'equation_of_state': 'slb3',
               'F_0': -848315.02,
               'V_0': 2.7260e-5,
               'K_0': 2233255e5,
               'Kprime_0': 4.1,
               'Debye_0': 755.38117,
               'grueneisen_0': 1.54222,
               'q_0': 0.84088,
               'G_0': 1500420.9e5,
               'Gprime_0': 1.73259,
               'eta_s_0': 2.55505,
               'n': 5.0, 'molar_mass': 0.13082353900000002, 'T_0': 300.0, 'E_0': 0.0, 'P_0': 0.0}

params_ffer = {'name': 'Fe-Ca-Ferrite', 'formula': {'Fe': 1.0, 'Al': 2.0, 'O': 4.0}, 'equation_of_state': 'slb3',
               'F_0': -1774197.05,
               'V_0': 3.7216e-5,
               'K_0': 2130000e5,
               'Kprime_0': 4.1,
               'Debye_0': 734.07527,
               'grueneisen_0': 1.56672,
               'q_0': 1,
               'G_0': 1597096.5e5,
               'Gprime_0': 1.93591,
               'eta_s_0': 2.34163,
               'n': 7.0, 'molar_mass': 173.804078e-3, 'T_0': 300.0, 'E_0': 0.0, 'P_0': 0.0}
# note - in burnman this is fecf

params_per = {'name': 'Periclase', 'formula': {'Mg': 4.0, 'O': 4.0}, 'equation_of_state': 'slb3',
              'F_0': -2278109.88,  # -569444.6,
              'V_0': 4.4976e-5,  # 1.1244e-05,
              'K_0': 1611439.3e5,  # 161383600000.0,
              'Kprime_0': 3.90838,  # 3.84045,
              'Debye_0': 770.90151,  # 767.0977,
              'grueneisen_0': 1.45033,  # 1.36127,
              'q_0': 1.54870,  # 1.7217,
              'G_0': 130900000000.0,
              'Gprime_0': 2.14668,  # 2.1438,
              'eta_s_0': 2.56123,  # 2.81765,
              'n': 8.0,  # 2.0
              'molar_mass': 0.16121760000000002,  # 0.040304400000000004
              'T_0': 300.0, 'E_0': 0.0, 'P_0': 0.0}
# params_per11 = minerals.SLB_2011.periclase().params

params_wus = {'name': 'Wuestite', 'formula': {'Fe': 4.0, 'O': 4.0}, 'equation_of_state': 'slb3',
              'F_0': -974607.49,  # -242146.0,
              'V_0': 4.9024e-05,  # 1.2264e-05,
              'K_0': 1607000e5,  # 179444200000.0,
              'Kprime_0': 4,  # 4.9376,
              'Debye_0': 454.17520,  # 454.1592,
              'grueneisen_0': 1.45033,  # 1.53047,
              'q_0': 1.54870,  # 1.7217,
              'G_0': 59000000000.0,
              'Gprime_0': 1.44764,  # 1.44673,
              'eta_s_0': 0.06776,  # -0.05731,
              'n': 8.0,  # 2.0,
              'molar_mass': 0.2873776,  # 0.0718444,
              'T_0': 300.0, 'E_0': 0.0, 'P_0': 0.0}
# params_wus11 = minerals.SLB_2011.wuestite().params

params_wuls = {'name': 'Wuestite_low_spin', 'formula': {'Fe': 4.0, 'O': 4.0}, 'equation_of_state': 'slb3',
              'F_0': -621968.16,
              'V_0': 4.3400e-05,
              'K_0': 1997000e5,
              'Kprime_0': 4,
              'Debye_0': 524.57881,
              'grueneisen_0': 1.45033,
              'q_0': 1.54870,
              'G_0': 59000000000.0,
              'Gprime_0': 1.44073,
              'eta_s_0': -0.13801,
              'n': 8.0,  # 2.0,
              'molar_mass': 0.2873776,  # 0.0718444,
              'T_0': 300.0, 'E_0': 0.0, 'P_0': 0.0}

params_seif = {'name': 'Seifertite', 'formula': {'Si': 1.0, 'O': 2.0}, 'equation_of_state': 'slb3',
               'F_0': -793366.84,  # -794335.4,
               'V_0': 1.367e-05,
               'K_0': 3271560.1e5,  # 327584300000.0,
               'Kprime_0': 4.01662,  # 4.01553,
               'Debye_0': 1128.94590,  # 1140.772,
               'grueneisen_0': 1.55674,  # 1.37466,
               'q_0': 2.20960,  # 2.83517,
               'G_0': 2274115.9e5,  # 227453200000.0,
               'Gprime_0': 1.77078,  # 1.76965,
               'eta_s_0': 4.55828,  # 4.97108,
               'n': 3.0, 'molar_mass': 0.0600843, 'T_0': 300.0, 'E_0': 0.0, 'P_0': 0.0}


# define mineral classes for lower mantle
# mppv = burnman.Mineral(params=params_mppv)
# fppv = burnman.Mineral(params=params_fppv)
# appv = burnman.Mineral(params=params_appv)
# capv = burnman.Mineral(params=params_capv)
# per = burnman.Mineral(params=params_per)
# wus = burnman.Mineral(params=params_wus)
# seif = burnman.Mineral(params=params_seif)

def get_properties(P, T, phase_name):  # Pa, K, str
    # phase name must match one of the datasets above
    # P, T can be arrays and must be same size
    if np.size(T) != np.size(P):
        raise Exception('T and P must match size')
    if np.size(T) == 1:
        T = np.array(T)
    if np.size(P) == 1:
        P = np.array(P)

    properties = ['molar_volume', 'density', 'thermal_expansivity', 'molar_heat_capacity_p', 'molar_mass']
    try:
        params = eval('params_' + phase_name)
    except NameError:  # if not given above, hopefully in BurnMan's SLB database?
        try:
            params = eval('minerals.SLB_2011.' + phase_name + '().params')
        except NameError:
            print('\n>>>>WARNING>>>> no entry for', phase_name, 'in SLB_2011')  # return 0?
            # return [0, 0, 0, 0]

    phase = burnman.Mineral(params=params)
    V, rho, alpha, cp_mol, M = phase.evaluate(properties, P, T)
    cp = cp_mol / M
    return [V, rho, alpha, cp]

# import matplotlib.pyplot as plt
# P = np.linspace(200e9, 400e9)
# T = 2000*np.ones_like(P)
# fig, axes = plt.subplots(1, 3)
# phases = ['wus', 'wus11', 'per', 'per11']
# ls = ['-', '--', '-', '--']
# for ii, phase in enumerate(phases):
#     _, rho, alpha, cp = get_properties(P, T, phase)
#     axes[0].plot(P, rho, label=phases[ii], alpha=0.5, ls=ls[ii])
#     axes[1].plot(P, alpha, alpha=0.5, ls=ls[ii])
#     axes[2].plot(P, cp, alpha=0.5, ls=ls[ii])
# axes[0].set_ylabel('density')
# axes[1].set_ylabel('alpha')
# axes[2].set_ylabel('cp')
# axes[0].legend()
# plt.show()