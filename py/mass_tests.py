""" look at maximum interior water budget as a function of mass for fixed composition """
import numpy as np
import perplexdata as px
import matplotlib.pyplot as plt
from useful_and_bespoke import colorize, dark_background, cornertext
from parameters import M_E, M_SiO2, M_MgO, TO, G
import matplotlib
from matplotlib import rc
from matplotlib.pyplot import rcParams
import saturation as sat

# # masses = np.logspace(np.log10(0.1), np.log10(5), num=10)*M_E
# masses = np.logspace(np.log10(1), np.log10(5), num=4)*M_E
# m_w = []
# v_peak = []
# for ii, mass in enumerate(masses):
#     dat = px.build_planet(name='M' + str(ii), Tp=1900, M_p=mass, test_CMF=0.33, test_oxides=px.wt_oxides_Earth,
#                           star=None, plot=False, profile='warm', tol=0.1)
#     df = dat.df_all
#     w_max = np.sum(df['mass_h2o(kg)'])   # max water in kg
#     m_w.append(w_max)
#     g = G * dat.M_p / dat.R_p ** 2
#     hpeak = 100e6 / (0.5 * 2700 * g)
#     v_peak.append(4 * np.pi * dat.R_p ** 2 * hpeak)
#
#
# print('masses', [mass/M_E for mass in masses], 'M_E')
# print('max water in planet mass fraction', [w/m for w, m in zip(m_w, masses)])
# print('max water in TO', [w/TO for w in m_w])
# print('h peak vol in TO', [v/(TO / 1000) for v in v_peak])
#
# # for ii, mass in enumerate(masses):
# #     print('mass', mass/M_E, 'max wmf', m_w[ii]/mass, '(', ii + 1, '/', len(masses), ')')
#
#

def plot_melting(p=np.linspace(1, 140)):
    p = self.pressure_m * 1e-9  # GPa
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

    fig, ax = plt.subplots(1, 1)
    ax.plot(T_a, p, label='solid adiabat')
    ax.plot(T_melt, p, label='Peridotite solidus (dry)')
    # ax.plot(T_melt_Dorn, p, label='solidus Dorn eq (dry)')
    ax.invert_yaxis()
    ax.set_xlabel('T (K)')
    ax.set_ylabel('p (GPa)')
    ax.legend()
    fig, *axes = dark_background(fig, axes, )
    plt.show()
#