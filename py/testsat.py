import numpy as np
import perplexdata as run
import parameters as p
import eos
import fugacity as fug
import saturation as sat
import matplotlib.pyplot as plt
from useful_and_bespoke import colorize


def test_geotherm(adiabat_file='geotherm_Tp1600.csv', path='', res=100):
    import pandas as pd
    T_geo, p_geo = sat.get_geotherm_file(adiabat_file, path=path, res=res)
    df = pd.DataFrame({'p': p_geo, 'T': T_geo})

    # test with modal water contents only
    df['X_ol'] = 0
    df['X_wd'] = 0
    df['X_rw'] = 0
    df['X_bdg'] = 0
    df.loc[df['p'] < 14.2e9, 'X_ol'] = 1
    df.loc[(df['p'] >= 14.2e9) & (df['p'] < 18e9), 'X_wd'] = 1
    df.loc[(df['p'] >= 18e9) & (df['p'] < 23e9), 'X_rw'] = 1
    df.loc[(df['p'] >= 23e9), 'X_bdg'] = 1
    df = sat.mineral_water_contents(p_geo, T_geo, df=df)
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    print(df)
    return df


""" test ppv partitioning """
T_geo, p_geo = sat.get_geotherm_file('geotherm_Tp1600.csv', res=100)#np.array([2264, 2437]), np.array([90, 135])*1e9
sat.ppv_partitioning(p_geo=p_geo, T_geo=T_geo, file='D_ppv_brg.csv', path='', res=100, plot=True)


# test_geotherm()

""" check some partitioning coeffs """

""" like figure A1 of Otsuka & Karato 2011 """
# # p, T = 5e9, 1873
# p = np.linspace(0.1, 15)*1e9
# T = np.array([1273, 1473, 1673, 1873])
# plt.figure()
# for TT in T:
#     f_h2o = fugacity(p=p, T=TT)  # Pa
#     plt.gca().plot(p*1e-9, f_h2o*1e-9, label=str(TT) + ' K')
# plt.xlabel('Pressure (GPa)')
# plt.ylabel('Fugacity (GPa)')
# plt.legend()
# plt.gca().set_yscale('log')
# plt.gca().set_yticks([1e-5, 1e0, 1e5, 1e10])
# plt.gca().set_xticks([0, 5, 10, 15])
# plt.gca().set_xlim([0, 15])


""" like figure 1 of Dong+ 2021 """
# fig = plt.figure()
# plt.title('Olivine')
# ax = plt.gca()
# T = np.linspace(1200, 2400)
# p = np.linspace(0.1, 15.5)*1e9
# TT, pp = np.meshgrid(T, p)
# f_h2o = sat.fugacity(p=pp, T=TT)  # Pa
# w_ol = sat.water_capacity(p=pp, T=TT, f_h2o=f_h2o, X_Fe=0, n=0.6447, a=0, b=4905.5403, c=0, d=0)
# CS = ax.contour(TT, pp*1e-9, np.log10(w_ol), cmap='coolwarm_r')
# plt.xlabel('Temperature (K)')
# plt.ylabel('Pressure (GPa)')
# CB = fig.colorbar(CS, shrink=0.8)
# CB.ax.set_ylabel('Water storage capacity (log10(ppm wt.))')
# fig, axes = plt.subplots(2, 1)
# T = np.linspace(1500, 2300)
# w_wd = water_capacity(p, T, f_h2o=0, X_Fe=0, n=0, a=-7.6356, b=13739.5371, c=0, d=0)   # in wt %
# w_rw = water_capacity(p, T, f_h2o=0, X_Fe=0, n=0, a=-6.8856, b=12206.2676, c=0, d=0)
# axes[0].plot(T, w_wd)
# axes[1].plot(T, w_rw)
# for ax in axes:
#     ax.set_ylabel('Water storage capacity (wt. %)')
#     ax.set_ylim([0, 4])
# axes[1].set_xlabel('Temperature (K)')
# axes[0].set_title('Wadsleyite')
# axes[1].set_title('Ringwoodite')
# plt.tight_layout()
#
# plt.show()


""" enstatite from Mierdel and Keppler (2004) - reproduce fig 6
also fig 4 keppler & b-f """
# # T, p = get_geotherm_file()
# # deltaH = -4.563  # kJ/mol
# # deltaV = 12.1  # cm3/mol
# # n = 1
# # A = 0.01354  # ppm/bar
# # p_bar = np.array([10e3, 20e3, 40e3, 60e3, 80e3, 100e3])
# # aluminious
# # deltaH = -79.685  # kJ/mol
# # deltaV = 11.3  # cm3/mol
# # n = 0.5
# # A = 0.042  # ppm/bar
# T = np.linspace(700, 1600) + 273
# fig, axes = plt.subplots(2, 1)
# p_bar = np.array([10e3, 20e3, 40e3])
# col = colorize(p_bar)[0]
# for p, c in zip(p_bar, col):
#     f_h2o = sat.fugacity(p=p*1e5, T=T, func=fug.PSfugacity)  # Pa
#     c_w_al = sat.water_capacity_thermo(p*1e5, T, f_h2o, A=0.042, n=0.5, deltaH=-79.685, deltaV=11.3)*1e6  # ppm
#     c_w_en = sat.water_capacity_thermo(p * 1e5, T, f_h2o, A=0.01354, n=1, deltaH=-4.563, deltaV=12.1) * 1e6  # ppm
#     # w_ol = sat.water_capacity(p=p*1e5, T=T, f_h2o=f_h2o, X_Fe=0, n=0.6447, a=0, b=4905.5403, c=0, d=0)
#     axes[0].plot(T - 273, c_w_al + c_w_en, c=c, label=str(p*1e-3)+' kbar')
#     # plt.plot(T - 273, w_ol, c=c, label=str(p*1e-3) + ' kbar, ol', ls='--')
# axes[0].set_ylim((0, 25000))
# axes[0].legend()
# p_bar = np.array([1, 60e3, 80e3, 100e3])
# col = colorize(p_bar)[0]
# for p, c in zip(p_bar, col):
#     f_h2o = sat.fugacity(p=p * 1e5, T=T, func=fug.PSfugacity)  # Pa
#     c_w_al = sat.water_capacity_thermo(p * 1e5, T, f_h2o, A=0.042, n=0.5, deltaH=-79.685, deltaV=11.3) * 1e6  # ppm
#     c_w_en = sat.water_capacity_thermo(p * 1e5, T, f_h2o, A=0.01354, n=1, deltaH=-4.563, deltaV=12.1) * 1e6  # ppm
#     # w_ol = sat.water_capacity(p=p * 1e5, T=T, f_h2o=f_h2o, X_Fe=0, n=0.6447, a=0, b=4905.5403, c=0, d=0)
#     axes[1].plot(T - 273, c_w_al + c_w_en, c=c, label=str(p * 1e-3) + ' kbar')
# axes[1].set_ylim((0, 3000))
# axes[1].legend()
# plt.xlabel('T (C)')
# plt.ylabel('water solubilty (ppm)')
# axes[0].set_title('opx')
# # plt.ylim((100, 1700))


""" test aluminious opx-ol partitioning from bolfan-cassanova and keppler review figure 11 """
# # deltaH = -79.685  # kJ/mol
# # deltaV = 11.3  # cm3/mol
# # n = 0.5
# # A = 0.042  # ppm/bar
# p_bar = np.linspace(20, 130)*1e3
# T = np.array([600, 800, 1000, 1200, 1400]) + 273
# plt.figure()
# for TT in T:
#     f_h2o = sat.fugacity(p=p_bar*1e5, T=TT, func=fug.PSfugacity)  # Pa
#     c_w_al = sat.water_capacity_thermo(p_bar * 1e5, TT, f_h2o, A=0.042, n=0.5, deltaH=-79.685, deltaV=11.3) * 1e6  # ppm
#     c_w_en = sat.water_capacity_thermo(p_bar * 1e5, TT, f_h2o, A=0.01354, n=1, deltaH=-4.563, deltaV=12.1) * 1e6  # ppm
#     c_opx = c_w_al + c_w_en
#     c_ol = sat.water_capacity_thermo(p_bar * 1e5, TT, f_h2o, A=0.0066, n=1, deltaH=0, deltaV=10.6) * 1e6
#     # c_ol_D21 = sat.water_capacity(p_bar*1e5, TT, f_h2o, X_Fe=0, n=0.6447, a=0, b=4905.5403, c=0, d=0)*1e-6
#     D_w = c_ol/c_opx
#     plt.plot(p_bar*1e-3, D_w, label=str(TT - 273)+' C')
#     # plt.plot(p_bar * 1e-3, c_ol, label=str(TT - 273) + ' C')
#     # plt.plot(p_bar * 1e-3, c_ol_D21, label=str(TT - 273) + ' C', ls='--')
# plt.xlabel('p (kbar)')
# plt.ylabel('$P_{ol/opx}$')
# # plt.ylabel('water content (ppm)')
# # plt.ylim((0, 1))
# plt.legend()


""" compare olivine parameterisations """
# T = np.linspace(700, 1600) + 273
# fig, axes = plt.subplots(1, 1)
# p_bar = np.array([10e3, 20e3, 40e3, 60e3, 80e3])
# col = colorize(p_bar)[0]
# for p, c in zip(p_bar, col):
#     f_h2o_PS = sat.fugacity(p=p * 1e5, T=T, func=fug.PSfugacity)  # Pa
#     f_h2o_FW = sat.fugacity(p=p * 1e5, T=T)  # Pa
#     w_ol_K = sat.water_capacity_thermo(p * 1e5, T, f_h2o_PS, A=0.0066, n=1, deltaH=0, deltaV=10.6) * 1e6
#     w_ol_D = sat.water_capacity(p=p * 1e5, T=T, f_h2o=f_h2o_FW, X_Fe=0, n=0.6447, a=0, b=4905.5403, c=0, d=0)  # ppm
#     plt.plot(T - 273, w_ol_D, c=c, label=str(p * 1e-3) + ' kbar, Dong', ls='-')
#     plt.plot(T - 273, w_ol_K, c=c, label=str(p * 1e-3) + ' kbar, Kohlstedt', ls='--')
# # axes[0].set_ylim((0, 25000))
# axes.legend()
# plt.xlabel('T (C)')
# plt.ylabel('water content (ppm)')

plt.show()
