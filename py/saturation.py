import numpy as np
import fugacity as fug
from scipy.interpolate import interp1d
from useful_and_bespoke import find_nearest_idx, colorize
from scipy import optimize as op
import matplotlib.pyplot as plt
import pandas as pd

Rb = 8.31446261815324


def water_capacity(p, T, f_h2o, X_Fe=0, n=0, a=0, b=0, c=0, d=0):
    """
    water storage capcity of NAMs in H2O-saturated conditions
    p in Pa, T in K, f_h2o (water fugacity) in Pa, X_Fe (iron content) in mole fraction
    """
    if n == 0:
        return np.exp(a + (b + c * p * 1e-9 + d * X_Fe) / T)
    # return np.exp(np.log(1 / (2 * n)) + a + n / 2 * np.log(f_h2o * 1e9) + (b + c * p * 1e9 + d * X_Fe) / T)
    return np.exp(a + n / 2 * np.log(f_h2o * 1e-9) + (b + c * p * 1e-9 + d * X_Fe) / T)

def water_capacity_thermo(p, T, f_h2o, A, n, deltaH, deltaV):
    """ delta H in kJ/mol, deltaV in cm3/mol, A in ppm/bar^n, p in Pa, T in K, fh2o in Pa"""
    f_h2o = f_h2o * 1e-5   # Pa to bar
    deltaV = deltaV / (100 ** 3)
    deltaH = deltaH * 1e3
    c = A * f_h2o ** n * np.exp(-(deltaH + deltaV * p) / (Rb * T))
    return c*1e-6  # weight fraction

def partitioning_coeff(p, T, f_h2o, A_j, deltaH_j, deltaV_j, n_j=0.5, A_i=0.0066, deltaH_i=0, deltaV_i=10.6, n_i=1):
    """
    A in ppm/bar, delta H in kJ/mol, delta V in cm3/mol, f, p and T in Pa and K
    returns weight fraction
    """
    f_h2o = f_h2o*1e-5  # to bar for use with A
    deltaV_j = deltaV_j / (100 ** 3)  # from cm3 to m3
    deltaV_i = deltaV_i / (100 ** 3)  # from cm3 to m3
    deltaH_i = deltaH_i * 1e3  # from kJ to J
    deltaH_j = deltaH_j * 1e3  # from kJ to J
    D_ij = A_i / A_j * f_h2o ** (n_i - n_j) * np.exp(
        -((deltaH_i - deltaH_j) + (deltaV_i - deltaV_j) * p) / (Rb * T))

    return D_ij


def fugacity(p, T, func=fug.fugacity_frost, **kwargs):
    if np.size(T) > 1:  # iterate over T
        if np.size(p) == 1:
            p = [p] * len(T)

        f = []
        for pp, TT in zip(p, T):
            f.append(func(pp, TT, **kwargs))
        return np.array(f)
    elif np.size(p) > 1:   # iterate over p
        if np.size(T) == 1:
            T = [T] * len(p)
        f = []
        for pp, TT in zip(p, T):
            f.append(func(pp, TT, **kwargs))
        return np.array(f)
    else:
        return func(p, T, **kwargs)


def ppv_partitioning(p_geo=None, T_geo=None, file='D_ppv_brg.csv', path='', res=10000):
    import pandas as pd

    # find where phase boundary intersects geotherm
    # if T_geo is None or p_geo is None:
    #     df1 = pd.read_csv(path + 'geotherm_Tp1600', names=['P_GPa', 'T_K'], index_col=False, skiprows=6)
    #     p_geo = df1['P_GPa'].to_numpy()
    #     T_geo = df1['T_K']
    TP_geo = interp1d(p_geo, T_geo, kind='linear')
    f = TP_geo(np.linspace(np.min(p_geo), np.max(p_geo), num=res))  # geotherm function

    df2 = pd.read_csv(path + 'ppv_phase_boundary.csv', names=['P_GPa', 'T_K'], index_col=False, skiprows=6)
    TP = interp1d(df2['P_GPa'].to_numpy(), df2['T_K'].to_numpy(), kind='linear')
    g = TP(np.linspace(np.min(df2['P_GPa'].to_numpy()), np.max(df2['P_GPa'].to_numpy()), num=res))  # phase bdy fn

    idx = np.argwhere(np.diff(np.sign(f - g))).flatten()
    T_trans = f[idx]

    # use this temperature to get partitioning coefficient
    df3 = pd.read_csv(path + file, names=['T_K', 'D'], index_col=False, skiprows=6)
    D = interp1d(df3['T_K'].to_numpy(), df3['D'].to_numpy(), kind='linear')
    D_trans = D(T_trans)
    print('T, D at phase boundary', T_trans, D_trans)
    return D_trans


def mineral_water_contents(p, T, X_Fe=0, df=None):  # p in Pa, T in K, X_Fe has no effect atch
    f_h2o = fugacity(p=p, T=T)  # Pa
    print()

    """
    OLIVINE
    The fitted value of n may represent the average number of hydroxyls through multiple substitution mecha-
    nisms over the wide range of pressure and water concentration for olivine (Ferot & Bolfan-Casanova, 2012;
    Otsuka & Karato, 2011). Temperature dependence of the water storage capacity in olivine (b/T) can be
    attributed to the decrease in the activity of the H 2 O component in the aqueous fluid or silicate melt with
    increasing temperature (Mibe et al., 2002; Stalder et al., 2001). 
    """
    w_ol = water_capacity(p, T, f_h2o, X_Fe, n=0.6447, a=0, b=4905.5403, c=0, d=0)*1e-6  # in ppm wt?

    """
    WADSLEYITE & RINGWOODITE Mg2SiO4
    """
    w_wd = water_capacity(p, T, f_h2o, X_Fe, n=0, a=-7.6356, b=13739.5371, c=0, d=0)*1e-2  # in wt %?
    w_rw = water_capacity(p, T, f_h2o, X_Fe, n=0, a=-6.8856, b=12206.2676, c=0, d=0)*1e-2

    """
    ORTHOPYROXENE, CLINOPYROXENE, GARNET
    from Ol partitioning
    """
    # Al-free enstatite MgSiO3
    D_ol_en = partitioning_coeff(p, T, f_h2o, A_j=0.01354, deltaH_j=-4.563, deltaV_j=12.1, n_j=1)  # their n is a typo based on Table
    w_en = w_ol / D_ol_en

    # Al-bearing opx
    D_ol_Alopx = partitioning_coeff(p, T, f_h2o, A_j=0.042, deltaH_j=-79.685, deltaV_j=11.3, n_j=0.5)
    w_Alopx = w_ol / D_ol_Alopx

    # cpx
    D_ol_cpx = partitioning_coeff(p, T, f_h2o, A_j=7.144, deltaH_j=0, deltaV_j=8.019, n_j=0.5)
    w_cpx = w_ol / D_ol_cpx

    # gt
    D_ol_gt = partitioning_coeff(p, T, f_h2o, A_j=0.679, deltaH_j=0, deltaV_j=5.71, n_j=0.5)
    w_gt = w_ol / D_ol_gt

    # print('D_ol_en', D_ol_en)
    # print('D_ol_alopx', D_ol_Alopx)
    # print('D_ol_cpx', D_ol_cpx)
    # print('D_ol_gt', D_ol_gt)

    # TODO: what do we do if we have no olivine? can you write a "theoretical" olivine water content?

    """
    SILICATE PEROVSKITES
    Calcium silicate perovskite may have a higher water storage capacity (Chen et al., 2020), especially at the P-T conditions
    of the uppermost lower mantle (Muir & Brodholt, 2018). However, this high-pressure mineral converts to an
    amorphous phase upon pressure release, so its water storage capacity has not been quantitatively constrained
    (Chen et al., 2020; Nemeth et al., 2017). The water storage capacity in ferropericlase is likely to be relatively low
    throughout the lower mantle (Bolfan-Casanova et al., 2002; Litasov, 2010).
    """
    # Bridgmanite MgSiO3 >>> this is half of Earth's mantle mass so maybe the most important?
    D_rw_bdg = 15  # +- 8, Inoue+ 2010, at 1873 K
    w_bdg = w_rw / D_rw_bdg

    # Ca-perovskite
    w_capv = 10e-6  # assumed in Dong+ 2021

    # ferropericlase
    w_fp = 10e-6  # assumed in Dong+ 2021

    """ 
    POSTPEROVSKITE
    partitioning with bridgmanite calculated in Townsend+ 2016
    water partitioning controlled by presence of Al because different substitution mechanisms
    >>> do you need to calculate intersection of partitioning T-dependence curve with phase boundary and use that value?
    """
    # Al-free postperovskite - same Mg substitution favours Bdg
    D_ppv_brg = ppv_partitioning(p_geo=p, T_geo=T, file='D_ppv_brg.csv', path='', res=10000)
    w_ppv = w_bdg * D_ppv_brg

    # Al-bearing ppv - Si-AlH defect favours Ppv
    D_Alppv_brg = ppv_partitioning(p_geo=p, T_geo=T, file='D_Alppv_brg.csv', path='', res=10000)
    w_Alppv = w_bdg * D_Alppv_brg

    """
    other phases by partitioning
    """
    # high pressure clinoenstatite
    D_ol_hpcpx = 0.8  # Withers+ 2007
    w_hpcpx = w_ol / D_ol_hpcpx

    # akimotoite?? - upper mantle, might not intersect

    # spinel?? MgAl2O4 - upper mantle, might want to find D

    # fill in blanks in partitioning - D_spinel/melt / D_olivine/melt

    # C2/c??

    # stishovite??

    # magnesiowustite??

    if df is not None:
        # compile
        df['f_h2o'] = f_h2o  # Pa
        df['c_h2o'] = 0
        filter_col = [col for col in df if col.startswith('X_')]  # get all columns corresponding to mole fractions
        for col in filter_col:
            phase = col[2:]
            df['w_' + phase] = eval('w_' + phase)
            df['c_h2o'] = df['c_h2o'] + df[col] * df['w_' + phase]
        return df


def get_geotherm_file(file='geotherm_Tp1600.csv', path='', res=100):
    df_read = pd.read_csv(path + file, names=['T_K', 'P_GPa'], index_col=False, skiprows=6)
    TP = interp1d(df_read['P_GPa'].to_numpy(), df_read['T_K'].to_numpy(), kind='linear')
    p_geo = np.linspace(np.min(df_read['P_GPa'].to_numpy()), np.max(df_read['P_GPa'].to_numpy()), num=res)
    T_geo = TP(p_geo)  # interpolated T as a function of p in GPa
    p_geo = p_geo * 1e9  # convert to Pa
    return T_geo, p_geo

# test
def test_geotherm(adiabat_file='geotherm_Tp1600.csv', path='', res=100):
    T_geo, p_geo = get_geotherm_file(adiabat_file, path=path, res=res)
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

    df = mineral_water_contents(p_geo, T_geo, df=df)

    pd.set_option("display.max_rows", None, "display.max_columns", None)
    print(df)
    return df


# test_geotherm()



# # like figure A1 of Otsuka & Karato 2011
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


# # like figure 1 of Dong+ 2021
#
fig = plt.figure()
plt.title('Olivine')
ax = plt.gca()
T = np.linspace(1200, 2400)
p = np.linspace(0.1, 15.5)*1e9
TT, pp = np.meshgrid(T, p)
f_h2o = fugacity(p=pp, T=TT)  # Pa
print(np.shape(f_h2o))
w_ol = water_capacity(p=pp, T=TT, f_h2o=f_h2o, X_Fe=0, n=0.6447, a=0, b=4905.5403, c=0, d=0)
print(np.shape(w_ol))
CS = ax.contour(TT, pp*1e-9, np.log10(w_ol), cmap='coolwarm_r')
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (GPa)')
CB = fig.colorbar(CS, shrink=0.8)
CB.ax.set_ylabel('Water storage capacity (log10(ppm wt.))')
#
#
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


# p = 4e9
# T = 1673
# f_h2o = fugacity(p=p, T=T)  # Pa
# D_ol_en = partitioning_coeff(p, T, f_h2o, A_j=0.01354, deltaH_j=-4.563, deltaV_j=12.1, n_j=1)
# print('ol-opx', D_ol_en)



# test enstatite from Mierdel and Keppler (2004) - reproduce fig 6 - approx the same
# T, p = get_geotherm_file()
# deltaH = -4.563  # kJ/mol
# deltaV = 12.1  # cm3/mol
# n = 1
# A = 0.01354  # ppm/bar
# aluminious
deltaH = -79.685  # kJ/mol
deltaV = 11.3  # cm3/mol
n = 0.5
A = 0.042  # ppm/bar
p_bar = np.array([1, 10e3, 20e3, 40e3, 60e3, 80e3, 100e3])
T = np.linspace(700, 1600) + 273
plt.figure()
col = colorize(p_bar)[0]
for p, c in zip(p_bar, col):
    f_h2o = fugacity(p=p*1e5, T=T)  # Pa
    c_w = water_capacity_thermo(p*1e5, T, f_h2o, A, n, deltaH, deltaV)*1e6  # ppm
    w_ol = water_capacity(p=p*1e5, T=T, f_h2o=f_h2o, X_Fe=0, n=0.6447, a=0, b=4905.5403, c=0, d=0)
    plt.plot(T - 273, c_w, c=c, label=str(p*1e-3)+' kbar, opx')
    # plt.plot(T - 273, w_ol, c=c, label=str(p*1e-3) + ' kbar, ol', ls='--')
plt.xlabel('T (C)')
plt.ylabel('water solubilty (ppm)')
# plt.ylim((100, 1700))
plt.ylim((100, 25000))
plt.legend()


# test aluminious opx from bolfan-cassanova and keppler review figure 11
deltaH = -79.685  # kJ/mol
deltaV = 11.3  # cm3/mol
n = 0.5
A = 0.042  # ppm/bar
p_bar = np.linspace(20, 130)*1e3
T = np.array([600, 800, 1000, 1200, 1400]) + 273
plt.figure()
for TT in T:
    f_h2o = fugacity(p=p_bar*1e5, T=TT)  # Pa
    c_opx = water_capacity_thermo(p_bar*1e5, TT, f_h2o, A, n, deltaH=deltaH, deltaV=deltaV)
    c_ol = water_capacity_thermo(p_bar * 1e5, TT, f_h2o, A=0.0066, n=1, deltaH=0, deltaV=10.6)
    c_ol_D21 = water_capacity(p_bar*1e5, TT, f_h2o, X_Fe=0, n=0.6447, a=0, b=4905.5403, c=0, d=0)*1e-6
    D_w = c_ol/c_opx
    # plt.plot(p_bar*1e-3, D_w, label=str(TT - 273)+' C')
    plt.plot(p_bar * 1e-3, c_ol, label=str(TT - 273) + ' C')
    plt.plot(p_bar * 1e-3, c_ol_D21, label=str(TT - 273) + ' C', ls='--')
plt.xlabel('p (kbar)')
# plt.ylabel('$P_{ol/opx}$')
plt.ylabel('water content (ppm)')
# plt.ylim((0, 1))
plt.legend()

plt.show()

