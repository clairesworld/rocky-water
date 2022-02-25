import numpy as np
import fugacity as fug
from scipy.interpolate import interp1d, splev
from useful_and_bespoke import find_nearest_idx, colorize, iterable_not_string
from scipy import optimize
import matplotlib.pyplot as plt
import pandas as pd

Rb = 8.31446261815324
TO = 1.335e21  # kg
python_path = '/home/claire/Works/rocky-water/py/'


def water_capacity(p, T, f_h2o, X_Fe=0, n=0, a=0, b=0, c=0, d=0):
    """
    water storage capcity of NAMs in H2O-saturated conditions
    p in Pa, T in K, f_h2o (water fugacity) in Pa, X_Fe (iron content) in mole fraction
    for alpha-olivine returns ppm and wd, ring returns wt % (annoying!)
    """
    if n == 0:
        return np.exp(a + (b + c * p * 1e-9 + d * X_Fe) / T)
    # return np.exp(np.log(1 / (2 * n)) + a + n / 2 * np.log(f_h2o * 1e9) + (b + c * p * 1e9 + d * X_Fe) / T)
    try:
        return np.exp(a + n / 2 * np.log(f_h2o * 1e-9) + (b + c * p * 1e-9 + d * X_Fe) / T)
    except Exception as e:
        print('f_h2o\n', type(f_h2o), np.shape(f_h2o), '\n', f_h2o)
        print('p\n', type(p), np.shape(p), '\n', p)
        print('T\n', type(T), np.shape(T), '\n', T)
        print(n, a, b, c, d)
        raise e


def water_capacity_thermo(p, T, f_h2o, A, n, deltaH, deltaV):
    """ delta H in kJ/mol, deltaV in cm3/mol, A in ppm/bar^n, p in Pa, T in K, fh2o in Pa"""
    f_h2o = f_h2o * 1e-5  # Pa to bar
    deltaV = deltaV / (100 ** 3)
    deltaH = deltaH * 1e3
    c = A * f_h2o ** n * np.exp(-(deltaH + deltaV * p) / (Rb * T))
    return c * 1e-6  # weight fraction


def partitioning_coeff(p, T, f_h2o, A_j, deltaH_j, deltaV_j, n_j=0.5, A_i=0.0066, deltaH_i=0, deltaV_i=10.6, n_i=1):
    """
    A in ppm/bar, delta H in kJ/mol, delta V in cm3/mol, f, p and T in Pa and K
    """
    f_h2o = f_h2o * 1e-5  # to bar for use with A
    deltaV_j = deltaV_j / (100 ** 3)  # from cm3 to m3
    deltaV_i = deltaV_i / (100 ** 3)  # from cm3 to m3
    deltaH_i = deltaH_i * 1e3  # from kJ to J
    deltaH_j = deltaH_j * 1e3  # from kJ to J
    D_ij = A_i / A_j * f_h2o ** (n_i - n_j) * np.exp(
        -((deltaH_i - deltaH_j) + (deltaV_i - deltaV_j) * p) / (Rb * T))
    return D_ij  # unitless


def D_brg_ppv_T16(T):
    # from spline fit to D_ppv_pv in Townsend+ 2016 fig 4, T in K
    # assumes aluminious phase is present
    tck = (np.array([1002.18, 1002.18, 1002.18, 1002.18, 1032.75, 1061.14, 1074.24,
                     1089.52, 1102.62, 1126.64, 1192.14, 1281.66, 1310.04, 1355.9,
                     1403.93, 1506.55, 1598.25, 1657.21, 1722.71, 1766.38, 1825.33,
                     1910.48, 1991.27, 2063.32, 2189.96, 2303.49, 2408.3, 2543.67,
                     2604.8, 2733.62, 2901.75, 2901.75, 2901.75, 2901.75]),
           np.array([18.549, 18.26903335, 17.44538721, 16.50087189, 15.40357154,
                     15.10045558, 14.08808662, 13.07856277, 11.12218703, 10.09923983,
                     9.03574161, 8.60649523, 7.57796125, 6.81966252, 6.15216311,
                     5.55654442, 5.35695363, 4.97222845, 4.77911912, 4.41518883,
                     4.21920876, 3.937155, 3.73907767, 3.40345837, 3.28508823,
                     3.12097595, 2.95470437, 2.81790598, 2.67048685, 2.60612,
                     0., 0., 0., 0.]),
           3)

    D = splev(T, tck)
    return D ** -1  # inverse partitioning


def fugacity(p, T, func=fug.fugacity_frost, **kwargs):
    if np.size(T) > 1:  # iterate over T
        if np.size(p) == 1:
            p = [p] * len(T)
        f = []
        for pp, TT in zip(p, T):
            try:
                f.append(func(pp, TT, **kwargs))
            except Exception as e:
                print(pp * 1e-9, 'GPa', TT, 'K')
                raise e
        return np.array(f)
    elif np.size(p) > 1:  # iterate over p
        if np.size(T) == 1:
            T = [T] * len(p)
        f = []
        for pp, TT in zip(p, T):
            f.append(func(pp, TT, **kwargs))
        return np.array(f)
    else:
        return func(float(p), float(T), **kwargs)


def last_coexisting_idx(phase1, phase2, df):
    # get largest idx at which 2 phases coexist, requires df of composition
    X1 = df['X_' + phase1]
    X2 = df['X_' + phase1]
    idx = X1.where((X1 > 0) & (X2 > 0)).last_valid_index()
    print('idx last overlap of', phase1, 'and', phase2, ':', idx)
    return idx


def wmf_phase(phase, sat, df=None):
    if df is None:
        return None
    try:
        df.replace([np.inf, -np.inf, float('inf'), float('NaN'), np.nan], 0,
                   inplace=True)  # in case of infinite fugacity for olivine e.g.
        w = df['X_' + phase].multiply(sat, fill_value=0)  # actual water mass fraction given phase concentration
    except KeyError:
        return np.zeros(len(df))  # phase not present
    return w.to_numpy()  # actual water mass fraction per phase per layer


def sat_partitioning(phase_i, phase_j, sat_j, D_ji, df=None, na_val=1e-6, **kwargs):
    # only allow partitioning in overlapping region, use last overlap otherwise
    # input phase names (str), weight fractions X, sat water content wmf w_j, and partitioning D_ji = cj/ci
    if df is None:
        return None
    sat_i = np.zeros_like(sat_j)
    try:
        D_ji[0]
    except TypeError:
        D_ji = D_ji * np.ones_like(sat_j)
    try:
        X_i = df['X_' + phase_i]
        X_j = df['X_' + phase_j]
        idx = X_i.where((X_i > 0) & (X_j > 0)).last_valid_index()  # index of deepest layer with coexisting j & i
        sat_i[:idx + 1] = sat_j[:idx + 1] / D_ji[:idx + 1]
        sat_i[idx + 1:] = sat_j[idx] / D_ji[idx]  # should be constant
    except KeyError:
        if 'X_' + phase_i not in df.columns:
            print('X_' + phase_i, 'not in df columns')
            sat_i = np.zeros(len(df))  # this phase not present anywhere
        elif 'X_' + phase_j not in df.columns:
            print('X_' + phase_j, 'not in df columns - using na val')
            sat_i = na_val  # if no coexisting levels, prescribe a value
            print('sat_i', np.shape(sat_i))
    except Exception as e:
        print(e)
        sat_i = na_val  # if no coexisting levels, prescribe a value
    # TODO: check for phase i existing above phase j
    # w_i = sat_i * X_i
    return sat_i  # returns actual water content in that phase ->


def mineral_water_contents(p, T, X_Fe=0, df=None):  # p in Pa, T in K, X_Fe has no effect atch
    # TODO: allow for aluminious phase partitioning when there is plagioclase, spinel, or garnet present

    # try:
    f_h2o = fugacity(p, T)  # Pa
    f_h2o_PS = fugacity(p, T, func=fug.PSfugacity)  # Pa, Pitzer and Sterner 1994
    # except Exception as e:
    #     print('p, T', p, T)
    #     raise e

    """
    OLIVINE
    The fitted value of n may represent the average number of hydroxyls through multiple substitution mecha-
    nisms over the wide range of pressure and water concentration for olivine (Ferot & Bolfan-Casanova, 2012;
    Otsuka & Karato, 2011). Temperature dependence of the water storage capacity in olivine (b/T) can be
    attributed to the decrease in the activity of the H 2 O component in the aqueous fluid or silicate melt with
    increasing temperature (Mibe et al., 2002; Stalder et al., 2001). 
    """
    sat_ol = water_capacity(p, T, f_h2o, X_Fe, n=0.6447, a=0, b=4905.5403, c=0, d=0) * 1e-6  # in ppm wt?
    sat_corr_ol = sat_ol

    """
    WADSLEYITE & RINGWOODITE Mg2SiO4
    """
    sat_wad = water_capacity(p, T, f_h2o, X_Fe, n=0, a=-7.6356, b=13739.5371, c=0, d=0) * 1e-2  # in wt %?
    sat_ring = water_capacity(p, T, f_h2o, X_Fe, n=0, a=-6.8856, b=12206.2676, c=0, d=0) * 1e-2
    sat_corr_wad = sat_wad
    sat_corr_ring = sat_ring

    """
    ORTHOPYROXENE
    from Ol partitioning for now - but it might make more sense to just calculate this directly?
    'Moreover, water solubility in orthopyroxene is particularly well understood and calibrated.' - Keppler & B-C
    """
    # this is the 1:1 sum of pure enstatite and aluminious phases
    sat_ol_thermo = water_capacity_thermo(p, T, f_h2o_PS, A=0.0066, n=1, deltaH=0,
                                          deltaV=10.6)  # explicit parameterisation from Kohlstedt 1996
    sat_al_thermo = water_capacity_thermo(p, T, f_h2o_PS, A=0.042, n=0.5, deltaH=-79.685, deltaV=11.3)  # Al-bearing
    sat_en_thermo = water_capacity_thermo(p, T, f_h2o_PS, A=0.01354, n=1, deltaH=-4.563, deltaV=12.1)  # pure enstatite
    sat_opx_thermo = sat_al_thermo + sat_en_thermo
    D_ol_opx = sat_ol_thermo / sat_opx_thermo
    sat_opx = sat_ol / D_ol_opx  # this and others are just hypothetical values without considering (co-)stability
    sat_corr_opx = sat_partitioning(phase_i='opx', phase_j='ol', sat_j=sat_ol, D_ji=D_ol_opx, df=df,
                                    na_val=sat_opx_thermo)

    """ 
    CLINOPYROXENE, GARNET
    from Ol partitioning
    """
    # cpx (jadeite NaAlSi2O6) - TODO omphacite? see Smythe+ 1991
    D_ol_cpx = partitioning_coeff(p, T, f_h2o_PS, A_j=7.144, deltaH_j=0, deltaV_j=8.019, n_j=0.5)
    sat_cpx = sat_ol / D_ol_cpx
    sat_corr_cpx = sat_partitioning(phase_i='cpx', phase_j='ol', sat_j=sat_ol, D_ji=D_ol_cpx, df=df,
                                    na_val=water_capacity_thermo(p, T, f_h2o, A=7.144, deltaH=0, deltaV=8.019, n=0.5))

    # gt (pyrope Mg₃Al₂Si₃O₁₂)
    D_ol_gt = partitioning_coeff(p, T, f_h2o_PS, A_j=0.679, deltaH_j=0, deltaV_j=5.71, n_j=0.5)
    sat_gt = sat_ol / D_ol_gt
    sat_corr_gt = sat_partitioning(phase_i='gt', phase_j='ol', sat_j=sat_ol, D_ji=D_ol_gt, df=df,
                                   na_val=water_capacity_thermo(p, T, f_h2o, A=0.679, n=0.5, deltaH=0, deltaV=5.71))

    # spinel?? MgAl2O4 - upper mantle, might want to find D - e.g. via - D_spinel/melt / D_olivine/melt
    # sp is "among the most water poor minerals" but could be important when co-existing with opx
    # sat_sp = 1e-6   # I made this up for now

    # print('D_ol_opx', D_ol_opx)
    # print('D_ol_cpx', D_ol_cpx)
    # print('D_ol_gt', D_ol_gt)

    # TODO: what do we do if we have no olivine? prob best to use en

    """
    SILICATE PEROVSKITES
    Calcium silicate perovskite may have a higher water storage capacity (Chen et al., 2020), especially at the P-T conditions
    of the uppermost lower mantle (Muir & Brodholt, 2018). However, this high-pressure mineral converts to an
    amorphous phase upon pressure release, so its water storage capacity has not been quantitatively constrained
    (Chen et al., 2020; Nemeth et al., 2017). The water storage capacity in ferropericlase is likely to be relatively low
    throughout the lower mantle (Bolfan-Casanova et al., 2002; Litasov, 2010).
    """
    # Bridgmanite MgSiO3 >>> this is half of Earth's mantle mass so maybe the most important?
    D_ring_bdg = 15  # +- 8, Inoue+ 2010, at 1873 K
    D_ring_pv = D_ring_bdg
    # use T at rw-bdg phase boundary
    sat_bdg = sat_ring / D_ring_pv
    sat_pv = sat_bdg  # for now not looking at Fe-perovskite
    sat_corr_pv = sat_partitioning(phase_i='pv', phase_j='ring', sat_j=sat_ring, D_ji=D_ring_pv, df=df,
                                   na_val=100e-6)  # na is guess for now, like FITR experiments in Dong+

    # Ca-perovskite
    sat_capv = 10e-6  # assumed in Dong+ 2021
    sat_corr_capv = sat_capv

    # ferropericlase i.e. magnesiowustite iron endmember
    sat_fp = 10e-6  # assumed in Dong+ 2021
    sat_wus = sat_fp  # stixrude wustite model is solution of magnesio-wuestite-na2al2o4
    sat_corr_wus = sat_wus

    """ 
    POSTPEROVSKITE
    partitioning with bridgmanite calculated in Townsend+ 2016
    water partitioning controlled by presence of Al because different substitution mechanisms
    Dong+ paper quotes D ~ 2--18 so presumably this means they assume Al is present - use this for now
    >>> do you need to calculate intersection of partitioning T-dependence curve with phase boundary and use that value?
    """
    # # Al-free postperovskite - same Mg substitution favours Bdg
    # D_ppv_brg = ppv_partitioning(p_geo=p, T_geo=T, file='D_ppv_brg.csv', res=10000)
    # sat_ppv = sat_bdg * D_ppv_brg

    # Al-bearing ppv - Si-AlH defect favours Ppv
    # D_ppv_brg = ppv_partitioning(p_geo=p, T_geo=T, file='D_Alppv_brg.csv', res=10000)
    D_pv_ppv = D_brg_ppv_T16(T)
    # print('D_pv_ppv', D_pv_ppv)
    sat_ppv = sat_pv / D_pv_ppv
    # need to use sat_corr_pv here because that's the actual capacity concentration available
    sat_corr_ppv = sat_partitioning(phase_i='ppv', phase_j='pv', sat_j=sat_corr_pv, D_ji=D_pv_ppv, df=df, na_val=1e-6)
    # TODO check if there is a stable Al-bearing phase


    """
    other transition zone and lower mantle phases by partitioning
    """
    # high pressure clinoenstatite - this is called C2/C in Stixrude solubility model
    # todo - in the region where Cen coexists with wd should you use that partitioning? (D_wd_cen ~ 3.8)
    D_ol_hpcpx = 0.8  # Withers+ 2007
    sat_hpcpx = sat_ol / D_ol_hpcpx
    sat_corr_hpcpx = sat_partitioning(phase_i='hpcpx', phase_j='ol', sat_j=sat_ol, D_ji=D_ol_hpcpx, df=df, na_val=1e-6)

    # akimotoite?? - upper mantle, might not intersect
    D_ring_aki = 21  # Keppler & Bolfan-Casanova Table 3
    sat_aki = sat_ring / D_ring_aki
    sat_corr_aki = sat_partitioning(phase_i='aki', phase_j='ring', sat_j=sat_ring, D_ji=D_ring_aki, df=df, na_val=1e-6)

    # stishovite??
    D_ring_st = 521  # Keppler & Bolfan-Casanova Table 3
    sat_st = sat_ring / D_ring_st
    sat_corr_st = sat_partitioning(phase_i='st', phase_j='ring', sat_j=sat_ring, D_ji=D_ring_st, df=df, na_val=1e-6)

    """ 
    compile 
    """
    if df is not None:
        df['f_h2o'] = f_h2o  # Pa
        df['c_h2o'] = 0
        filter_col = [col for col in df if col.startswith('X_')]  # get all columns corresponding to mole fractions
        for col in filter_col:
            phase = col[2:]
            try:
                df['sat_' + phase] = eval('sat_' + phase)  # (hypothetical) saturation water mass fraction at this T, p
                df['sat_corr_' + phase] = eval('sat_corr_' + phase)  # accounting for co-stability when partitioning
                df['w_' + phase] = wmf_phase(phase, df['sat_corr_' + phase],
                                             df)  # actual water mass fraction given phase concentration
                df['c_h2o'] = df['c_h2o'].add(
                    df['w_' + phase], fill_value=0)  # total water concentration, sum of phase concentrations adds to 1
            except NameError:
                print('>>>MISSING PHASE IN SATURATION MODEL:', phase)
                # set missing to v small for now
                df['w_' + phase] = 1e-9
            except SyntaxError as e:
                phase = phase.replace('.', '')
                print('>>>MISSING PHASE IN SATURATION MODEL:', phase)
                # set missing to v small for now
                df['w_' + phase] = 1e-9
        for col in filter_col:  # do this after once you've calculated total layer water mass
            phase = col[2:]
            df['frac_h2o_' + phase] = df['w_' + phase] / df['c_h2o']  # proportion of layer's water in this phase
        df['mass_h2o(kg)'] = df['c_h2o'] * df['mass(kg)']  # convert to kg

        return df


def total_water(df):
    return np.sum(df['mass_h2o(kg)']) / np.sum(df['mass(kg)'])  # denominator is mass of that layer


def get_geotherm_file(file='geotherm_Tp1600.csv', path=python_path, res=100):
    df_read = pd.read_csv(path + file, names=['T_K', 'P_GPa'], index_col=False, skiprows=6)
    TP = interp1d(df_read['P_GPa'].to_numpy(), df_read['T_K'].to_numpy(), kind='linear')
    p_geo = np.linspace(np.min(df_read['P_GPa'].to_numpy()), np.max(df_read['P_GPa'].to_numpy()), num=res)
    T_geo = TP(p_geo)  # interpolated T as a function of p in GPa
    p_geo = p_geo * 1e9  # convert to Pa
    return T_geo, p_geo

#
# def ppv_partitioning(p_geo=None, T_geo=None, file='D_ppv_brg.csv', path=python_path, res=10000, plot=False):
#     # p in Pa, T in K
#     import pandas as pd
#
#     # find where phase boundary intersects geotherm
#     df_phase = pd.read_csv(path + 'ppv_phase_boundary.csv', names=['P_GPa', 'T_K'], index_col=False, skiprows=6)
#     p_prime = np.linspace(np.min(p_geo), np.max(p_geo), num=res)*1e-9
#     TP_geo = interp1d(np.array(p_geo)*1e-9, T_geo, kind='linear')
#     # f = TP_geo(p_prime)  # geotherm function
#
#     TP_phase = interp1d(df_phase['P_GPa'].to_numpy(), df_phase['T_K'].to_numpy(), kind='linear')
#     # g = TP_phase(np.linspace(np.min(df_phase['P_GPa'].to_numpy()), np.max(df_phase['P_GPa'].to_numpy()), num=res))  # phase bdy fn
#
#     df_D = pd.read_csv(path + file, names=['T_K', 'D'], index_col=False, skiprows=6)
#     D = interp1d(df_D['T_K'].to_numpy(), df_D['D'].to_numpy(), kind='linear')
#
#     try:
#         p_intersect = optimize.fsolve(lambda x: TP_geo(x) - TP_phase(x), p_prime[int(len(p_prime)/2)])
#         T_trans = TP_geo(p_intersect)
#         D_trans = D(T_trans)  # use this temperature to get partitioning coefficient
#         # print('T, D at phase boundary', T_trans, D_trans)
#     except ValueError as e:
#         # no intersection
#         T_trans = 0
#         D_trans = 0
#         print('no ppv phase!     - > ', e)
#         plot = True
#
#     if plot:
#         fig, axes = plt.subplots(1, 2)
#         axes[0].plot(p_prime, TP_geo(p_prime), 'k', label='geotherm')
#         axes[0].plot(df_phase['P_GPa'].to_numpy(), df_phase['T_K'].to_numpy(), 'k--', label='phase boundary')
#         axes[0].axhline(T_trans, c='r', label='intersection')
#         axes[0].set_xlabel('p (GPa)')
#         axes[0].set_ylabel('T (K)')
#         axes[1].plot(df_D['T_K'].to_numpy(), df_D['D'].to_numpy(), 'b')
#         axes[1].scatter(T_trans, D_trans, c='r', label='intersection')
#         axes[1].set_xlabel('T (K)')
#         axes[1].set_ylabel('D ppv-brg')
#         for ax in axes:
#             ax.legend()
#         plt.tight_layout()
#     return D_trans
#
