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


def partitioning_coeff(p, T, f_h2o, A_j, deltaH_j, deltaV_j, n_j=0.5, A_i=0.0066, deltaH_i=0, deltaV_i=10.6, n_i=1,
                       correction=1):
    """
    A in ppm/bar, delta H in kJ/mol, delta V in cm3/mol, f, p and T in Pa and K
    correction would refer to olivine (default i phase) which is from Kohlstedt 96, underestimated
    """
    f_h2o = f_h2o * 1e-5  # to bar for use with A
    deltaV_j = deltaV_j / (100 ** 3)  # from cm3 to m3
    deltaV_i = deltaV_i / (100 ** 3)  # from cm3 to m3
    deltaH_i = deltaH_i * 1e3  # from kJ to J
    deltaH_j = deltaH_j * 1e3  # from kJ to J
    # D_ij = A_i / A_j * f_h2o ** (n_i - n_j) * np.exp(
    #     -((deltaH_i - deltaH_j) + (deltaV_i - deltaV_j) * p) / (Rb * T))

    c_i = A_i * f_h2o ** n_i * np.exp(-(deltaH_i + deltaV_i * p) / (Rb * T)) * correction
    c_j = A_j * f_h2o ** n_j * np.exp(-(deltaH_j + deltaV_j * p) / (Rb * T))
    D_ij = c_i / c_j
    return D_ij  # unitless


def D_brg_ppv_T16(T):
    # from spline fit to D_ppv_pv in Townsend+ 2016 fig 4 bottom right, T in K
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
    return D ** -1  # inverse partitioning because they plot D_ppv_brg


def interp_PS03(p, T):  # Pa, K
    """ extract water capacity of stishovite from Panero & Stixrude 2003 figure 4
    note p and T will either both be scalars,  or both 1d arrays of same length, to return 1D array
    """
    if np.size(p) == 1 and np.size(T) == 1:
        p = [p]
        T = [T]
    elif np.size(p) != np.size(T):
        raise NotImplementedError('stishovite water capacity interpolation: p and must be same length')
    c = np.zeros(len(p))  # in wt%

    tck_0GPa = (np.array([1162.88, 1162.88, 1162.88, 1162.88, 1984.85, 1984.85, 1984.85, 1984.85]),
                np.array([-0.00124267,  0.00718062,  0.01839298,  0.14147842,  0., 0.,  0.,  0.]), 3)  # 0 GPa
    tck_25GPa = (np.array([1000.  , 1000.  , 1000.  , 1000.  , 1986.11, 1986.11, 1986.11, 1986.11]),
                 np.array([0.01921425, 0.04020424, 0.32797362, 1.13612814, 0.        , 0.        , 0.        , 0.        ]), 3)
    tck_60GPa =  (np.array([ 998.737,  998.737,  998.737,  998.737, 1614.9  , 1614.9  , 1614.9  , 1614.9  ]),
                  np.array([0.1303189 , 0.28071091, 0.76647274, 1.49226089, 0.        ,0.        , 0.        , 0.        ]), 3)

    # # fit spline
    # from scipy.interpolate import splrep
    # xy0 = np.loadtxt('PS03_fig4_60GPa.csv', delimiter=",")
    # x0 = xy0[:, 0]
    # y0 = xy0[:, 1]
    # spl0 = splrep(x0, y0, s=1)
    # print('spl0\n', spl0)

    for ii, pi in enumerate(p):
        if pi <= 25e9:  # between 0 and 25 GPa
            tck0, tck1 = tck_0GPa, tck_25GPa
            p0, p1 = 0, 25
        elif pi <= 60e9:
            tck0, tck1 = tck_25GPa, tck_60GPa
            p0, p1 = 25, 60
        else:
            c[ii] = 3.1111082685244025  # don't extrapolate beyond 2000 K  - otherwise quickly get to 10s of wt% at higher p
            continue

        # evaluate spline
        y0 = splev(T[ii], tck0)  # value at p0, T
        y1 = splev(T[ii], tck1)  # value at p1, T

        # linear interpolation
        m = (y1 - y0) / (p1 - p0)
        b = y1 - m * p1
        c[ii] = m * (pi * 1e-9) + b  # p in GPa
        # print(pi*1e-9, 'GPa', T[ii], 'K', c[ii], 'wt%')

    if np.size(c) == 1:
        return np.squeeze(c) * 1e-2  # wt frac
    else:
        return c * 1e-2  # wt frac


def fugacity(p, T, func=fug.fugacity_frost, **kwargs):
    # p in Pa, T in K?
    if np.size(T) > 1:  # iterate over T
        if np.size(p) == 1:
            p = [p] * len(T)
        f = []
        for pp, TT in zip(p, T):
            try:
                fh2o = func(pp, TT, **kwargs)
                if np.isinf(fh2o):
                    fh2o = 1e307  # note fugacity won't come into partitioning or water capacity equations at these pressures, this is just to avoid NoneType warnings
                f.append(fh2o)
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


def check_partitioning(phase_i, phase_j, sat_j, D_ji, df=None, na_val=1e-6, **kwargs):
    # only allow partitioning in overlapping region, use last overlap otherwise
    # input phase names (str), weight fractions X, sat water content wmf w_j, and partitioning D_ji = cj/ci
    if df is None:
        raise NotImplementedError('function saturation.check_partitioning requires input df')
    try:
        D_ji[0]
    except TypeError:
        D_ji = D_ji * np.ones(len(df))
    try:
        sat_j[0]
    except TypeError:
        sat_j = sat_j * np.ones(len(df))
    sat_i = np.zeros_like(sat_j)
    try:
        X_i = df['X_' + phase_i]
        X_j = df['X_' + phase_j]
        idx = X_i.where((X_i > 0) & (X_j > 0)).last_valid_index()  # index of deepest layer with coexisting j & i
        # print(phase_j, phase_i, 'idx', idx)
        sat_i[:idx + 1] = sat_j[:idx + 1] / D_ji[:idx + 1]
        sat_i[idx + 1:] = sat_j[idx] / D_ji[idx]  # should be constant
    except KeyError:
        if 'X_' + phase_i not in df.columns:
            sat_i = np.zeros(len(df))  # desired phase not present anyways, doesn't get used
        elif 'X_' + phase_j not in df.columns:
            if phase_j in ('ol', 'ring'):
                print(phase_j, 'not in Perple_X composition. Using fictive value to partition with', phase_i)
                # note hypothetical saturation will be given in sat_j
                if phase_j == 'ol':
                    idx_j_last = np.argmax(df['P(bar)'] > 14.5e4)  # don't extrapolate olivine beyond Dong+ fit range
                elif phase_j == 'ring':
                    idx_j_last = np.argmax(df['P(bar)'] > 23e4)  # don't extrapolate olivine beyond Dong+ fit range
                sat_i[:idx_j_last] = sat_j[:idx_j_last] / D_ji[:idx_j_last]
                sat_i[idx_j_last:] = sat_i[idx_j_last - 1]  # extend constant value
            else:
                print('Problem attempting to partition water:', phase_j, 'not in Perple_X composition. Prescribing fixed value for', phase_i)
                sat_i = na_val  # if no coexisting levels, prescribe a value
    except TypeError as e:
        if phase_i == 'pv':  # ring-pv transition may be abrupt and not overlap, in which use boundary values
            idx_j_last = df['X_' + phase_j].to_numpy().nonzero()[0][-1]  # deepest ring layer
            idx_i_first = df['X_' + phase_i].to_numpy().nonzero()[0][0]  # shallowest pv layer
            sat_i[:] = sat_j[idx_j_last]/ D_ji[idx_j_last]  # should be constant
        elif idx is None and phase_j == 'ol':
            # in rare cases can get a tiny olivine shell but doesn't overlap with cpx, gt etc.
            print('No overlap of', phase_j, 'and', phase_i, 'in Perple_X composition. Using fictive value for', phase_j, 'to partition with', phase_i)
            idx_ol_last = np.argmax(df['P(bar)'] > 14.5e4)  # don't extrapolate olivine beyond Dong+ fit range
            sat_i[:idx_ol_last] = sat_j[:idx_ol_last] / D_ji[:idx_ol_last]
            sat_i[idx_ol_last:] = sat_i[idx_ol_last - 1]  # extend constant value
        elif idx is None and phase_i == 'aki':
            # rogue aki deep layer
            idx_j_last = df['X_' + phase_j].to_numpy().nonzero()[0][-1]  # deepest ring layer
            sat_i[:] = sat_j[idx_j_last]/ D_ji[idx_j_last]  # should be constant
        elif idx is None and phase_i == 'ppv':
            # pv-ppv transition is abrupt
            idx_j_last = df['X_' + phase_j].to_numpy().nonzero()[0][-1]  # deepest pv layer
            idx_i_first = df['X_' + phase_i].to_numpy().nonzero()[0][0]  # shallowest ppv layer
            sat_i[:] = sat_j[idx_j_last]/ D_ji[idx_j_last]  # should be constant
        else:
            print(phase_i, phase_j, 'idx', idx)
            raise e
    except IndexError as e:
        print('IndexError line 199-203 (phase_i =', phase_i, ', phase_j =', phase_j)
        print('sat_i', sat_i)
        print('sat_i', sat_j)
        print('D_ji', D_ji)
        raise e

    return sat_i  # returns coexistance-corrected water saturation in phase i ->


def mineral_water_contents(p, T, X_Fe=0, df=None):  # p in Pa, T in K, X_Fe has no effect atch
    K96_factor = 3  # Mookherjee & Karato use 3.5; use 1 to match Dong & Fischer

    # double-checked and fugacity calculations are accurate, 11/05/22
    # FW is undefined above ~112 GPa (deep in Pv stability field, doesn't matter for calcs)
    f_h2o = fugacity(p, T)  # Pa
    f_h2o_PS = fugacity(p, T, func=fug.PSfugacity)  # Pa, Pitzer and Sterner 1994 - used in K96 for ol

    """
    OLIVINE
    The fitted value of n may represent the average number of hydroxyls through multiple substitution mecha-
    nisms over the wide range of pressure and water concentration for olivine (Ferot & Bolfan-Casanova, 2012;
    Otsuka & Karato, 2011). Temperature dependence of the water storage capacity in olivine (b/T) can be
    attributed to the decrease in the activity of the H 2 O component in the aqueous fluid or silicate melt with
    increasing temperature (Mibe et al., 2002; Stalder et al., 2001). 
    TODO try Ferot & Bolfan-Casanova 2012 param
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
    Mierdel 2007 experiments, 1.5 to 3.5 GPa and 1073 to 1373 K
    """
    # enstatite from Rauch & Keppler, aluminous from Mierdel+
    # both use PS fugacity, and both Bell95 & Paterson calib reported in Tables but seems like Bell is the fit
    # this is the 1:1 sum of pure enstatite and aluminious phases
    sat_ol_thermo = water_capacity_thermo(p, T, f_h2o_PS, A=0.0066, n=1, deltaH=0,
                                          deltaV=10.6) * K96_factor  # explicit parameterisation from Kohlstedt 1996
    sat_al_thermo = water_capacity_thermo(p, T, f_h2o_PS, A=0.042, n=0.5, deltaH=-79.685, deltaV=11.3)  # Al-bearing
    sat_en_thermo = water_capacity_thermo(p, T, f_h2o_PS, A=0.01354, n=1, deltaH=-4.563, deltaV=12.1)  # pure enstatite - assume other end-member same
    sat_opx_thermo = sat_al_thermo + sat_en_thermo
    D_ol_opx = sat_ol_thermo / sat_opx_thermo
    sat_opx = sat_ol / D_ol_opx  # this and others are just hypothetical values without considering (co-)stability
    sat_corr_opx = check_partitioning(phase_i='opx', phase_j='ol', sat_j=sat_ol, D_ji=D_ol_opx, df=df,
                                      na_val=sat_opx_thermo)

    """ 
    CLINOPYROXENE, GARNET, other upper mantle minerals
    from Ol partitioning
    """
    # cpx (tehnically for jadeite NaAlSi2O6 but fits other cpx - note Bromily & Keppler uses Peterson82 calib for ol)
    # K96 uses Pitzer & Sterner fh2o, B&K no mention...
    D_ol_cpx = partitioning_coeff(p, T, f_h2o_PS, A_j=7.144, deltaH_j=0, deltaV_j=8.019, n_j=0.5, correction=1)
    sat_cpx = sat_ol / D_ol_cpx
    sat_corr_cpx = check_partitioning(phase_i='cpx', phase_j='ol', sat_j=sat_ol, D_ji=D_ol_cpx, df=df,
                                      na_val=water_capacity_thermo(p, T, f_h2o, A=7.144, deltaH=0, deltaV=8.019, n=0.5))

    # high pressure clinoenstatite - note this is called C2/C in Stixrude solubility model
    D_ol_hpcpx = 0.7  # Withers & hirschmann 2007 figure 8 "weighted average"
    D_wd_hpcpx = 3.4  # Withers & hirschmann 2007 figure 8 "weighted average" todo use this for overlap
    sat_hpcpx = sat_ol / D_ol_hpcpx
    sat_corr_hpcpx = check_partitioning(phase_i='hpcpx', phase_j='ol', sat_j=sat_ol, D_ji=D_ol_hpcpx, df=df,
                                        na_val=None)
    # At 13.4 GPa, we observe a single value of DWd/CEn of 2.8 ± 0.4 (Table 2). This is of similar magnitude but
    # slightly smaller than values between 3.3 and 5.2 observed for DWd/CEn in previous studies at 15–15.5 GPa
    # (Inoue et al. 1995; Bolfan-Casanova et al. 2000) (Fig. 8). - Withers & Hirschmann 2007 - because only a few layers
    # don't use this and introduce discontinuity - or maybe - note B-C uses Paterson calib.
    # safe to use single value because just a few GPa

    # gt (pyrope-rich, Mg₃Al₂Si₃O₁₂)
    # this will result in transition zone values of about 1000 ppm at 1600 K Tp. Liu+ 2021 get a bit lower.
    # D_ol_gt = partitioning_coeff(p, T, f_h2o_PS, A_j=0.679, deltaH_j=0, deltaV_j=5.71, n_j=0.5, correction=K96_factor)
    D_ol_gt = 1  # Novella+ 14, Ferot & B-C,
    sat_gt = sat_ol / D_ol_gt
    sat_corr_gt = check_partitioning(phase_i='gt', phase_j='ol', sat_j=sat_ol, D_ji=D_ol_gt, df=df,
                                     na_val=water_capacity_thermo(p, T, f_h2o, A=0.679, n=0.5, deltaH=0, deltaV=5.71))

    # spinel?? MgAl2O4 - narrow in upper mantle, might want to find D - e.g. via - D_spinel/melt / D_olivine/melt
    # sp is "among the most water poor minerals" (Keppler & B-C)
    sat_sp = 1e-6
    sat_corr_sp = sat_sp

    # plagioclase - make something tiny up because it will literally make up 1% of a couple of layers at the top
    sat_pl = 1e-6
    sat_an = sat_pl  # not sure why not showing up in solution?
    sat_corr_an = sat_pl

    # akimotoite - how often is this important? rarely--- often like one tiny layer
    D_ring_aki = 21  # Keppler & Bolfan-Casanova Table 3
    sat_aki = sat_ring / D_ring_aki
    sat_corr_aki = check_partitioning(phase_i='aki', phase_j='ring', sat_j=sat_ring, D_ji=D_ring_aki, df=df,
                                      na_val=None)

    """
    SILICATE PEROVSKITES
    """
    # Dong+ note SIMS and FTIR for pv don't mesh (factor of 2?). argues to use partitioning with ring rather than experiment as
    # most experiments are only at ~25 GPa, single value with no TP-dependence#
    # Inoue direct partitioning experiments (D~15) gives around 500 ppm c_pv through the lower mantle (~1800 K) - SIMS so maybe hydrous inclusions
    # Panero+2020 give ~40 ppm at Tp=1600 K
    # problem with past works is FTIR/SIMS picking up fluid/hydrous phase inclusions/impurities. need large crystals for
    # good FTIR data. Liu argues that higher solubilities (Murakami+ 2002) are this because they don't see broad bands
    # in their large pure crystals
    # recently Liu+ 2021 try to address discrepancies with good quality single crystal at 24-26 GPa, 1700-1900 K
    # "Hernández et al. (2013), using first-principles calculations, computed the H2O partition coefficient between pure
    # MgSi2O4 ringwoodite, MgO periclase, and MgSiO3 bridgmanite, and found a ringwoodite–bridgmanite H2O partition
    # coefficient of about 10:1." - Townsend+

    # D = 21 gives 300 ppm at 1850 K ish, meanwhile Liu21 get 15-30 at 1900
    # D_ring_pv = 21  # 15 +- 8, Inoue+ 2010, at 1873 K (bdg) - Dong use 15, 21 is more in line with Liu
    # sat_pv = sat_ring / D_ring_pv
    # sat_corr_pv = check_partitioning(phase_i='pv', phase_j='ring', sat_j=sat_ring, D_ji=D_ring_pv, df=df,
    #                                  na_val=None)  # na is guess for now, like FITR experiments in Dong+
    sat_pv = 30e-6  # Liu+21, 1900 K - upper limit at 1900 K - higher temps would be drier??
    sat_corr_pv = sat_pv

    # Ca-perovskite
    # Calcium silicate perovskite may have a higher water storage capacity (Chen et al., 2020), especially at the P-T conditions
    #     of the uppermost lower mantle (Muir & Brodholt, 2018). However, this high-pressure mineral converts to an
    #     amorphous phase upon pressure release, so its water storage capacity has not been quantitatively constrained
    #     (Chen et al., 2020; Nemeth et al., 2017). - do they mean hard to see contamination of hydrous phases?
    # Chen, Litasov&Ohtani2003, Murakami+2002 say important but Panero+2020 finds it dry
    # Murakami+2003 - 0.3-0.4 wt%
    # Litasov & Ohtani: approx 5 wt% at 1900 K, 25 GPa
    # Chen+2020: synthesised 19-120 GPa and 1400-2200 K - get approx 0.5-1 wt% - FTIR measurements and can't rule out contamination - propse Ca defect mechanism?
    # uncertainty because dry cubic capv structure (500 K to mantle temps) not seen in these exps, instead colder tetragonal structure still stable in presence of H2O, might fuck up thermodynamic eq.
    sat_capv = 0.5e-2  # low estimate from Chen+2020, 10e-6 assumed in Dong+ 2021 - probably want to test this
    sat_corr_capv = sat_capv

    # ferropericlase i.e. magnesiowustite iron endmember
    # note stixrude wustite model is solution of magnesio-wuestite-na2al2o4
    # Panero+2020 finds it dry but Bolfan-Casanova+ 2000 says wetter than Pv
    # Liu+2021 also find it <100 ppm; "Ferropericlase, another major mineral in this region, can contain no more than 50
    # ppm wt. H2O (Fig. 6) (Bolfan-Casanova et al., 2002; Litasov, 2010)."
    # Merli+ 2016: first principles, D_fp_pv = 0.31–0.56 would be in line with ~30 ppm brg and ~10 ppm fp as above
    # Litaov & Ohtani - 11 ppm periclase, 14 ppm ferropericlase
    sat_wus = 10e-6  # assumed in Dong+ 2021 - this is in Litasov 2010
    sat_corr_wus = sat_wus
    sat_wuls = sat_wus  # low spin state phase - assume the same
    sat_corr_wuls = sat_corr_wus

    # from Liu2021:  Even if we consider the case of some amount of water in davemaoite (CaSiO3 perovskite)
    # (Chen et al., 2020) and hydrous phases such as phase D and phase H in this region (e.g., Ohtani et al., 2004;
    # Nishi et al., 2014), ambient lower-mantle temperatures (Katsura et al., 2010) are too high for the stability of
    # hydrous minerals (Walter et al., 2015). The presence of metallic iron would cause the top of the lower mantle to
    # be more dry because iron would react with hydroxyl in the hydrous and nominally anhydrous minerals (Zhu et al.,
    # 2019). We therefore conclude that the majority of the top region of a pyrolitic lower mantle is nearly dry."

    """ 
    POSTPEROVSKITE
    partitioning with bridgmanite calculated in Townsend+ 2016
    water partitioning controlled by presence of Al because different substitution mechanisms
    Dong+ paper quotes D ~ 2--18 so presumably this means they assume Al-free
    
    todo: do you need to sum both defects in fig 4? otherwise, assuming just a single type of defect. 
    
    "The above formulation assumes that defects are partitioning between two solid phases in the absence of any fluids. 
    As with previous work (cf. Hernández et al., 2013), the results do not represent water storage capacity which is 
    defined as the amount of water in a nominally anhydrous mineral in equilibrium with a hydrous melt (Kohlstedt et al., 1996)."
    """
    # Al-free postperovskite - Mg-H substitution favours Bdg

    # Al-bearing ppv - Si-AlH defect favours Ppv
    D_pv_ppv = D_brg_ppv_T16(T)  # get partition coefficient at all T
    sat_ppv = sat_corr_pv / D_pv_ppv
    # need to use sat_corr_pv here because that's the actual capacity concentration available (but currently constant/equal)
    # function below will use T at intersection to find D
    sat_corr_ppv = check_partitioning(phase_i='ppv', phase_j='pv', sat_j=sat_corr_pv, D_ji=D_pv_ppv, df=df, na_val=None)

    # fe-al perovskite - this is always in lowest mantle so assume same as aluminous ppv? - not significant, <1% wt of layer
    sat_fapv = sat_ppv
    sat_corr_fapv = sat_corr_ppv

    """
    SiO2 polymorphs
    """
    # stishovite
    # in the most extreme cases, st goes to about 120 GPa (2250 K at 1600 K Tp) then changes to seif
    # D_ring_st = 521  # Keppler & Bolfan-Casanova Table 3, this is for Al-free
    # sat_st = sat_ring / D_ring_st
    # sat_corr_st = sat_partitioning(phase_i='st', phase_j='ring', sat_j=sat_ring, D_ji=D_ring_st, df=df, na_val=1e-6)
    # sat_st = (4.1 - 0.042 * (p * 1e-9) + 45 * 1/T) * 1e-2  # Lin 2020 unpublished preprint, 32-52 GPa, to 1835 K
    # ^ should not use "hydrous SiO2" as this isn't captured by thermodynamic model.

    # also see Litasov+ 2007 who finds it rather high - 0.3 wt% at 20 GPa 1673 K (strong Al dependence)
    # Litasov+2007 find <30 ppm for Al-free (but let's assume this would not be natural) - note this matches B-C2000
    # Liu+21 and Panero+2003 is 300-400 ppm (25 GPa/1900 K
    sat_st = interp_PS03(p, T)  # to max of 3 wt% at 60 GPa
    sat_corr_st = sat_st

    # seifertite - above 120 GPa - Lin 2020 preprint Table 1 has data - not sure about T,P regression just take avg
    # sat_seif = (1.04 + 0.63 + 0.99 + 0.76)/ 4 * 1e-2
    sat_seif = 1e-2  # assume 1 wt% (line Lin 2020 preprint) - so more like hydrous phase but at these pressures it's all fucked anyways
    sat_corr_seif = sat_seif

    # coesite - various experiments covering its entire stability field (4-12 GPa and 900-2300 K) fit by Yan+2021 yay
    sat_coes = (-105 + 5.2 * (p * 1e-9) + 0.112 * (T - 273.15)) * 1e-6
    sat_corr_coes = sat_coes

    # note on thermodynamic stability of hydrous silica phases: if anything, seems like H2O reduces st-seif transition by "at least 20 GPa" -Nisr+2020

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
        for col in filter_col:  # do this after once you've calculated total layer water mass
            phase = col[2:]
            df['frac_h2o_' + phase] = df['w_' + phase] / df['c_h2o']  # proportion of layer's water in this phase
        df['mass_h2o(kg)'] = df['c_h2o'] * df['mass(kg)']  # convert to kg
        return df


def total_water_frac(df, i_min=0, i_max=-1):
    """ returns mass fraction of water WITH RESPECT TO layers involved (e.g., only mantle mass) """
    subset = df.iloc[i_min:i_max]  # select a chunk of rows if asked
    return np.sum(subset['mass_h2o(kg)']) / np.sum(subset['mass(kg)'])  # denominator is mass of that layer


def total_water_mass(df, i_min=0, i_max=None):
    """ returns mass of water in kg """
    # print('total_water_mass(), i_min', i_min, 'i_max', i_max)
    if i_max is None:
        subset = df.iloc[i_min:]  # select a chunk of rows if asked
    else:
        subset = df.iloc[i_min:i_max]  # select a chunk of rows if asked
    # print(subset)
    return np.sum(subset['mass_h2o(kg)'])




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
# def get_geotherm_file(file='geotherm_Tp1600.csv', path=python_path, res=100):
#     df_read = pd.read_csv(path + file, names=['T_K', 'P_GPa'], index_col=False, skiprows=6)
#     TP = interp1d(df_read['P_GPa'].to_numpy(), df_read['T_K'].to_numpy(), kind='linear')
#     p_geo = np.linspace(np.min(df_read['P_GPa'].to_numpy()), np.max(df_read['P_GPa'].to_numpy()), num=res)
#     T_geo = TP(p_geo)  # interpolated T as a function of p in GPa
#     p_geo = p_geo * 1e9  # convert to Pa
#     return T_geo, p_geo