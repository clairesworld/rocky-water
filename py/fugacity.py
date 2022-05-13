# -*- coding: utf-8 -*-
"""
Created on Thu Feb 03 12:12:46 2017

@author: tony.withers@uwo.ca

Functions to calculate H2O molar volume (PSvolume) and fugacity (PSfug)
using the Pitzer and Sterner equation of state.

Pitzer, K.S. and Sterner, S.M., 1994. Equations of state valid
continuously from zero to extreme pressures for H2O and CO2.
Journal of Chemical Physics. 101: 3111-3116.
"""

import math
from scipy import optimize
import numpy as np

coeff = []
coeff.append([0, 0, 0.24657688e6, 0.51359951e2, 0, 0])
coeff.append([0, 0, 0.58638965e0, -0.28646939e-2, 0.31375577e-4, 0])
coeff.append([0, 0, -0.62783840e1, 0.14791599e-1, 0.35779579e-3, 0.15432925e-7])
coeff.append([0, 0, 0, -0.42719875e0, -0.16325155e-4, 0])
coeff.append([0, 0, 0.56654978e4, -0.16580167e2, 0.76560762e-1, 0])
coeff.append([0, 0, 0, 0.10917883e0, 0, 0])
coeff.append([0.38878656e13, -0.13494878e9, 0.30916564e6, 0.75591105e1, 0, 0])
coeff.append([0, 0, -0.65537898e5, 0.18810675e3, 0, 0])
coeff.append([-0.14182435e14, 0.18165390e9, -0.19769068e6, -0.23530318e2, 0, 0])
coeff.append([0, 0, 0.92093375e5, 0.12246777e3, 0, 0])


def PSeos(volume, temperature, targetP):  # cc/mol, Kelvins, bars
    R = 8314510  # Pa.cc/K/mol
    den = 1 / volume  # mol/cc
    c = []
    for i in range(10):
        c.insert(i, coeff[i][0] * temperature ** -4 + coeff[i][1] * temperature ** -2
                 + coeff[i][2] * temperature ** -1 + coeff[i][3]
                 + coeff[i][4] * temperature + coeff[i][5] * temperature ** 2)
    try:
        pressure = (den + c[0] * den ** 2 - den ** 2 * ((c[2] + 2 * c[3] * den + 3 * c[4] * den ** 2
                                                         + 4 * c[5] * den ** 3) / (
                                                                c[1] + c[2] * den + c[3] * den ** 2 + c[4] * den ** 3
                                                                + c[5] * den ** 4) ** 2) + c[6] * den ** 2 * math.exp(
            -c[7] * den)
                    + c[8] * den ** 2 * math.exp(-c[9] * den)) * R * temperature / 1e5
    except OverflowError as e:
        return np.inf
    return pressure - targetP  # bars


def PSvolume(pressure, temperature):  # bars, Kelvins
    volume = optimize.root(PSeos, 10, args=(temperature, pressure))
    return volume.x


def PSfugacity(pressure, temperature):  # Pa, Kelvins
    pressure = pressure * 1e-5  # to bars
    R = 8314510  # Pa.cc/K/mol
    c = []
    for i in range(10):
        try:
            c.insert(i, coeff[i][0] * temperature ** -4 + coeff[i][1] * temperature ** -2
                 + coeff[i][2] * temperature ** -1 + coeff[i][3]
                 + coeff[i][4] * temperature + coeff[i][5] * temperature ** 2)
        except ValueError:
            temperature = float(temperature)
            c.insert(i, coeff[i][0] * temperature ** -4 + coeff[i][1] * temperature ** -2
                     + coeff[i][2] * temperature ** -1 + coeff[i][3]
                     + coeff[i][4] * temperature + coeff[i][5] * temperature ** 2)
    volume = PSvolume(pressure, temperature)
    den = 1 / volume  # mol/cc
    try:
        fug = math.exp(math.log(den) + c[0] * den + (1 / (c[1] + c[2] * den + c[3] * den ** 2
                                                          + c[4] * den ** 3 + c[5] * den ** 4) - 1 / c[1])
                       - c[6] / c[7] * (math.exp(-c[7] * den) - 1)
                       - c[8] / c[9] * (math.exp(-c[9] * den) - 1)
                       + pressure * 1e5 / (den * R * temperature)
                       + math.log(R * temperature) - 1)
    except (OverflowError, ValueError):
        # print('PS: undefined fugacity at', pressure * 1e5 * 1e-9, 'GPa')
        return np.inf
    return fug  # Pa


""" MODIFIED REDLICH-KWONG EOS from Frost & Wood 1997 """
import numpy.polynomial.polynomial as poly


def saturation_pressure(T):
    # T in K, returns Pa, from Holland & Powell 1991 (originally in kbar)
    return (-13.627e-3 + 7.29395e-7 * T ** 2 - 2.34622e-9 * T ** 3 + 4.83607e-15 * T ** 5) * 1e3 * 1e5


def fugacity_holland(p, T, p0=2.0, a0=1113.4, a4=-0.22291, a5=-3.8022e-4, a6=1.779e-7, a7=5.8487, a8=-2.1370e-2,
                     a9=6.8133e-5,
                     c0=-3.02565e-2, c1=-5.343144e-6, d0=-3.229755e-3, d1=2.221522e-6, b=1.465):
    """ T in K, p in Pa, Holland & Powell 1991 - all constants are in kJ and kbar so convert Rb"""
    p_kbar = p * 1e-5 * 1e-3
    Rb_kbar = 8.3145 * 1e-3
    psat = saturation_pressure(T) * 1e-5 * 1e-3
    # p_kbar = psat  # ??

    # RTlnf = Rb * T * np.log(1000 * p_kbar) + b * p_kbar + a / (b * np.sqrt(T)) * (
    #             np.log(Rb_kbar * T + b * p_kbar) - np.log(Rb_kbar * T + 2 * b * p_kbar)) \
    #         + 2 / 3 * c * p * np.sqrt(p_kbar) + d / 2 * p_kbar ** 2
    # return np.exp(RTlnf / (Rb_kbar * T))  # oops that was for co2

    # A = a / (b*Rb_kbar*T**1.5)
    # B = b*p_kbar/(Rb_kbar*T)
    # ln_gamma =
    print('   calculating fugacity with p =', p_kbar, 'kbar, T =', T, 'K')
    # return fug.PSfugacity(p_kbar*1e3, T)*1e5  # returns Pa

    if T <= 695:
        # a_gas
        a = a0 + a7 * (673 - T) + a8 * (673 - T) ** 2 + a9 * (673 - T) ** 3
    else:
        a = a0 + a4 * (T - 673) + a5 * (T - 673) ** 2 + a6 * (T - 673) ** 3

    # solve MRK equation
    zeros = poly.polyroots([p_kbar,
                            -Rb_kbar * T,
                            -(b * Rb_kbar * T + b ** 2 * p_kbar - a / np.sqrt(T)),
                            -a * b / np.sqrt(T)])
    zeros = zeros[~np.iscomplex(zeros)]  # should be only one root above critical point 647 K
    zeros = zeros[zeros > 0]  # should be only one root above critical point 647 K
    print('V_MRK =', zeros)
    V_init = Rb_kbar * T / p_kbar + b
    print('   V_init', V_init)
    if np.size(zeros) > 1:
        raise Exception('fugacity: too many roots for V! (p =', p_kbar, 'kbar, T =', T, 'K)')
    else:
        V_MRK = np.squeeze(zeros)

    # calculate molar volume
    c = c0 + c1 * T
    d = d0 + d1 * T
    V = V_MRK + c * np.sqrt(p_kbar - p0) + d * (p_kbar - p0)
    print('V =', V)

    # calculate fugacity coefficient
    z = p_kbar * V / (Rb_kbar * T)  # compressibility factor
    A = a / (b * Rb_kbar * T ** 1.5)
    B = b * p_kbar / (Rb_kbar * T)
    if p_kbar > p0:
        ln_gamma_virial = 1 / (Rb_kbar * T) * (2 / 3 * c * (p_kbar - p0) ** 1.5 + (d / 2 * (p_kbar - p0) ** 2))
    else:
        ln_gamma_virial = 0
    ln_gamma = z - 1 - np.log(z - B) - A * np.log(1 + B / z) + ln_gamma_virial
    print('ln_gamma_virial', ln_gamma_virial)
    print('ln_gamma', ln_gamma)
    print('z', z, 'A', A, 'B', B)
    print('z-B', z - B)
    print('1 + B/z', 1 + B / z)
    f = np.exp(ln_gamma) * p_kbar
    return f * 1e3 * 1e5  # Pa


def fugacity_holland_SI(p, T, p0=2.0, a0=1113.4, a4=-0.22291, a5=-3.8022e-4, a6=1.779e-7, a7=5.8487, a8=-2.1370e-2,
                        a9=6.8133e-5,
                        c0=-3.02565e-2, c1=-5.343144e-6, d0=-3.229755e-3, d1=2.221522e-6, b=1.465):
    """
    same as above but attempting to convert to SI first
    T in K, p in Pa, Holland & Powell 1991 - all constants are in kJ and kbar"""

    Rb = 8.3145
    psat = saturation_pressure(T)
    print('   calculating fugacity with p =', p * 1e-9, 'GPa, T =', T, 'K')

    if T <= 695:
        # a_gas
        a = a0 + a7 * (673 - T) + a8 * (673 - T) ** 2 + a9 * (673 - T) ** 3
    else:
        a = a0 + a4 * (T - 673) + a5 * (T - 673) ** 2 + a6 * (T - 673) ** 3

    # convert to SI
    a = a * 1e-2
    b = b * 1e-5
    p0 = p0 * 1e8

    # solve MRK equation
    zeros = poly.polyroots([-a * b / np.sqrt(T),
                            -(b * Rb * T + b ** 2 * p - a / np.sqrt(T)),
                            -Rb * T,
                            p])
    zeros = zeros[~np.iscomplex(zeros)]  # should be only one root above critical point 647 K
    zeros = zeros[zeros > 0]  # should be only one root above critical point 647 K
    print('V_MRK =', zeros, '=', zeros * ((10 ** 2) ** 3), 'cm3/mol')
    V_init = Rb * T / p + b
    print('   V_init', V_init)
    if np.size(zeros) > 1:
        raise Exception('fugacity: too many roots for V! (p =', p * 1e-9, 'GPa, T =', T, 'K)')
    else:
        V_MRK = np.squeeze(zeros)

    # calculate molar volume
    c = c0 + c1 * T
    d = d0 + d1 * T
    c = c * 1e-4  # convert to SI
    d = d * 1e-8
    print('c', c, 'd', c)
    V = V_MRK + c * np.sqrt(p - p0) + d * (p - p0)
    print('V_virial', c * np.sqrt((p - p0)) + d * (p - p0))
    print('V =', V)

    # calculate fugacity coefficient
    z = p * V / (Rb * T)  # compressibility factor
    A = a / (b * Rb * T ** 1.5)
    B = b * p / (Rb * T)
    if p > p0:
        ln_gamma_virial = 1 / (Rb * T) * (2 / 3 * c * (p - p0) ** 1.5 + (d / 2 * (p - p0) ** 2))
    else:
        ln_gamma_virial = 0
    ln_gamma = z - 1 - np.log(z - B) - A * np.log(1 + B / z) + ln_gamma_virial
    print('ln_gamma_virial', ln_gamma_virial)
    print('ln_gamma', ln_gamma)
    print('z', z, 'A', A, 'B', B)
    print('z-B', z - B)
    print('1 + B/z', 1 + B / z)
    f = np.exp(ln_gamma) * p
    return f  # Pa


def fugacity_frost(p, T,
                   a0=5.3957e8, a1=-6.362e5, a2=236.81,
                   b0=27.732, b1=-2.0179e-2, b2=9.2125e-6,
                   c0=-0.1244, c1=1.79e-4, c2=-7.8587e-8,
                   d0=2.186e-4, d1=-3.6836e-7, d2=1.6127e-10):
    """
    input/output units: Pa, K
    working units: cm3, bar, K, mol
     Frost & Wood 1997, geochem cosmochem acta"""
    p = p*1e-5

    R = 83.1446261815324  # cm3 bar /K /mol
    a = a0 + a1 * T + a2 * T ** 2
    b = b0 + b1 * T + b2 * T ** 2
    c = c0 + c1 * T + c2 * T ** 2
    d = d0 + d1 * T + d2 * T ** 2

    V = (R * T / p) + b - (a * R * np.sqrt(T)) / ((R * T + b * p) * (R * T + 2 * b * p)) + c * np.sqrt(p) + d * p
    RTlnf = R * T * np.log(p) + b * p + a / (b * np.sqrt(T)) * (
            np.log(R * T + b * p) - np.log(R * T + 2 * b * p)) + 2 / 3 * c * p * np.sqrt(p) + d / 2 * p ** 2
    f = np.exp(RTlnf / (R * T))
    # if np.isnan(f) or np.isinf(f):
    #     print('FW: undefined fugacity at', p*1e5*1e-9, 'GPa')
    return f*1e5  # Pa

