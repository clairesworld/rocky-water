import numpy as np
import math
import warnings

R_b = 8.3144598  # universal gas constant in J mol −1 K −1
G = 6.67408e-11

def debye(thetaT):
    """ from Lena Noack """
    val = 3 / thetaT ** 3 * IntEth(thetaT)
    return val


def IntEth(x):
    """ from Lena Noack """
    # int_0^x  tau^3/(exp tau - 1) dtau
    nr = 100  # number of intergation steps
    dx = x / nr
    IntVal = 0.0
    tau = 0.5 * dx
    for i in range(nr):
        IntVal = IntVal + dx * (tau ** 3 / (np.exp(tau) - 1))
        tau = tau + dx
    return IntVal


def BM3(V, T, n, T0, V0, M, K0, KP0, G0, GP0, theta0, gamma0, q0, etaS0):
    """ from Lena Noack """
    x = V0 / V
    f = 0.5 * (x ** (2.0 / 3.0) - 1.0)

    theta = theta0 * np.sqrt(1.0 + 6.0 * f * gamma0 + 3.0 * f ** 2 * gamma0 * (-2.0 + 6.0 * gamma0 - 3.0 * q0))
    gamma = (1.0 + 2.0 * f) * gamma0 * (1.0 + f * (-2.0 + 6.0 * gamma0 - 3.0 * q0)) * (theta0 / theta) ** 2

    Eth = 3 * n * R_b * T * debye(theta / T)
    Eth0 = 3 * n * R_b * T * debye(theta / T0)
    E = gamma / V * (Eth - Eth0)
    P = 1.5 * K0 * (x ** (7.0 / 3.0) - x ** (5.0 / 3.0)) * (1 + 0.75 * (KP0 - 4) * (x ** (2.0 / 3.0) - 1)) + E

    return P


def Holz(V, T, Z, T0, V0, M, K0, KP0, theta0, gamma0, gammaInf, beta, a0, m, g):
    """ from Lena Noack """
    PFG0 = 0.10036e9 * (Z / V0) ** (5.0 / 3.0)  # 1003.6e9 -> 0.10036e9 due to V0*1000^5/3=value*10^5
    c0 = -np.log(3.0 * K0 / PFG0)
    c2 = 3.0 / 2.0 * (KP0 - 3.0) - c0
    zeta = (V / V0) ** (1.0 / 3.0)  # zetax(p,T,errorCode) ! zeta=(V/V0)**(1./3.)
    zeta3 = zeta ** 3

    theta = theta0 * zeta3 ** (-gammaInf) * np.exp((1 - zeta3 ** beta) * (gamma0 - gammaInf) / beta)
    gamma = (gamma0 - gammaInf) * zeta3 ** beta + gammaInf

    E = gamma / V * (3 * R_b * (0.5 * theta + theta / (np.exp(theta / T) - 1)) - 3 * R_b * (
                0.5 * theta + theta / (np.exp(theta / T0) - 1)))

    P = 3.0 * K0 * zeta ** (-5) * (1 - zeta) * np.exp(c0 * (1 - zeta)) * (1 + c2 * (zeta - zeta ** 2)) + E
    P = P + 1.5 * R_b / V * m * a0 * zeta ** (3.0 * m) * T ** 2

    return P


def Holz_fail(V, T, Z, T0, V0, M, KT0, KP0, theta0, gamma0, gammaInf, beta, a0, m, g, c0=None, c2=None, P0=0, V0_Bou=None):
    if c0 is None and c2 is None:
        PFG0 = 0.10036e9 * (Z / V0) ** (5.0 / 3.0)  # 1003.6e9 -> 0.10036e9 due to V0*1000^5/3=value*10^5
        c0 = -np.log(3.0 * KT0 / PFG0)
        c2 = 3.0 / 2.0 * (KP0 - 3.0) - c0

    x = V / V0
    x_corr = x * V0 / V0_Bou  # thermal part corrected to match Bouchet V0

    theta = theta0 * x ** (-gammaInf) * np.exp((1 - x_corr ** beta) * (gamma0 - gammaInf) / beta)
    gamma = (gamma0 - gammaInf) * x_corr ** beta + gammaInf

    P_iso = P0 + 3.0 * KT0 * x ** (-5/3) * (1 - x ** (1/3)) * np.exp(c0 * (1 - x ** (1/3))) * (1 + c2 * x ** (1/3) * (1 - x ** (1/3)))
    P_harm = 3 * R_b * gamma / V * (theta / 2 + theta / (np.exp(theta / T) - 1))
    P_ae = 3 * R_b / (2 * V) * m * a0 * x_corr ** m * T ** 2
    P = P_iso + P_harm + P_ae

    return P



def get_V(P, T, n, T0, V0, M, K0, KP0, G0, GP0, theta0, gamma0, gammaInf, q0, etaS0, a0, m, g, EOS, **kwargs):
    """ from Lena Noack """
    # get V that fits to P and T

    tau = 0.01
    maxIter = 100

    V_i = V0
    it = 1
    dV = 0.5 * V0
    forward = -1

    if EOS == 3:
        P_i = BM3(V_i, T, n, T0, V0, M, K0, KP0, G0, GP0, theta0, gamma0, q0, etaS0)
    elif EOS == 2:
        P_i = BM3(V_i, T, n, T0, V0, M, K0, 4.0, G0, GP0, theta0, gamma0, q0, etaS0)
    else:
        P_i = Holz(V_i, T, n, T0, V0, M, K0, KP0, theta0, gamma0, gammaInf, q0, a0, m, g, **kwargs)

    while (abs(P_i - P) > tau) and (it < maxIter):
        it = it + 1
        if P_i > P:
            V_i = V_i + dV
        else:
            V_i = V_i - dV
        dV = 0.8 * dV

        if EOS == 3:
            P_i = BM3(V_i, T, n, T0, V0, M, K0, KP0, G0, GP0, theta0, gamma0, q0, etaS0)
        elif EOS == 2:
            P_i = BM3(V_i, T, n, T0, V0, M, K0, 4.0, G0, GP0, theta0, gamma0, q0, etaS0)
        else:
            P_i = Holz(V_i, T, n, T0, V0, M, K0, KP0, theta0, gamma0, gammaInf, q0, a0, m, g, **kwargs)

    V = V_i
    return V


def EOS_all(P, T, material):
    """ from Lena Noack, P in GPa, T in K """

    # Computes thermodynamic properties molar volume V, density rho, thermal expansivity alpha,
    # and heat capacity Cp for a given pressure, temperature and composition (material)
    # material = 1 -> Forsterite
    # material = 2 -> Mg-Perovskite + Periclase
    # material = 3 -> Mg-Post-Perovskite + Periclase
    # material = 4 -> Iron hcp phase
    #
    # volume V(p,T) is obtained via get_V(...,EOS)
    # EOS = 1: Holzapfel EOS
    # EOS = 2: Second-order Birch-Murnaghan EOS
    # EOS = 3: Third-order Birch-Murnaghan EOS

    # set EOS parameters
    V = 0
    rho = 0
    alpha = 0
    Cp = 0
    if material == 1:  # V0 is in mm^3/mol
        # Forsterite
        n = 7
        T0 = 298
        V0 = 43600
        M = 140.693
        K0 = 128
        KP0 = 4.2
        G0 = 82
        GP0 = 1.5
        theta0 = 809
        gamma0 = 0.99
        gammaInf = 0
        q0 = 2.1
        etaS0 = 2.3
        a0 = 0
        m = 0
        g = 0
        EOS = 3  # 3rd-order Birch-Murnaghan
        V = get_V(P, T, n, T0, V0, M, K0, KP0, G0, GP0, theta0, gamma0, gammaInf, q0, etaS0, a0, m, g, EOS)
        x = V0 / V
        f = 0.5 * (x ** (2 / 3) - 1.0)
        theta = theta0 * np.sqrt(1.0 + 6.0 * f * gamma0 + 3.0 * f ** 2 * gamma0 * (-2.0 + 6.0 * gamma0 - 3.0 * q0))
        gamma = (1.0 + 2.0 * f) * gamma0 * (1.0 + f * (-2.0 + 6.0 * gamma0 - 3.0 * q0)) * (theta0 / theta) ** 2

        Cv = 3 * n * R_b * (4 * debye(theta / T) - (theta / T) * 3 / (np.exp(theta / T) - 1))
        Cv0 = 3 * n * R_b * (4 * debye(theta / T0) - (theta / T0) * 3 / (np.exp(theta / T0) - 1))
        q = (-2 * gamma + 6. * gamma ** 2 + (1 + 2 * f) ** 2 * gamma0 * (2 - 6 * gamma0 + 3 * q0) * (
                theta0 / theta) ** 2) / (3 * gamma)
        KT = (1 + 2 * f) ** (5 / 2) * K0 * (1 + f * (-5 + 3 * KP0) + 27 / 2 * f ** 2 * (-4 + KP0))
        KT = KT + (-gamma ** 2 / V * (Cv * T - Cv0 * T0) + gamma / V * (1 - q + gamma) * (
                3 * n * R_b * T * debye(theta / T) - 3 * n * R_b * T0 * debye(theta / T0)))
        alpha = gamma * Cv / (KT * V)
        Cp = Cv * (1 + alpha * gamma * T) / (1e-3 * M)  # divide by mol mass
        rho = 1e6 * M / V
    elif material == 2:
        # Mg-Perovskite
        n = 5
        T0 = 298
        V0 = 24450
        M = 100.389
        K0 = 251
        KP0 = 4.1
        G0 = 173
        GP0 = 1.7
        theta0 = 905
        gamma0 = 1.5
        gammaInf = 0
        q0 = 1.1
        etaS0 = 2.6
        a0 = 0
        m = 0
        g = 0
        EOS = 3  # 3rd-order Birch-Murnaghan
        V = get_V(P, T, n, T0, V0, M, K0, KP0, G0, GP0, theta0, gamma0, gammaInf, q0, etaS0, a0, m, g, EOS)

        x = V0 / V
        f = 0.5 * (x ** (2 / 3) - 1.0)
        theta = theta0 * np.sqrt(1.0 + 6.0 * f * gamma0 + 3.0 * f ** 2 * gamma0 * (-2.0 + 6.0 * gamma0 - 3.0 * q0))
        gamma = (1.0 + 2.0 * f) * gamma0 * (1.0 + f * (-2.0 + 6.0 * gamma0 - 3.0 * q0)) * (theta0 / theta) ** 2

        Cv = 3 * n * R_b * (4 * debye(theta / T) - (theta / T) * 3 / (np.exp(theta / T) - 1))
        Cv0 = 3 * n * R_b * (4 * debye(theta / T0) - (theta / T0) * 3 / (np.exp(theta / T0) - 1))
        q = (-2 * gamma + 6 * gamma ** 2 + (1 + 2 * f) ** 2 * gamma0 * (2 - 6 * gamma0 + 3 * q0) * (
                theta0 / theta) ** 2) / (3 * gamma)
        KT = (1 + 2 * f) ** (5 / 2) * K0 * (1 + f * (-5 + 3 * KP0) + 27 / 2 * f ** 2 * (-4 + KP0))
        KT = KT + (-gamma ** 2 / V * (Cv * T - Cv0 * T0) + gamma / V * (1 - q + gamma) * (
                3 * n * R_b * T * debye(theta / T) - 3 * n * R_b * T0 * debye(theta / T0)))
        alpha = gamma * Cv / (KT * V)
        Cp = Cv * (1 + alpha * gamma * T) / (1e-3 * M)  # divide by mol mass
        rho = 1e6 * M / V

        Pv_Cp = Cp
        Pv_alpha = alpha
        Pv_rho = rho

        # Periclase
        n = 2
        T0 = 298
        V0 = 11240
        M = 40.3044
        K0 = 161
        KP0 = 3.8
        G0 = 131
        GP0 = 2.1
        theta0 = 767
        gamma0 = 1.36
        gammaInf = 0
        q0 = 1.7
        etaS0 = 2.8
        a0 = 0
        m = 0
        g = 0
        EOS = 3  # 3rd-order Birch-Murnaghan
        V = get_V(P, T, n, T0, V0, M, K0, KP0, G0, GP0, theta0, gamma0, gammaInf, q0, etaS0, a0, m, g, EOS)

        x = V0 / V
        f = 0.5 * (x ** (2 / 3) - 1.0)
        theta = theta0 * np.sqrt(1.0 + 6.0 * f * gamma0 + 3.0 * f ** 2 * gamma0 * (-2.0 + 6.0 * gamma0 - 3.0 * q0))
        gamma = (1.0 + 2.0 * f) * gamma0 * (1.0 + f * (-2.0 + 6.0 * gamma0 - 3.0 * q0)) * (theta0 / theta) ** 2

        Cv = 3 * n * R_b * (4 * debye(theta / T) - (theta / T) * 3 / (math.exp(theta / T) - 1))
        Cv0 = 3 * n * R_b * (4 * debye(theta / T0) - (theta / T0) * 3 / (math.exp(theta / T0) - 1))
        q = (-2 * gamma + 6 * gamma ** 2 + (1 + 2 * f) ** 2 * gamma0 * (2 - 6 * gamma0 + 3 * q0) * (
                theta0 / theta) ** 2) / (3 * gamma)
        KT = (1 + 2 * f) ** (5 / 2) * K0 * (1 + f * (-5 + 3 * KP0) + 27 / 2 * f ** 2 * (-4 + KP0))
        KT = KT + (-gamma ** 2 / V * (Cv * T - Cv0 * T0) + gamma / V * (1 - q + gamma) * (
                3 * n * R_b * T * debye(theta / T) - 3 * n * R_b * T0 * debye(theta / T0)))
        alpha = gamma * Cv / (KT * V)
        Cp = Cv * (1 + alpha * gamma * T) / (1e-3 * M)  # divide by mol mass
        rho = 1e6 * M / V

        Pe_Cp = Cp
        Pe_alpha = alpha
        Pe_rho = rho

        mf_pv = 0.7135
        mf_pe = 0.2865

        rho = (mf_pv / Pv_rho + mf_pe / Pe_rho) ** (-1)
        Cp = mf_pv * Pv_Cp + mf_pe * Pe_Cp
        alpha = mf_pv * Pv_alpha * rho / Pv_rho + mf_pe * Pe_alpha * rho / Pe_rho
    elif material == 3:
        # Mg-Post-Perovskite
        n = 5
        T0 = 298
        V0 = 24420
        M = 100.389
        K0 = 231
        KP0 = 4.0
        G0 = 150
        GP0 = 2.0
        theta0 = 855
        gamma0 = 1.89
        gammaInf = 0
        q0 = 1.1
        etaS0 = 1.2
        a0 = 0
        m = 0
        g = 0
        EOS = 3  # 3rd-order Birch-Murnaghan
        V = get_V(P, T, n, T0, V0, M, K0, KP0, G0, GP0, theta0, gamma0, gammaInf, q0, etaS0, a0, m, g, EOS)

        x = V0 / V
        f = 0.5 * (x ** (2 / 3) - 1.0)
        theta = theta0 * math.sqrt(1.0 + 6.0 * f * gamma0 + 3.0 * f ** 2 * gamma0 * (-2.0 + 6.0 * gamma0 - 3.0 * q0))
        gamma = (1.0 + 2.0 * f) * gamma0 * (1.0 + f * (-2.0 + 6.0 * gamma0 - 3.0 * q0)) * (theta0 / theta) ** 2

        Cv = 3 * n * R_b * (4 * debye(theta / T) - (theta / T) * 3 / (math.exp(theta / T) - 1))
        Cv0 = 3 * n * R_b * (4 * debye(theta / T0) - (theta / T0) * 3 / (math.exp(theta / T0) - 1))
        q = (-2 * gamma + 6 * gamma ** 2 + (1 + 2 * f) ** 2 * gamma0 * (2 - 6 * gamma0 + 3 * q0) * (
                theta0 / theta) ** 2) / (3 * gamma)
        KT = (1 + 2 * f) ** (5 / 2) * K0 * (1 + f * (-5 + 3 * KP0) + 27 / 2 * f ** 2 * (-4 + KP0))
        KT = KT + (-gamma ** 2 / V * (Cv * T - Cv0 * T0) + gamma / V * (1 - q + gamma) * (
                3 * n * R_b * T * debye(theta / T) - 3 * n * R_b * T0 * debye(theta / T0)))
        alpha = gamma * Cv / (KT * V)
        Cp = Cv * (1 + alpha * gamma * T) / (1e-3 * M)  # divide by mol mass
        rho = 1e6 * M / V

        Ppv_Cp = Cp
        Ppv_alpha = alpha
        Ppv_rho = rho

        # Periclase
        n = 2
        T0 = 298
        V0 = 11240
        M = 40.3044
        K0 = 161
        KP0 = 3.8
        G0 = 131
        GP0 = 2.1
        theta0 = 767
        gamma0 = 1.36
        gammaInf = 0
        q0 = 1.7
        etaS0 = 2.8
        a0 = 0
        m = 0
        g = 0
        EOS = 3  # 3rd-order Birch-Murnaghan
        V = get_V(P, T, n, T0, V0, M, K0, KP0, G0, GP0, theta0, gamma0, gammaInf, q0, etaS0, a0, m, g, EOS)

        x = V0 / V
        f = 0.5 * (x ** (2 / 3) - 1.0)
        theta = theta0 * math.sqrt(1.0 + 6.0 * f * gamma0 + 3.0 * f ** 2 * gamma0 * (-2.0 + 6.0 * gamma0 - 3.0 * q0))
        gamma = (1.0 + 2.0 * f) * gamma0 * (1.0 + f * (-2.0 + 6.0 * gamma0 - 3.0 * q0)) * (theta0 / theta) ** 2

        Cv = 3 * n * R_b * (4 * debye(theta / T) - (theta / T) * 3 / (math.exp(theta / T) - 1))
        Cv0 = 3 * n * R_b * (4 * debye(theta / T0) - (theta / T0) * 3 / (math.exp(theta / T0) - 1))
        q = (-2 * gamma + 6 * gamma ** 2 + (1 + 2 * f) ** 2 * gamma0 * (2 - 6 * gamma0 + 3 * q0) * (
                theta0 / theta) ** 2) / (3 * gamma)
        KT = (1 + 2 * f) ** (5 / 2) * K0 * (1 + f * (-5 + 3 * KP0) + 27 / 2 * f ** 2 * (-4 + KP0))
        KT = KT + (-gamma ** 2 / V * (Cv * T - Cv0 * T0) + gamma / V * (1 - q + gamma) * (
                3 * n * R_b * T * debye(theta / T) - 3 * n * R_b * T0 * debye(theta / T0)))
        alpha = gamma * Cv / (KT * V)
        Cp = Cv * (1 + alpha * gamma * T) / (1e-3 * M)  # divide by mol mass
        rho = 1e6 * M / V

        Pe_Cp = Cp
        Pe_alpha = alpha
        Pe_rho = rho

        mf_pv = 0.7135
        mf_pe = 0.2865

        rho = (mf_pv / Ppv_rho + mf_pe / Pe_rho) ** (-1)
        Cp = mf_pv * Ppv_Cp + mf_pe * Pe_Cp
        alpha = mf_pv * Ppv_alpha * rho / Ppv_rho + mf_pe * Pe_alpha * rho / Pe_rho
    elif material == 4:  # hcp-iron
        # n = 26  # Z in EOS
        # T0 = 300.0
        # q0 = 0.826  # beta in EOS
        # a0 = 0.0002121
        # M = 55.845e6  # molar mass
        # m = 1.891
        # g = 1.339
        # G0 = 0
        # GP0 = 0
        # etaS0 = 0
        # gamma0 = 1.408
        # gammaInf = 0.827
        # theta0 = 44.574
        # V0_Bou = 6290.0
        # V0 = 4285.75

        # V0 = 6290.0
        # KP0 = 4.719
        # K0 = 253.844

        # P0 = 234.4
        # c0 = 3.19
        # c2 = 2.40
        # KT0 = 1145.7

        if P > 10e3:
            raise Exception('EXTRAPOLATION ERROR: Hakim+ Holzapfel EoS for Fe')
        # V = get_V(P, T, n, T0, V0, M, KT0, None, G0, GP0, theta0, gamma0, gammaInf, q0, etaS0, a0, m, g, EOS,
        #           P0=P0, c0=c0, c2=c2, V0_Bou=V0_Bou)

        n = 26  # Z in EOS
        T0 = 300.0
        V0 = 6290.0  # 6.290 cm3/mol
        M = 55.845e6
        K0 = 253.844
        KP0 = 4.719
        theta0 = 44.574
        gamma0 = 1.408
        gammaInf = 0.827
        q0 = 0.826  # beta in EOS
        a0 = 0.0002121
        m = 1.891
        g = 1.339
        G0 = 0
        GP0 = 0
        etaS0 = 0
        EOS = 1  # Holzapfel EOS
        V = get_V(P, T, n, T0, V0, M, K0, KP0, G0, GP0, theta0, gamma0, gammaInf, q0, etaS0, a0, m, g, EOS)

        Z = n
        beta = q0

        PFG0 = 0.10036e9 * (Z / V0) ** (5.0 / 3.0)  # e9 or e5?
        c0 = -math.log(3.0 * K0 / PFG0)
        c2 = 3.0 / 2.0 * (KP0 - 3.0) - c0
        zeta = (V / V0) ** (1.0 / 3.0)  # zetax(p,T,errorCode) ! zeta=(V/V0)**(1./3.)
        zeta3 = zeta ** 3  # V/V0
        rho = M / V
        gamma = (gamma0 - gammaInf) * zeta3 ** beta + gammaInf
        theta = theta0 * zeta3 ** (-gammaInf) * math.exp((1 - zeta3 ** beta) * (gamma0 - gammaInf) / beta)
        Cv = 3.0 * R_b * theta ** 2 * math.exp(theta / T) / ((-1.0 + math.exp(theta / T)) ** 2 * T ** 2)
        Cv0 = 3.0 * R_b * theta ** 2 * math.exp(theta / T0) / ((-1.0 + math.exp(theta / T0)) ** 2 * T0 ** 2)
        Eth = 3.0 * R_b * (theta / 2.0 + theta / (math.exp(theta / T) - 1.0))
        Eth0 = 3.0 * R_b * (theta / 2.0 + theta / (math.exp(theta / T0) - 1.0))
        dEeadx = 1.5 * R_b * m ** 2 * a0 * zeta3 ** (m - 1) * T ** 2
        dEeadx0 = 1.5 * R_b * m ** 2 * a0 * zeta3 ** (m - 1) * T0 ** 2
        Eea = 1.5 * R_b * m * a0 * zeta3 ** m * T ** 2
        Eea0 = 1.5 * R_b * m * a0 * zeta3 ** m * T0 ** 2
        KT = (math.exp(c0 - c0 * zeta) * K0 * (5.0 + zeta * (
                -4.0 + 2.0 * c2 * (-2.0 + zeta) * (-1.0 + zeta) + c0 * (-1.0 + zeta) * (
                -1.0 + c2 * (-1.0 + zeta) * zeta)))) / zeta ** 5
        KT = KT + (-(gamma ** 2 * (T * Cv - T0 * Cv0)) + (gamma * (1.0 - beta + gamma) + beta * gammaInf) * (
                Eth - Eth0)) / V
        KT = KT - (dEeadx - dEeadx) / V0 + (Eea - Eea0) / V
        Cv = Cv + 3.0 * R_b * m * a0 * zeta3 ** m * T
        alpha = gamma * Cv / (KT * V)
        Cp = Cv * (1.0 + alpha * gamma * T) / (M * 1.0e-9)  # from J/mol K to J/kg K: division by mol mass

    return [V, rho, alpha, Cp]


def g_profile(n, radius, density):
    # gravity has boundary condition 0 in center (and not at surface)
    # analogously, a surface gravity can be determind from M and
    # guessed R, and gravity can be interpolated from surface downwards
    # Problem with that approach: neg. gravity values in center possible
    gravity = np.zeros(n)
    for i in range(1, n):
        dr = radius[i] - radius[i - 1]
        gravity[i] = gravity[i - 1] + dr * (
                4 * math.pi * G * density[i - 1] - 2 * gravity[i - 1] / radius[i - 1])
    return gravity


def pt_profile(n, radius, density, gravity, alpha, cp, psurf, tsurf, i_cmb=None, deltaT_cmb=0, extrapolate_nan=False):
    """ input psurf in bar, pressure and temperature are interpolated from surface downwards """
    pressure = np.zeros(n)
    temperature = np.zeros(n)
    pressure[-1] = psurf * 1e5  # surface pressure in Pa
    temperature[-1] = tsurf  # potential surface temperature, not real surface temperature
    # last_good_lapse_rate = alpha[-1] / cp[-1] * gravity[-1]  # temp
    # flag = True
    for i in range(2, n + 1):  # M: 1:n-1; n-i from n-1...1; P: from n-2...0 -> i from 2...n
        dr = radius[n - i + 1] - radius[n - i]
        pressure[n - i] = pressure[n - i + 1] + dr * gravity[n - i] * density[n - i]
        lapse_rate = alpha[n - i] / cp[n - i] * gravity[n - i]
        # if lapse_rate == 0 or np.isnan(lapse_rate):
        #     print('n - i', n - i, 'pressure[n - i]', pressure[n - i]*1e-9, 'GPa', 'alpha', alpha[n - i], 'cp', cp[n - i], 'g', gravity[n - i])
        #     if extrapolate_nan:
        #             lapse_rate = last_good_lapse_rate
        #             if flag:
        #                 print('bad lapse rate from p =', pressure[n - i]*1e-9, 'GPa')
        #                 flag = False
        #     else:
        #         raise ZeroDivisionError()
        # last_good_lapse_rate = lapse_rate
        ## ^^tried to implement constant extrapolation of adiabatic gradient but still need something for bad density to get pressure
        temperature[n - i] = temperature[n - i + 1] + dr * lapse_rate * temperature[n - i + 1]
        if n - i == i_cmb:
            # add temperature jump across cmb (discontinuity)
            temperature[n - i] = temperature[n - i] + deltaT_cmb
    return pressure, temperature  # Pa, K

def pt_profile_dynamic(n, radius, density, gravity, alpha, cp, psurf, tsurf, i_cmb=None, deltaT_cmb=0):
    # input psurf in bar
    # pressure and temperature are interpolated from surface downwards with possibly different array lengths
    pressure = np.zeros(n)
    temperature = np.zeros(n)
    pressure[-1] = psurf * 1e5  # surface pressure in Pa
    temperature[-1] = tsurf  # potential surface temperature, not real surface temperature
    for i in range(2, n + 1):  # M: 1:n-1; n-i from n-1...1; P: from n-2...0 -> i from 2...n
        dr = radius[n - i + 1] - radius[n - i]
        pressure[n - i] = pressure[n - i + 1] + dr * gravity[n - i] * density[n - i]
        temperature[n - i] = temperature[n - i + 1] + dr * alpha[n - i] / cp[n - i] * gravity[n - i] * temperature[n - i + 1]
        if n - i == i_cmb:
            # add temperature jump across cmb (discontinuity)
            temperature[n - i] = temperature[n - i] + deltaT_cmb
        if math.isinf(temperature[n-i]):
            print('i', i)
            print('alpha', alpha[n - i])
            print('cp', cp[n - i])
            print('gravity', gravity[n - i])
            raise ZeroDivisionError
    return pressure, temperature  # Pa, K


def core_planet_radii(CMF, Mp, rho_c, rho_m):  # from bulk i.e. average values
    if CMF > 0:
        Rc = (3 * CMF * Mp / (4 * math.pi * rho_c)) ** (1 / 3)
    else:
        Rc = 0
    Rp = (Rc ** 3 + 3 * (1 - CMF) * Mp / (4 * math.pi * rho_m)) ** (1 / 3)
    return Rc, Rp


def EOS_lm(P, T, phases, proportions):
    """ P in GPa, T in K call burnman stuff to use their built in eos for lower mantle
    # proportions in wt % (same as perple_x compositional output """
    import burnman_eos
    alphas = []
    cps = []
    rhos = []
    for (wt, phase) in zip(proportions, phases):
        [V, rho, alpha, Cp] = burnman_eos.get_properties(P*1e9, T, phase)
        rhos.append(rho)
        alphas.append(alpha)  # need total density first
        cps.append(Cp)

    inv_rho_tot = 0
    cp_tot = 0
    alpha_tot = 0
    for (wt, r, c) in zip(proportions, rhos, cps):
        inv_rho_tot += wt * 1e-2 / r
        cp_tot += wt * 1e-2 * c

    rho_tot = inv_rho_tot**-1

    for (wt, r, a) in zip(proportions, rhos, alphas):
        alpha_tot += wt * 1e-2 * a * rho_tot / r

    return None, rho_tot, alpha_tot, cp_tot


# import matplotlib.pyplot as plt
# P = np.linspace(0.1, 300)
# T = 300*np.ones_like(P)
# V, rho, alpha, cp = [], [], [], []
# for pp, TT in zip(P, T):
#     V1, rho1, alpha1, cp1 = EOS_all(pp, TT, material=4)
#     rho.append(rho1)
#     alpha.append(alpha1)
#     cp.append(cp1)
#     V.append(V1)
# # _, rho, alpha, cp = EOS_lm(P, T, ['ppv', 'per'], [80, 20])
# fig, axes = plt.subplots(1, 4)
# axes[0].plot(P, rho)
# axes[1].plot(P, alpha)
# axes[2].plot(P, cp)
# axes[3].plot(V, P)
# axes[0].axvline(234)
# axes[0].set_ylabel('density')
# axes[1].set_ylabel('alpha')
# axes[2].set_ylabel('cp')
# plt.show()
#
