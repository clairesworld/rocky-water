import matplotlib.pyplot as plt
import numpy as np

M_E = 5.972e24  # earth mass in kg
R_E = 6371e3  # earth radius in m
rho_av = 5513.258738589093  # earth bulk density in kg/m3 using above
rho_m = 4500
G = 6.67408e-11
P0 = 25e9
Psfc = 1000e5

def um_mass_scaling(M_p=None, R_p=None, plot=False):
    if R_p is None and M_p is not None:
        R_p = (3 * M_p / (4 * np.pi * rho_av)) ** (1/3)
    elif M_p is None and R_p is not None:
        M_p = rho_av * 4/3 * np.pi * R_p ** 3
    else:
        raise NotImplementedError('um_mass_scaling() needs M_p or R_p')

    R0 = (Psfc - P0 + (4/3 * np.pi * G * rho_av * rho_m * R_p**2)) / (4/3 * np.pi * G * rho_av * rho_m * R_p)
    M_um = 4/3 * rho_m * np.pi * (R_p ** 3 - R0 ** 3)

    if plot:
        plt.plot(R_p / R_E, M_um)
        plt.xlabel('R_p (R_E)')
        plt.ylabel('upper mantle mass (kg)')
        plt.show()

    return M_um ## kg

