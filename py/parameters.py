###### PHYSICAL CONSTANTS ######
M_E = 5.972e24  # earth mass in kg
R_E = 6371e3  # earth radius in m
rho_E = 5513.258738589093  # earth bulk density in kg/m3 using above
TO = 1.335e21  # earth ocean mass in kg
G = 6.67408e-11
years2sec = 31557600
sec2Gyr = 1e-9 / years2sec
AU2m = 1.5e11
R_b = 8.3144598  # universal gas constant in J mol −1 K −1
sb = 5.67e-8  # Stefan Boltzmann constant in W m^-2 K^-4

# molar masses
M_O = 15.999
M_Si = 28.0855
M_Mg = 24.305
M_Ca = 40.078
M_Al = 26.981539
M_Fe = 55.845
M_Na = 22.989769

# molar masses of oxides
M_SiO2 = M_Si + 2 * M_O
M_CaO = M_Ca + M_O
M_MgO = M_Mg + M_O
M_Al2O3 = 2 * M_Al + 3 * M_O
M_FeO = M_Fe + M_O
M_Na2O = 2 * M_Na + M_O

# solar composition from Lodders 09 in log(n/H)
ca_sol = 6.33 - 12
al_sol = 6.47 - 12
fe_sol = 7.45 - 12
si_sol = 7.52 - 12
mg_sol = 7.54 - 12
na_sol = 6.30 - 12

# print(10 ** mg_sol / 10 ** si_sol)   # sun = 1.04
