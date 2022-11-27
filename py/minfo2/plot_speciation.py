# plot molar gas fraction dependence on fo2

import numpy as np
import matplotlib.pyplot as plt
from py.useful_and_bespoke import colorize, dark_background
from oxygen_fugacity_plots import figpath
from py.minfo2.perplexfugacitydata import read_qfm_os
from py.perplexdata import perplex_path_default

# parameterised IW
initial_H2O_mole_fraction = 0.5
initial_CO2_mole_fraction = 0.5
T = 1800  # K
P = 1  # bar
a = -6
b = 4
buffer = np.linspace(a, b, 100)
LogfO2_IW = (-27714 / T) + 6.899 + ((0.05 * (P - 1) / T) + buffer)
O2_fugacity_IW = 10 ** LogfO2_IW

# oli's QFM
LogfO2_QFM = read_qfm_os(T=T, P=P, perplex_path=perplex_path_default, fin='data_tables/fmqNl_fo2_oli.dat', verbose=False) + buffer
O2_fugacity_QFM = 10 ** LogfO2_QFM

Gf_CH4 = -77437 + 22.5098 * T * np.log10(T) + 29.5967 * T
g_co = -214104 + 25.2183 * T * np.log10(T) - 262.1545 * T
Gf_CO = g_co / 2
Gf_CO2 = -392647 + 4.5855 * T * np.log10(T) - 16.9762 * T
g_h2o = -483095 + 25.3687 * T * np.log10(T) + 21.9563 * T
Gf_H2O = g_h2o / 2

# 1) H2O
K1 = np.exp(((-Gf_H2O) * 2) / (8.314 * T))

# oxygen fugacity at IW buffer
ratio_H2O_H2_IW = np.sqrt((O2_fugacity_IW) * K1)
H2_mole_fraction_IW = (initial_H2O_mole_fraction / (
            ratio_H2O_H2_IW + 1))  # (initial H2O mole fraction is from the solubility model)
H2O_mole_fraction_IW = initial_H2O_mole_fraction - H2_mole_fraction_IW

# at QFM
ratio_H2O_H2_QFM = np.sqrt((O2_fugacity_QFM) * K1)
H2_mole_fraction_QFM = (initial_H2O_mole_fraction / (
            ratio_H2O_H2_QFM + 1))  # (initial H2O mole fraction is from the solubility model)
# H2_wt_percent_IW=H2_mole_fraction_IW*100
H2O_mole_fraction_QFM = initial_H2O_mole_fraction - H2_mole_fraction_QFM

# 2)   CO + 1/2O2 = CO2
deltaG_reaction_CO_CO2 = Gf_CO2 - Gf_CO
K2 = np.exp(-deltaG_reaction_CO_CO2 / (8.314 * T))
# oxygen fugacity at IW buffer
ratio_CO2_CO_IW = ((O2_fugacity_IW) ** 0.5) * K2
CO_mole_fraction_IW = (initial_CO2_mole_fraction / (ratio_CO2_CO_IW + 1))
CO2_mole_fraction_IW = initial_CO2_mole_fraction - CO_mole_fraction_IW

# at QFM
ratio_CO2_CO_QFM= ((O2_fugacity_QFM) ** 0.5) * K2
CO_mole_fraction_QFM = (initial_CO2_mole_fraction / (ratio_CO2_CO_QFM + 1))
CO2_mole_fraction_QFM = initial_CO2_mole_fraction - CO_mole_fraction_QFM

# plot
# plt.style.use('dark_background')
fig, ax = plt.subplots(1, 1, figsize=(11, 1))
# fig, ax = plt.subplots(1, 1, figsize=(6, 3))
lw = 3
fontsize = 18  # 22
colors = [colorize(np.arange(4), cmap='Spectral_r')[0][ii] for ii in range(4)]
colors[0] = 'xkcd:pale blue'
ls = ['-', '--', '-.', ':']
# gases = [H2O_mole_fraction_IW, H2_mole_fraction_IW, CO2_mole_fraction_IW, CO_mole_fraction_IW]
# gases = [H2O_mole_fraction_QFM, H2_mole_fraction_QFM, CO2_mole_fraction_QFM, CO_mole_fraction_QFM]
gases = [H2O_mole_fraction_QFM, H2_mole_fraction_QFM]
labels = [r'H$_2$O', r'H$_2$', r'CO$_2$', r'CO']

for ii, gas in enumerate(gases):
    ax.plot(buffer, gas, lw=lw, label=labels[ii], color=colors[ii], ls=ls[ii])

ax.set_xlabel(r'log$f{\rmO}_2$ ($\Delta$QFM)', fontsize=fontsize + 4, labelpad=10)
ax.set_ylabel('Volatile mole fraction', fontsize=fontsize + 4)
ax.set_title(r'1 bar, 1800 K', fontsize=fontsize + 4)
ax.set_xlim(buffer.min(), buffer.max())
ax.tick_params(axis='both', which='major', labelsize=fontsize - 2)
ax.legend(frameon=False, fontsize=fontsize, bbox_to_anchor=(1.01, 1), loc='upper left')
ax.set_xlim(-4, 3)
# ax.set_xlim(-6, 1)

fig, ax = dark_background(fig, ax)

fig.savefig(figpath + 'speciation_delta_qfm2.png', bbox_inches='tight')
plt.show()
