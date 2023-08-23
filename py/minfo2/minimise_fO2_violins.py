import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import FormatStrFormatter
from py.useful_and_bespoke import dark_background
from oxygen_fugacity_plots import figpath
from matplotlib import rc
import datetime

today = str(datetime.date.today())
# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

datapath = '/home/claire/Works/min-fo2/sean_fastchem/'
datapath_sj = './fO2_out/claire_project_violins_4GPa/'

pad_size = 12
label_size = 24
tick_size = 18  # cg
major_length = 4
minor_length = 2
tick_width = 2
axes_label = 18
leg_font = 16
spec_font = 20

fig, ax = plt.subplots(figsize=(12, 10))
gs1 = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs1[0])
ax2 = plt.subplot(gs1[1])
ax3 = plt.subplot(gs1[2])
ax4 = plt.subplot(gs1[3])

for axis in [ax1, ax2, ax3, ax4]:

    if axis == ax1:
        chem_in = 'table_FMQ_solar_COSH_chem.dat'
        norm_in = 'table_FMQ_solar_COSH_norm.dat'
    if axis == ax2:
        chem_in = 'table_FMQ_Venus_atm_COSH_chem.dat'
        norm_in = 'table_FMQ_Venus_atm_COSH_norm.dat'
    if axis == ax3:
        chem_in = 'table_IW-2_solar_COSH_chem.dat'
        norm_in = 'table_IW-2_solar_COSH_norm.dat'
    if axis == ax4:
        chem_in = 'table_IW-2_Venus_atm_COSH_chem.dat'
        norm_in = 'table_IW-2_Venus_atm_COSH_norm.dat'

    df = pd.read_csv(datapath + chem_in, delimiter='\t')
    df_norm = pd.read_csv(datapath + norm_in, delimiter='\t')
    D = []
    # D2=[]
    species_list = ['CO', 'CO2', 'CH4', 'H2O', 'H2S', 'SO2']
    colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:pink', 'tab:brown', 'tab:olive',
               'tab:cyan', 'tab:grey']
    p = 0
    pos = [3, 6, 9, 12, 15, 18]  # [4,8,12,16,20,24]#[3,6,9,12,15,18]#[2,4,6,8,10]
    for species in species_list:

        norm = df_norm[species].values
        print(species, norm)
        observable = 1e-6 / norm
        logobservable = np.log10(observable)
        print(observable)
        data = df[species].values / norm
        logdata = np.log10(df[species].values / norm)
        random_scatter = np.ones_like(logdata)
        for i in range(len(logdata)):
            random_scatter[i] = np.random.rand()
        # axis.scatter(random_scatter+pos[p]-2.5, logdata, s=1, marker='x', color='grey', alpha=0.5)
        xline = np.linspace(pos[p] - 1, pos[p] + 1, 2)
        D.append(logdata)

        # axis.plot(xline,np.ones_like(xline)*logobservable, c='r')
        axis.fill_between([pos[p] - 1.5, pos[p] + 1.5], y1=logobservable, y2=-4, color='tab:grey', edgecolor=None,
                          alpha=0.2, zorder=2)

        # test for plotting unobservable region
        #    unobservable = [x for x in logdata if x <= logobservable]
        #    if len(unobservable)==0: unobservable =[0]
        #    D2.append(unobservable)

        p += 1

    vp = axis.violinplot(D, positions=pos, widths=2, showmeans=False, showextrema=False,
                         showmedians=False)  # , points=1000)
    # vp2= axis.violinplot(D2,positions = pos, widths=2, showmeans=False, showextrema=False, showmedians=False)#, points=1000)

    col = 0
    for body in vp['bodies']:
        body.set_alpha(0.6)
        body.set_facecolor(colours[col])
        body.set_edgecolor(colours[col])
        col += 1
    # for body in vp2['bodies']:
    #    body.set_alpha(0.9)
    #    body.set_facecolor('tab:grey')

    # [spine.set_linewidth(tick_width) for spine in axis.spines.values()]
    axis.tick_params(axis='both', which='major', pad=pad_size, labelsize=0, direction='in', top=True, right=True)
    axis.tick_params(axis='y', which='both', pad=pad_size, labelsize=0, direction='in', top=True, right=True)
    # axis.tick_params(axis='both', denominator='major', length=major_length, width=tick_width)
    # axis.tick_params(axis='both', denominator='minor', length=minor_length, width=tick_width)
    # axis.xaxis.set_minor_locator(AutoMinorLocator())
    axis.yaxis.set_minor_locator(AutoMinorLocator())
    axis.set_xticks(pos)

    axis.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    axis.set_ylim([-4, 4])
    axis.set_xlim([1.5, pos[p - 1] + 1.5])

    # axis.set_yscale('log')
    # axis.set_ylim([1e-5,1e4])

ax1.tick_params(axis='y', which='both', pad=pad_size, labelsize=tick_size, direction='in', top=True, right=True)

ax3.tick_params(axis='y', which='both', pad=pad_size, labelsize=tick_size, direction='in', top=True, right=True)
ax3.tick_params(axis='x', which='both', pad=pad_size, labelsize=label_size, direction='in', top=True, right=True)
ax3.set_xticklabels([x.replace('2', '$_2$').replace('4', '$_4$') for x in species_list])

ax4.tick_params(axis='x', which='both', pad=pad_size, labelsize=label_size, direction='in', top=True, right=True)
ax4.set_xticklabels([x.replace('2', '$_2$').replace('4', '$_4$') for x in species_list])

ax1.set_title('Primordial', fontsize=label_size)
ax2.set_title('Volcanic', fontsize=label_size)
ax1.set_ylabel(r'log($\frac{X_i}{X_{i, FMQ}}$)', fontsize=label_size)
ax3.set_ylabel(r'log($\frac{X_i}{X_{i, IW-2}}$)', fontsize=label_size)

ax2.legend([r'$X_i < 10^{-6}$'], loc='upper right', fontsize=leg_font, frameon=False)

gs1.tight_layout(fig)
gs1.update(wspace=0.1, hspace=0.1)

# fig, *axs = dark_background(fig, [ax1, ax2, ax3, ax4])
fig.savefig(figpath + 'violins.pdf', bbox_inches='tight', dpi=300)

plt.show()










