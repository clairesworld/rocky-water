import oxygen_fugacity_plots as pltfo2
import numpy as np
import perplexfugacitydata as perplex
import meltsfugacitydata as melts
import matplotlib.pyplot as plt
from py.useful_and_bespoke import dark_background, cornertext
import matplotlib

T_iso = 1373
p_min, p_max = 1e4, 4e4
T_min, T_max = 1372.5, 1900.5  # endpoint can't equal T_of_interest
pressures_of_interest = np.linspace(p_min, p_max, 12)  # for alphaMELTS
X_ferric = 0.03
star = None
name = 'Stolper'
plot_kwargs = {}
labelsize = 22
ticksize = 16

c_phases = [plt.cm.tab20b(14),  # gt
            plt.cm.tab20b(10),  # cpx
            plt.cm.tab20b(11),  # opx
            # plt.cm.tab20b(9),  # hpcpx
            plt.cm.tab20b(7),  # ol
            # plt.cm.tab20b(5),  # wad
            # plt.cm.tab20b(4),  # ring
            # plt.cm.tab20c(19),  # pv
            plt.cm.tab20b(19),  # qtz
            # plt.cm.tab20b(18),  # coes
            # plt.cm.tab20b(17),  # st
            # plt.cm.tab20b(3),  # wus
            # plt.cm.tab20b(0),  # capv
            ]
cmap = matplotlib.colors.ListedColormap(c_phases)

" try to recreate Stolper fig 3 - MELTS "

wt_oxides_DMM_ext = {'SiO2': 44.71, 'MgO': 38.73, 'CaO': 3.17, 'Al2O3': 3.98, 'FeO': 8.17,
                     'Na2O': 0.28, 'TiO2': 0.13, 'Cr2O3': 0.57}  # Workman & Hart depleted mantle, extended

ii = 0
opp = melts.output_parent_default + 'mgsi_from_earth/'

P = np.linspace(1e4, 4e4)  # bar
T = 1373.15  # K
LogfO2_QFM = (-24014 / T) + 8.555 + ((0.092 * (P - 1) / T))
O2_fugacity_QFM = 10 ** LogfO2_QFM

for ii, mgsi in enumerate(np.linspace(0.75, 1.5, num=25)[:1]):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    name = 'mgsi' + str(mgsi).replace('.', ',')
    fig, axes = pltfo2.stolper_subplot(name=name, output_parent_path=opp, show_buffer=True,
                                       save=False,
                                       # p_min=1.1,  # in GPa
                                       lw=2, linec='xkcd:british racing green', c_buf='w', cmap=cmap,
                                       model='melts', fname='dmm_fo2_subplot_mlt_dark', fig=fig, axes=axes,
                                       labelsize=labelsize, legsize=16,
                                       **plot_kwargs)

    axes[1].set_ylim(0, 100)
    axes[0].set_ylim(-11, -6)
    axes[1] = cornertext(axes[1], '1000$^\circ$C isotherm', color='k', size=labelsize, pos='top right')
    fig.suptitle('Mg/Si = {:.2f}'.format(mgsi), color='xkcd:off white', fontsize=labelsize)

    for ax in axes:
        ax.set_xlim(1, 4)
        ax.set_xlabel('Pressure (GPa)', fontsize=labelsize, labelpad=12)
        ax.tick_params(axis='both', labelsize=ticksize)
        fig, ax1 = dark_background(fig, ax)

    print('QFM', LogfO2_QFM)
    axes[0].plot(P, LogfO2_QFM, 'r')
    plt.subplots_adjust(wspace=0.2)
    fig.savefig(pltfo2.figpath + 'anim/modes' + str(ii) + '.png', bbox_inches='tight')
    plt.show()
