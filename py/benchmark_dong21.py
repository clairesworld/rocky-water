""" make plots like Dong+ 2021 for benchmarking """

import numpy as np
import perplexdata as px
import matplotlib.pyplot as plt
from useful_and_bespoke import colorize
from parameters import M_E
import matplotlib

labelsize = 14
oxides = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'Na2O', 'FeO']
solution_phases = ['O', 'Sp', 'Cpx', 'Wad', 'Ring', 'Pv', 'Wus', 'C2/c', 'Opx', 'Aki', 'Ppv', 'Gt', 'CF',
                   'Pl', 'NAl'
                   ]


def figure_3(num=9, save=True):
    cmap = 'coolwarm'
    Ts_vec = np.linspace(1550, 1950, num=num)
    c = colorize(Ts_vec, cmap=cmap)[0]
    names = ['Ts' + str(int(s)) for s in Ts_vec]
    dat_list = px.build_multi_planets(loop_var_name='Tsurf', loop_var_vals=Ts_vec, plot=False, names=names,
                                      get_saturation=True, M_p=M_E, oxides=oxides, solution_phases=solution_phases,
                                      test_CMF=0.33, test_oxides=px.wt_oxides_Earth, overwrite=True)

    fig, axes = plt.subplots(1, 2)
    for ii, dat in enumerate(dat_list):
        fig, axes[0] = dat.profile('T(K)', yvar='p', xlabel='Temperature (K)', xscale=1,
                                   lw=2, c=c[ii], logx=False, reverse_y=True, labelsize=labelsize,
                                   label='Tp = ' + str(Ts_vec[ii]) + ' K',
                                   y2var=None, y2label=None, y2scale=None, save=False, fig=fig, ax=axes[0])
        fig, axes[1] = dat.profile('c_h2o', yvar='p', xlabel='H$_2$O storage capacity (wt ppm)', xscale=1e6,
                                   lw=2, c=c[ii], logx=True, reverse_y=True, labelsize=labelsize,
                                   label='Tp = ' + str(Ts_vec[ii]) + ' K',
                                   y2var=None, y2label=None, y2scale=None, save=False, fig=fig, ax=axes[1])
    axes[0].legend()
    plt.tight_layout()
    if save:
        fig.savefig(px.fig_path + 'D21_fig3.png')

    ## check out some compositional profiles
    # for ii in [0, -1]:
    #     dat_list[ii].plot_composition(denominator='pressure', fig_path=px.fig_path, save=True)


def figure_4(Tp, labelsize=12, cmap='tab20', plot_composition=True, save=True,
             test_oxides=px.wt_oxides_Earth, test_CMF=0.33, M_p=M_E, **kwargs):
    col = matplotlib.cm.get_cmap(cmap)
    fig, axes = plt.subplots(2, 1)
    axes[0].set_title('$T_p$ = ' + str(Tp) + ' K', fontsize=labelsize)
    dat = px.build_planet(name='Tp' + str(Tp), Tsurf=Tp, M_p=M_p, test_CMF=test_CMF, test_oxides=test_oxides,
                          oxides=oxides, solution_phases=solution_phases, plot=False, store_all=True,
                          get_saturation=True, overwrite=True, **kwargs)
    x = dat.data['P(bar)'].to_numpy() * 1e-4  # GPa
    w = dat.data['c_h2o'].to_numpy() * 1e6  # ppm
    n_phases = len(dat.phases_px)

    axes[0].plot(x, w, c='xkcd:navy')
    axes[0].set_yscale('log')
    axes[0].set_ylabel('Storage capacity\n(ppm)', fontsize=labelsize)
    axes[0].set_ylim(1e0, 1e5)
    axes[0].set_yticks([1e1, 1e3, 1e5])

    water_col = [col for col in dat.data if col.startswith('w_')]
    water_fractions = dat.data[water_col].to_numpy()
    mode_arr = water_fractions / np.expand_dims(dat.data['c_h2o'].to_numpy(), axis=1)

    axes[1].stackplot(x, mode_arr.T * 100, labels=[col[2:] for col in dat.data if col.startswith('X_')],
                      colors=[col(j) for j in np.arange(1, n_phases + 1)/(n_phases + 1)])
    axes[1].set_ylabel('Water modality\nin NAMs', fontsize=labelsize)
    axes[1].set_xlabel('Pressure (GPa)', fontsize=labelsize)
    axes[1].set_ylim(0, 100)
    leg = axes[1].legend(bbox_to_anchor=(1.01, 2.3), loc='upper left', frameon=False, fontsize=10)

    for ax in axes:
        ax.set_xlim(x.min(), x.max())
    # plt.tight_layout()
    if save:
        fig.savefig(px.fig_path + 'exo_D21_fig4_' + str(Tp) + '.png', bbox_extra_artists=(leg,), bbox_inches='tight')

    if plot_composition:
        dat.plot_composition(save=save)



# figure_3()
figure_4(1600, star="HIP 54035", test_oxides=None)
# figure_4(1900)
plt.show()
