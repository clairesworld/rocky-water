""" test variations with stellar composition"""
import numpy as np
import perple_x as px
import matplotlib.pyplot as plt
from useful_and_bespoke import colorize, dark_background, cornertext
from parameters import M_E, M_SiO2, M_MgO
import matplotlib
from matplotlib import rc
from matplotlib.pyplot import rcParams
import saturation as sat

rc('text', usetex=True)  # turn off for running over ssh
# rcParams['font.family'] = 'serif'
# rcParams['font.serif'] = 'CMU Serif'

labelsize = 25
oxides = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'Na2O', 'FeO']
solution_phases = ['O', 'Sp', 'Cpx', 'Wad', 'Ring', 'Pv', 'Wus', 'C2/c', 'Opx', 'Aki', 'Ppv', 'Gt', 'CF',
                   'Pl', 'NAl']

all_phases = ['ol', 'opx', 'cpx', 'pl', 'gt', 'sp', 'neph', 'qtz', 'coes', 'hpcpx', 'st', 'cf', 'seif', 'wus', 'wad',
              'ring', 'capv', 'pv', 'ppv']

# wt_oxides_Earth = {'SiO2': 44.48, 'MgO': 39.22, 'CaO': 3.44, 'Al2O3': 3.59, 'FeO': 8.10, 'Na2O': 0.36}
def update_MgSi(MgSi=None, oxides=px.wt_oxides_Earth):
    # vary Mg/Si, use original otherwise
    # update m_MgO and m_SiO2, preserving total mass
    m_tot = oxides['MgO'] + oxides['SiO2']  # total mass (%) of MgO and SiO2 - conserve this
    n_MgO = oxides['MgO'] / M_MgO  # convert to moles
    n_SiO2 = oxides['SiO2'] / M_SiO2
    n_MgSi_old = n_MgO / n_SiO2  # n_Mg = n_MgO etc
    if MgSi is None:
        MgSi = n_MgSi_old  # no update

    m_ratio_new = MgSi * (M_MgO / M_SiO2)  # m = nM, in % | new mass ratio
    x = m_ratio_new
    y = m_tot
    m_MgO_new = x*y/(1 + x)
    m_SiO2_new = m_tot - m_MgO_new
    oxides['MgO'] = m_MgO_new
    oxides['SiO2'] = m_SiO2_new
    return oxides


def plot_water_composition_profiles(name=None, fig_name=None, title=None, star=None,
                                    Tp=1600, M_p=M_E, CMF=0.33, wt_oxides=px.wt_oxides_Earth,
                                    labelsize=12, legsize=10, ticksize=12, xpad=10, ypad=10,
                                    cmap_phases='tab20', linec='xkcd:navy', linew=2, xlims=None,
                                    save=True, reorder=False,
                                    show_saturation=True, vertical=False):
    if name is None:
        name = 'Tp' + str(Tp)
    if fig_name is None:
        fig_name = name
    if title is None:
        title = '$T_p$ = ' + str(Tp) + ' K'

    cmap_vals = matplotlib.cm.get_cmap(cmap_phases)

    dat = px.build_planet(name=name, Tsurf=Tp, M_p=M_p, test_CMF=CMF, test_oxides=wt_oxides, star=star,
                          oxides=oxides, solution_phases=solution_phases, plot=False, store_all=True,
                          get_saturation=True, overwrite=True)

    x = dat.df_all['P(bar)'].to_numpy()*1e-4  # GPa
    w = dat.df_all['c_h2o'].to_numpy()*1e6  # ppm
    df_frac = dat.df_all[[col for col in dat.df_all if col.startswith('frac_h2o_')]]
    df_X = dat.df_all[[col for col in dat.df_all if col.startswith('X_')]]

    if reorder:   # so consistent between runs
        phase_names_orig = [col[2:] for col in dat.df_all if col.startswith('X_')]

        cols_frac_all = []
        for ph in all_phases:
            this_col = str('frac_h2o_' + ph)
            cols_frac_all.append(this_col)
            if ph not in phase_names_orig:  # put 0s in df
                df_frac[this_col] = 1e-10
        df_frac = df_frac[cols_frac_all]

        cols_X_all = []
        for ph in all_phases:
            this_col = str('X_' + ph)
            cols_X_all.append(this_col)
            if ph not in phase_names_orig:  # put 0s in df
                df_X[this_col] = 1e-20
        df_X = df_X[cols_X_all]

        n_phases = len(all_phases)
        colours = [cmap_vals(j) for j in np.arange(0, n_phases) / (n_phases)]
        phase_names = all_phases

    else:
        phase_names = [col[2:] for col in dat.df_all if col.startswith('X_')]
        n_phases = len(phase_names)
        colours = [cmap_vals(j) for j in np.arange(1, n_phases + 1) / (n_phases + 1)]

    # set up figure
    if show_saturation:
        fig, axes = plt.subplots(3, 1, figsize=(20, 20))
        ax_sat = axes[1]
        if vertical:
            print('vertical axes not implemented with stackplot')
    else:
        if not vertical:
            fig, ax = plt.subplots(1, 1, figsize=(20, 8))
        else:
            fig, ax = plt.subplots(1, 1, figsize=(8, 20))
        ax_sat = ax
        axes = [ax]
    axes[0].set_title(title, fontsize=labelsize, c=linec)

    # plot storage capacity per layer
    if show_saturation:
        axes[0].plot(x, w, c=linec, lw=linew)
        axes[0].set_yscale('log')
        axes[0].set_ylabel('Storage capacity\n(ppm)', fontsize=labelsize, labelpad=ypad)
        axes[0].set_ylim(1e0, 1e5)
        axes[0].set_yticks([1e1, 1e3, 1e5])
        m_w_tot = sat.total_water(dat.df_all) * dat.M_p / sat.TO
        axes[0] = cornertext(axes[0], 'total = {0:3.1f} earth oceans'.format(m_w_tot), size=labelsize, c=linec)

        # plot modal water
        axes[2].stackplot(x, df_frac.to_numpy().T*100, labels=phase_names,
                          colors=colours)
        axes[2].set_ylabel('Water modality\nin NAMs', fontsize=labelsize, labelpad=ypad)
        axes[2].set_ylim(0, 100)
        leg = axes[2].legend(bbox_to_anchor=(1.01, 0), loc='lower left', frameon=False, fontsize=legsize)

    # plot composition profile
    for ii, col in enumerate(df_X.columns):
        y = df_X[col].to_numpy() * 100
        if vertical and not show_saturation:
            ax_sat.plot(y, x, c=colours[ii], lw=linew, label=col[2:])
        else:
            ax_sat.plot(x, y, c=colours[ii], lw=linew)
    if vertical and not show_saturation:
        ax_sat.set_xlabel('Molar abundance', fontsize=labelsize, labelpad=ypad)
        axes[-1].set_ylabel('Pressure (GPa)', fontsize=labelsize, labelpad=xpad)
        ax_sat.set_ylim(x.min(), x.max())  # pressure min, max
        ax_sat.invert_yaxis()
        leg = ax_sat.legend(bbox_to_anchor=(1.01, 0), loc='lower left', frameon=False, fontsize=legsize)
    else:
        ax_sat.set_ylabel('Molar abundance', fontsize=labelsize, labelpad=ypad)
        ax_sat.set_ylim(2, 100)
        axes[-1].set_xlabel('Pressure (GPa)', fontsize=labelsize, labelpad=xpad)

    if xlims is None:
        xlims = (x.min(), x.max())
    for ii, ax in enumerate(axes):
        ax.set_xlim(xlims)
        ax.tick_params(axis='both', which='major', labelsize=ticksize)
        if ii != len(axes) - 1:
            ax.set_xticks([])
    # plt.tight_layout()
    fig, *axes = dark_background(fig, axes, )
    if save:
        fig.savefig(px.fig_path + fig_name + '.png', bbox_extra_artists=(leg,), bbox_inches='tight',
                    facecolor=fig.get_facecolor())

    return fig, axes


# MgSi = np.linspace(0.5, 1.5, num=50)
# for ii, rat in enumerate(MgSi):
#     wt_oxides = update_MgSi(rat, px.wt_oxides_Earth)  # modify oxide list
#     plot_water_composition_profiles(name='MgSi' + str(ii), fig_name=None,
#                                         Tp=1600, M_p=M_E, CMF=0.33, wt_oxides=wt_oxides,
#                                         labelsize=45, legsize=35, ticksize=30, xpad=30, ypad=30,
#                                     cmap_phases='tab20', title='Mg/Si = {0:2.2f}'.format(rat),
#                                     # xlims=(0, 140),
#                                         save=True, linec='w', linew=3, reorder=True)

for ii, star in enumerate(['sun', "HIP 54035"]):
    plot_water_composition_profiles(name=star.replace(" ", ""), star=star, fig_name=None,
                                        Tp=1600, M_p=M_E, CMF=0.33, wt_oxides=None,
                                        labelsize=45, legsize=35, ticksize=30, xpad=30, ypad=30,
                                    cmap_phases='tab20', title=['Solar', 'HD 95735'][ii],
                                    xlims=(0.01, 100),
                                        save=True, linec='w', linew=3, reorder=True, show_saturation=False,
                                    vertical=True)

# plot_water_composition_profiles(name='Earth1600', fig_name=None,
#                                     Tp=1600, M_p=M_E, CMF=0.33, wt_oxides=px.wt_oxides_Earth,
#                                     labelsize=45, legsize=35, ticksize=30, xpad=30, ypad=30,
#                                 cmap_phases='tab20', title='Earth',
#                                 xlims=(0.01, 100),
#                                     save=True, linec='w', linew=3, reorder=True, show_saturation=False,
#                                 vertical=True)