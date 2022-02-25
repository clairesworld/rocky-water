""" test variations with temperature profile"""
import numpy as np
import perple_x as px
import matplotlib.pyplot as plt
from useful_and_bespoke import colorize, dark_background, cornertext
from parameters import M_E, M_SiO2, M_MgO, TO
import matplotlib
from matplotlib import rc
from matplotlib.pyplot import rcParams
import saturation as sat
from matplotlib.gridspec import GridSpec

# rc('text', usetex=True)  # turn off for running over ssh
# rcParams['font.family'] = 'serif'
# rcParams['font.serif'] = 'CMU Serif'

labelsize = 25
oxides = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'Na2O', 'FeO']
solution_phases = ['O', 'Sp', 'Cpx', 'Wad', 'Ring', 'Pv', 'Wus', 'C2/c', 'Opx', 'Aki', 'Ppv', 'Gt', 'CF',
                   'Pl', 'NAl']

all_phases = ['ol', 'opx', 'cpx', 'pl', 'gt', 'sp', 'neph', 'qtz', 'coes', 'hpcpx', 'st', 'cf', 'seif', 'wus', 'wad',
              'ring', 'capv', 'pv', 'ppv']



def plot_temperature_dependence(star=None,
                                    Tp_min=1500, Tp_max=2000, res=10, M_p=M_E, CMF=0.33, wt_oxides=px.wt_oxides_Earth,
                                wt_percent=True, labelsize=20, cmap='Pastel1', plot_all=False,
                                   legsize=16, ticksize=14, xpad=20, ypad=20,
                                    cmap_phases='tab20', linec='xkcd:navy', linew=2, xlims=None,
                                    save=True, **kwargs):
    """ make plot for several masses showing mantle structure and relative water in each zone"""
    fig, axes = plt.subplots(1, 1)
    ax = axes
    cmap_vals = matplotlib.cm.get_cmap(cmap)
    w = []
    w_upper = []
    w_trans = []
    w_lower = []
    Tp = np.linspace(Tp_min, Tp_max, num=res)
    for ii, Tpi in enumerate(Tp):

        dat = px.build_planet(name='Tp' + str(ii), Tp=Tpi, M_p=M_p, test_CMF=CMF, test_oxides=wt_oxides, star=star,
                          oxides=oxides, solution_phases=solution_phases, plot=plot_all, store_all=True, profile='adiabat',
                          **kwargs)

        i_wad = np.argmax(dat.df_all['X_wad'].to_numpy() > 0)
        i_pv = np.argmax(dat.df_all['X_pv'].to_numpy() > 0)
        m_h2o = dat.df_all['mass_h2o(kg)']
        m_tot = dat.df_all['mass(kg)']
        w_max = np.sum(m_h2o)

        if wt_percent:
            w.append(w_max / np.sum(m_tot) * 1e2)
            w_upper.append(np.sum(m_h2o[:i_wad]) / np.sum(m_tot[:i_wad]) * 1e2)
            w_trans.append(np.sum(m_h2o[i_wad:i_pv]) / np.sum(m_tot[i_wad:i_pv]) * 1e2)
            w_lower.append(np.sum(m_h2o[i_pv:]) / np.sum(m_tot[i_pv:]) * 1e2)
            ylabel = 'Water capacity (wt.%)'
        else:
            w.append(w_max / TO)
            w_upper.append(np.sum(m_h2o[:i_wad]) / TO)
            w_trans.append(np.sum(m_h2o[i_wad:i_pv]) / TO)
            w_lower.append(np.sum(m_h2o[i_pv:]) / TO)
            ylabel = 'Water capacity (Earth oceans)'
            ax.set_ylim(0, 2.7)

    ax.set_xlim(np.min(Tp), np.max(Tp))
    ax.plot(Tp, np.array(w), label='Total', c='w', lw=5)
    ax.plot(Tp, np.array(w_upper), label='Upper mantle', c='xkcd:lilac', lw=2)
    ax.plot(Tp, np.array(w_trans), label='Transition zone', c='xkcd:pale blue', lw=2)
    ax.plot(Tp, np.array(w_lower), label='Lower mantle', c='xkcd:gold', lw=2)
    ax.set_xlabel('Potential temperature (K)', fontsize=labelsize, labelpad=xpad)
    ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ypad)
    # ax.set_title(str(M_p/M_E) + ' $M_E$')
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.legend(frameon=False, fontsize=legsize)
    plt.tight_layout()
    fig, *ax = dark_background(fig, ax, )
    fig.savefig(px.fig_path + 'Tp_dependence.png', bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.show()

def sat_difference(name=None, star=None, Tp=1800, M_p=M_E, CMF=0.33, wt_oxides=px.wt_oxides_Earth, **kwargs):
    # get difference in water saturation content between solidus and some potential temperature
    if name is None:
        name = 'Earth'
    dat0 = px.build_planet(name=name + '_Ts', Tp=Tp, M_p=M_p, test_CMF=CMF, test_oxides=wt_oxides, star=star,
                          plot=True, profile='solidus', **kwargs)
    w_max_Ts = dat0.c_h2o_tot
    dat1 = px.build_planet(name=name + '_Tp', Tp=Tp, M_p=M_p, test_CMF=CMF, test_oxides=wt_oxides, star=star,
                          plot=True, **kwargs)
    w_max_Tp = dat1.c_h2o_tot
    print('mass at Ts', w_max_Ts/TO, 'TO')
    print('mass at Tp', w_max_Tp/TO, 'TO')
    print('diff', (w_max_Ts - w_max_Tp)/TO, 'TO')

plot_temperature_dependence(star=None, wt_percent=False,
                            Tp_min=1500, Tp_max=2000, res=10, M_p=M_E, CMF=0.33, wt_oxides=px.wt_oxides_Earth,
                            tol=0.01, maxIter=20, plot_all=False)

# sat_difference(Tp=1800, tol=0.01)

