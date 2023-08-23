# recreate adibekyan fig 2
# vary core iron partitioning randomly and see what intrinsic scatter you get for CMF

from pathlib import Path
import sys

path_root = Path('/home/claire/Works/')
sys.path.append(str(path_root))
print(sys.path)
import numpy as np
import py.perplexdata as px
import py.parameters as p
import py.main as rw
import py.ask_hypatia as hyp
import py.bulk_composition as bulk
import eos
import os
import matplotlib.pyplot as plt
import py.plot_perplex as plotpx
import pickle as pkl
import matplotlib.colors
from useful_and_bespoke import colourbar
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.lines as mlines
import matplotlib
from useful_and_bespoke import hex_to_rgb
from adibekyan2021.main import planet_composition
import pandas as pd


def planets_from_hypatia_random_coreeff(core_eff_min=0.7, core_eff_max=0.99, n_sample=-1, Mp=1, Tp=1600,
                                        restart=None,
                                        stopafter=None, skip_stars=[], maxIter=30, n='auto', tol=1e-4, **kwargs):
    """ restart is string name of last star that worked"""

    sample_names = hyp.random_star(n_sample, **kwargs)  # returns list of strings being star ids
    # print('skip stars')
    for st in skip_stars:
        print(st, 'in skip_stars, removing from sample')
        sample_names.remove(st)
    planets = []

    if restart is not None and n_sample == -1:
        # must be looping all stars for this to work, otfherwise can't work as will be random names
        ii_start = sample_names.index(restart)
    else:
        ii_start = 0
    if stopafter is not None and n_sample == -1:
        ii_last = sample_names.index(stopafter)
        subsample_names = sample_names[ii_start:ii_last + 1]
    else:
        subsample_names = sample_names[ii_start:]

    for ii, star in enumerate(subsample_names):
        print(ii + 1, '/', len(subsample_names))

        core_eff = np.random.uniform(core_eff_min, core_eff_max, 1)[0]  # np.random.uniform returns an array
        planet_dict = {'M_p': Mp * p.M_E, 'Tp': Tp, 'core_efficiency': core_eff,
                       'maxIter': maxIter, 'tol': tol, 'n': n,
                       }

        pl = rw.build_planet(star=star, plot_all=False, **planet_dict, verbose=True,
                             vertex_data='stx21ver', option_file='perplex_option_claire', **kwargs)
        if pl is not None:
            pl.write_star_composition(fname='nH_star.txt',
                                      path=pl.output_path)  # save parameter nH_star (abundances) to file
            planets.append(pl)
    return planets


def run_random_coreeff(core_eff_min=0.7, core_eff_max=0.99, Tp=1600, Mp=1, n_sample=-1, n='auto', restart=None,
                       perplex_path='/raid1/cmg76/perple_x/', tail='',  # apollo hires defaults
                       use_local_composition=False, **kwargs):
    # run at higher res over masses (primarily to get upper mantle)
    dist_name = 'U' + str(int(core_eff_min * 100)) + '-' + str(int(core_eff_max * 100))
    output_parent_path = perplex_path + 'output/random_coreeff/' + dist_name + '/'

    planets = planets_from_hypatia_random_coreeff(core_eff_min=core_eff_min, core_eff_max=core_eff_max,
                                                  n_sample=n_sample, pass_all=True,
                                                  restart=restart,
                                                  use_local_composition=use_local_composition,
                                                  perplex_path=perplex_path,
                                                  output_parent_path=output_parent_path,
                                                  n=n, Mp=Mp, Tp=Tp, **kwargs)
    return planets


def load_random_coreeff(output_path, verbose=False, **kwargs):
    print('reading', output_path)
    subfolders = [f.path for f in os.scandir(output_path) if f.is_dir()]
    dats = []
    count = 0
    for dir in subfolders:
        try:
            with open(dir + '/dat.pkl', "rb") as pfile:
                dat = pkl.load(pfile)
            dats.append(dat)
            count += 1
        except FileNotFoundError:
            if verbose:
                print('warning:', dir, 'contains no dat.pkl')
        except PermissionError as e:
            print(dir, e)
    return dats


def plot_density_mgfe(pl_list, fig=None, ax=None, rho_scale=p.rho_E, save=True, fig_path=plotpx.fig_path,
                      fformat='.pdf', vmin=0.1, vmax=0.4, cmap='pink', oxide_list=None, show_earth=True, show_exo=False,
                      labelsize=16, ticksize=14, f_iron_type='abundances', xlim=None, ylim=None, xlabelpad=None,
                      typical_uncert='a23', err_c='0.85', fname='adibekyan', ylabel=None, xlabel=None, add_cbar=True,
                      y_var='rho',
                      **kwargs):
    # like Adibekyan fig 2
    cmf = []
    m, r = [], []
    ys = []

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    for pl in pl_list:
        if y_var == 'rho':
            rho = pl.M_p / (4 / 3 * np.pi * pl.R_p ** 3)
        elif y_var == 'mass':
            rho = pl.M_p / p.M_E

        stellar_abundances = bulk.get_stellar_percent_abundance(pl.nH_star[0:5], which='mass', oxide_list=oxide_list)
        if f_iron_type == 'abundances':
            f_star_iron = stellar_abundances[pl.oxide_list.index('FeO')]
        elif f_iron_type == 'stoichiometric':
            mg = 10 ** (pl.nH_star[pl.oxide_list.index('MgO')] + 12)
            si = 10 ** (pl.nH_star[pl.oxide_list.index('SiO2')] + 12)
            fe = 10 ** (pl.nH_star[pl.oxide_list.index('FeO')] + 12)
            f_star_iron = planet_composition(None, None, None, Mg_abs=mg, Si_abs=si, Fe_abs=fe)

        ax.scatter(f_star_iron, rho / rho_scale, c=pl.CMF, vmin=vmin, vmax=vmax, cmap=cmap,
                   s=80, alpha=0.8, edgecolors='k', zorder=199, )  # plot on top

        cmf.append(pl.CMF)
        m.append(pl.M_p)
        r.append(pl.R_p)
        ys.append(rho / rho_scale)

    if show_earth:
        earth = rw.build_planet(M_p=1 * p.M_E, test_CMF=0.325, test_oxides=px.wt_oxides_MD95,
                                maxIter=30, tol=tol, n=n, Tp=Tp,
                                plot_all=False, get_saturation=True, verbose=True, clean=True,
                                vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                                name='Earth_' + str(Tp) + 'K',
                                )
        rho = earth.M_p / (4 / 3 * np.pi * earth.R_p ** 3)
        if f_iron_type == 'abundances':
            f_star_iron = stellar_abundances[earth.oxide_list.index('FeO')]
        elif f_iron_type == 'stoichiometric':
            mg = 10 ** (p.mg_sol + 12)
            si = 10 ** (p.si_sol + 12)
            fe = 10 ** (p.fe_sol + 12)
            f_star_iron = planet_composition(None, None, None, Mg_abs=mg, Si_abs=si, Fe_abs=fe)
        print('Earth x, y =', f_star_iron, rho / rho_scale)
        print('Earth M, R', earth.M_p, earth.R_p)
        ax.scatter(f_star_iron, 1, c=0.325, vmin=vmin, vmax=vmax, cmap=cmap, marker='$\oplus$',
                   s=300, zorder=200)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if typical_uncert:
        # show typical uncertainty
        yerr_a21 = np.mean(
            [0.42, 0.51, 0.33, 0.28, 0.35, 0.06, 0.14, 0.12, 0.13, 0.15, 0.24, 0.14, 0.25, 0.2, 0.25, 0.14,
             0.17, 0.30, 0.17, 0.15, 0.21, 0.14])  # adibekyan density error table s4
        print('yerr_a21', yerr_a21)

        # or use cayman densities
        yerr_u23_rocky_gcm3 = np.mean([0.3, 0.33, 0.23, 0.25, 0.24, 0.79, 0.74, 0.22, 0.88])
        yerr_u23_all_gcm3 = np.mean(
            [0.3, 0.33, 0.23, 0.25, 0.24, 0.79, 0.74, 0.22, 0.88, 1.37, 2.3, 1.7, 1.29, 1.33, 1.21,
             2.23, 1.09, 2.04, 1.19, 3.11, 2.16, 2.70, 2.79, 0.6, 0.67, 2.12, 0.41, 2.8, 1.51, 1.22,
             0.87, 1.51, 1.09, 1.77, 1.12, 2.17, 1.92, 2.89, 2.74, 2.47, 0.91, 0.88, 0.84, 0.99,
             1.07, 0.98, 0.81, 5.88, 0.92, 7.14, 0.93, 2.46])
        yerr_u23 = yerr_u23_rocky_gcm3 / (p.rho_E * 1e-3)  # scale by Earth bulk density in g/cm3
        print('yerr_u23', yerr_u23)

        if typical_uncert == 'a23':
            xa, ya = 37, 0.97
            yerr = yerr_a21
            serr = 'typical uncert.\n(Adibekyan+21)'
        elif typical_uncert == 'u23':
            xa, ya = 37, 0.97
            yerr = yerr_u23
            serr = 'typical uncert.\n(Unterborn+23)'

        # xa = np.ptp(ax.get_xlim()) * 0.9 + ax.get_xlim()[0]  # 90% of axis
        xa = 25  # ax.get_xlim()[0] + 1
        ya = 1.3  # np.mean(ys)
        print('xa, ya, yerr', xa, ya, yerr)
        ax.errorbar(x=xa, y=ya, yerr=yerr, xerr=None, ecolor=err_c, capsize=3)
        ax.text(xa + 0.3, ya + yerr, serr, c=err_c, va='top', ha='left', fontsize=14)

    # label and colour
    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    if add_cbar:
        cbar = colourbar(vector=cmf, cmap=cmap, vmin=vmin, vmax=vmax, label='Core mass fraction', ax=ax,
                         labelsize=labelsize, ticksize=ticksize,  # ticks=[0.01, 0.1, 0.2, 0.3],
                         labelpad=25)
        print('cmf', cmf)
    # print('m (kg)', m)
    # print('r (m)', r)
    # print('rho', [mm / (4 / 3 * np.pi * rr ** 3) for mm, rr in zip(m, r)])
    # print('p.rho_E', rho_scale, 'kg/m3')

    if show_exo:
        fig, ax = plot_real_planets(fig, ax, y_var=y_var, **kwargs)

    if save:
        fig.savefig(fig_path + fname + fformat, bbox_inches='tight', dpi=100)
    return fig, ax


def plot_real_planets(fig, ax, c_a21='g', f_iron_type='stoichiometric', y_var='rho', **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    # add 22 planets from Adibekyan+ 2021 Table S4
    list_a21_SM = [
        # super Mercuries
        {'star': 'K2-38', 'planet': 'b', 'rho': 1.35, 'rho_err': 0.42, 'f_star_iron': 33.40, 'f_err': 2.14, 'mass': 7.31},
        {'star': 'K2-106', 'planet': 'b', 'rho': 1.54, 'rho_err': 0.51, 'f_star_iron': 35.35, 'f_err': 1.36, 'mass': 8.36},
        {'star': 'K2-229', 'planet': 'b', 'rho': 1.41, 'rho_err': 0.33, 'f_star_iron': 33.44, 'f_err': 1.28, 'mass': 2.59},
        {'star': 'Kepler-107', 'planet': 'c', 'rho': 1.45, 'rho_err': 0.28, 'f_star_iron': 34.79, 'f_err': 2.41, 'mass': 9.39},
        {'star': 'Kepler-406', 'planet': 'b', 'rho': 1.50, 'rho_err': 0.35, 'f_star_iron': 33.05,  'f_err': 1.31, 'mass': 6.39},
    ]

    list_a21_UE = [
        # likey volatile-rich or molten according to Unterborn+ 2023 (unlike t to be rocky)
        {'star': '55 Cnc', 'planet': 'e', 'rho': 0.75, 'rho_err': 0.06, 'f_star_iron': 30.22, 'f_err': 1.91, 'mass': 8.59},  # molten
        {'star': 'WASP-47', 'planet': 'e', 'rho': 0.89, 'rho_err': 0.21, 'f_star_iron': 30.79, 'f_err':1.86, 'mass': 9.22},  # volatile-rich
        {'star': 'TOI-561', 'planet': 'b', 'rho': 0.51, 'rho_err': 0.14, 'f_star_iron': 25.23, 'f_err':1.06, 'mass': 1.59},  # volatile-rich
    ]
    list_a21_rocky = [
        # or not included in Unterborn+ 2023 sample
        {'star': 'HD 213885', 'planet': 'b', 'rho': 1.06, 'rho_err': 0.12, 'f_star_iron': 33.95, 'f_err': 0.97, 'mass': 8.83},
        {'star': 'HD 219134', 'planet': 'b', 'rho': 0.96, 'rho_err': 0.13, 'f_star_iron': 32.45,  'f_err':1.68, 'mass': 4.27},
        {'star': 'HD 219134', 'planet': 'c', 'rho': 1.08, 'rho_err': 0.15, 'f_star_iron': 32.45, 'f_err': 1.68, 'mass': 3.96},  # U23 confirm rocky
        {'star': 'HD 3167', 'planet': 'b', 'rho': 0.74, 'rho_err': 0.24, 'f_star_iron': 32.15, 'f_err': 1.27, 'mass': 5.02},
        {'star': 'K2-141', 'planet': 'b', 'rho': 1.07, 'rho_err': 0.24, 'f_star_iron': 34.33, 'f_err': 3.51, 'mass': 5.09},
        {'star': 'K2-216', 'planet': 'b', 'rho': 0.89, 'rho_err': 0.35, 'f_star_iron': 34.68, 'f_err': 5.23, 'mass': 7.91},
        {'star': 'K2-265', 'planet': 'b', 'rho': 0.90, 'rho_err': 0.20, 'f_star_iron': 33.73, 'f_err': 1.15, 'mass': 6.54},
        {'star': 'K2-291', 'planet': 'b', 'rho': 1.12, 'rho_err': 0.25, 'f_star_iron': 34.10, 'f_err': 1.03, 'mass': 6.41},
        {'star': 'Kepler-10', 'planet': 'b', 'rho': 0.85, 'rho_err': 0.14, 'f_star_iron': 29.19, 'f_err': 0.82, 'mass': 3.33},
        {'star': 'Kepler-20', 'planet': 'b', 'rho': 0.92, 'rho_err': 0.17, 'f_star_iron': 32.42, 'f_err': 1.26, 'mass': 9.69},
        {'star': 'Kepler-93', 'planet': 'b', 'rho': 0.95, 'rho_err': 0.17, 'f_star_iron': 33.64, 'f_err': 0.96, 'mass': 4.0},
        {'star': 'TOI-402', 'planet': 'b', 'rho': 0.99, 'rho_err': 0.15, 'f_star_iron': 30.87, 'f_err': 1.38, 'mass': 7.21},
        {'star': 'EPIC 249893012', 'planet': 'b', 'rho': 0.76, 'rho_err': 0.14, 'f_star_iron': 31.33, 'f_err': 1.17, 'mass': 8.75},
        # super-Mercury, molten
        {'star': 'Kepler-78', 'planet': 'b', 'rho': 0.90, 'rho_err': 0.30, 'f_star_iron': 31.80, 'f_err': 2.40, 'mass': 1.69},
        # super-Mercury, molten
    ]

    alpha_unearth = 0.3
    if f_iron_type == 'stoichiometric':
        for d in list_a21_rocky:
            ax.errorbar(x=d['f_star_iron'], y=d[y_var], xerr=d['f_err'], yerr=d['rho_err'], fmt='s',
                        c=c_a21 + [1], ecolor=c_a21 + [1], markeredgecolor=c_a21 + [1],
                        elinewidth=1, markersize=7, zorder=2)
        for d in list_a21_UE:
            ax.errorbar(x=d['f_star_iron'], y=d[y_var], xerr=d['f_err'], yerr=d['rho_err'], fmt='D',
                        c=c_a21 + [alpha_unearth], ecolor=c_a21 + [alpha_unearth], markeredgecolor=c_a21 + [1],
                        elinewidth=1, markersize=7, zorder=1)
        for d in list_a21_SM:
            ax.errorbar(x=d['f_star_iron'], y=d[y_var], xerr=d['f_err'], yerr=d['rho_err'], fmt='v',
                        c=c_a21 + [alpha_unearth], ecolor=c_a21 + [alpha_unearth], markeredgecolor=c_a21 + [1],
                        elinewidth=1, markersize=7, zorder=1)
    else:
        raise NotImplementedError('need to calculate stellar abundance Fe ratio')

    # # add planets from Unterborn+ 2023
    # df = pd.read_csv('/home/claire/Works/unterborn_2023_apj.dat', delimiter=r"\s+", header=0, index_col=None)
    # print(df.head())
    #
    # # nominally rocky
    # rows = range(10)
    #
    # for ii in rows:
    #     row = df.iloc[ii]
    #
    # none of these planets have measured host-star compositons

    # unlikely rocky

    return fig, ax


def subplot_density_mgfe(pl_list, save=True, fig_path=plotpx.fig_path, fformat='.pdf', ylim0=None, ylim1=None, cmf=[],
                         typical_uncert='a23', fname='adibekyan', ylabel=None, xlabel=None, labelsize=16, c_a21='g',
                         legsize=12, vmin=0.0, vmax=0.35, cmap='copper', **kwargs):
    fig, axes = plt.subplots(2, 1, figsize=(9, 10))

    # zoom in fake planets
    fig, axes[0] = plot_density_mgfe(pl_list, fig=fig, ax=axes[0], save=False, show_exo=False,
                                     vmin=vmin, vmax=vmax, cmap=cmap,
                                     typical_uncert=None, labelsize=labelsize,
                                     ylim=ylim0, ylabel='', xlabel='', add_cbar=False, **kwargs)

    # compare to real
    fig, axes[1] = plot_density_mgfe(pl_list, fig=fig, ax=axes[1], save=False, show_exo=True,
                                     vmin=vmin, vmax=vmax, cmap=cmap,
                                     typical_uncert=None, labelsize=labelsize, c_a21=c_a21,
                                     ylim=ylim1, ylabel='', xlabel=xlabel, add_cbar=False, **kwargs)

    # fix labels and adjust
    axes[0].set_xlabel('')
    axes[0].set_xticklabels([])
    fig.supylabel(ylabel, fontsize=labelsize, x=0.03)
    plt.subplots_adjust(hspace=0.05)

    # add double cbar
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.95, 0.25, 0.03, 0.5])  # left, bottom, width, height
    cbar = colourbar(vector=np.array(cmf) * 100, label=r'Core mass fraction (%)', ax=axes[0], cax=cbar_ax, cmap=cmap,
                     vmin=vmin * 100, vmax=vmax * 100,
                     labelsize=labelsize,  # ticks=[0.01, 0.1, 0.2, 0.3],
                     labelpad=25, **kwargs)

    # add legend
    alpha_unearth = 0.4
    handles = [
        mlines.Line2D([], [], color=matplotlib.cm.get_cmap(cmap)(1), marker='o', lw=0, alpha=0.6, mec='k',
                      markersize=9, label='Modelled pure rock-iron planets'),
        mlines.Line2D([], [], color=c_a21, marker='s', lw=0, alpha=1,
                      markersize=6, label='Known exoplanets, indistinguishable'),
        mlines.Line2D([], [], marker='v', lw=0, c=c_a21 + [alpha_unearth], markeredgecolor=c_a21 + [1],
                      markersize=6, label='Known exoplanets, likely super-Mercuries'),
        mlines.Line2D([], [],  marker='D', lw=0, c=c_a21 + [alpha_unearth], markeredgecolor=c_a21 + [1],
                      markersize=5, label='Known exoplanets, likely volatile-rich'),
        mlines.Line2D([], [],
                      color=matplotlib.cm.get_cmap(cmap)(matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)(0.325)),
                      marker='$\oplus$', lw=0, alpha=1,
                      markersize=13, label='Earth')
    ]
    axes[1].legend(frameon=False, fontsize=legsize, handles=handles, loc='upper left',
                   # bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left"
                   )

    if save:
        fig.savefig(fig_path + fname + fformat, bbox_inches='tight', dpi=100)
    return fig, axes


fig_path = '/home/claire/Works/thesis/v0/Introduction/Figs/'  #plotpx.fig_path
oxide_list = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO']
# tol = 1e-2  # quick
# n = 800  # quick
tol = 1e-4  # better
n = 'auto'  # better
Mp = 1
Tp = 1600
cmf = [0.0738050909226403, 0.03619148074451505, 0.12007514003807998, 0.17378562263487016, 0.06462788368743537,
       0.13899073613366894, 0.12630556656140685, 0.0452846447032519, 0.07308845481200626, 0.23371944605002945,
       0.11768004561630445, 0.045465581082614975, 0.10741652527389006, 0.07339678740152966, 0.12495548027361834,
       0.24368645412869738, 0.03225559705720619, 0.236044201462323, 0.2876529926355693, 0.2925143294950069,
       0.28412323793044864, 0.13905755245882184, 0.30638864972386315, 0.21257844297329453, 0.2878616551602129,
       0.03856549154597606, 0.22239498247836942, 0.03890263328841183, 0.2033648769320962, 0.031420722930920145,
       0.2655670164410848, 0.0936672642985304, 0.30680143945078425, 0.24724128750216268, 0.13728714658094113,
       0.12169659784100642, 0.18041300182269193, 0.22317810246857914, 0.10301629784049325, 0.1958789620515503,
       0.05318138097824966, 0.23851621679092497, 0.16320364072575289, 0.11529597685823674, 0.29954266139744623,
       0.20048128906383827, 0.05396980936322031, 0.29419113580436745, 0.10829907284599814, 0.06144546158034872,
       0.27577863325637947, 0.0325625694822019, 0.20970907194817198, 0.042844947679461606, 0.24607989992343196,
       0.1220053677858886, 0.19554813765270182, 0.029090116333161622, 0.2366612099014872, 0.11782501352948863,
       0.22510129033089907, 0.1911459539712127, 0.1521548787472902, 0.05008225032995158, 0.14271196225698948,
       0.1761936927411093, 0.24303400262865765, 0.2662443100513557, 0.11804136716693296, 0.1532150940378515,
       0.27214227370328303, 0.15406087507510066, 0.2575242224761499, 0.22774464494582727, 0.08823320215061571,
       0.13003305232717122, 0.13770047271716127, 0.03816920751678421, 0.15000158170590439, 0.16639156961510715,
       0.3212696846235311, 0.26578858750974005, 0.04628854485991189, 0.2505638247297602, 0.22489281383280438,
       0.23686285940038043, 0.052285416455418426, 0.038956259566535555, 0.041954413843387835, 0.04214022947863382,
       0.1910822791672293, 0.28105373123641636, 0.2562737086125053, 0.30448148940036274, 0.23592436246254156,
       0.1817938898769503, 0.08433560268727777, 0.09928640046560544, 0.3311224065638243, 0.32017285914822863,
       0.2802586249355036, 0.21226180284327983, 0.30751008398837615, 0.2539162467174997, 0.20510342755275932,
       0.20225924757612992, 0.2894786850705112, 0.08929975534238914, 0.17874871306071183, 0.21935679370810235,
       0.15242616090720296, 0.29110853843828993, 0.1910616513989296, 0.10795674679597661, 0.2297568125089314,
       0.2760412238229401, 0.24445937374126372, 0.07343821249840539, 0.13099028441090232, 0.18544289264136174,
       0.157355978950786, 0.27832555561413547, 0.21928879947623262, 0.19546206388025664, 0.190597713574639,
       0.22866314815154268, 0.18999215947911316, 0.06714323195370532, 0.19785987306468172, 0.26699321433264595,
       0.06443498646593061, 0.09902383218891715, 0.27986041943294315, 0.08420006820384619, 0.14691631109655018,
       0.16397074862370942, 0.23321103490932024, 0.1352657426511494, 0.27002468021941856, 0.23344331362655985,
       0.21850304365144638, 0.26894078343950584, 0.11313342695047586, 0.08294324439483239, 0.3942745597129296,
       0.24460507243152563, 0.17823066205693994, 0.1761897233407403, 0.23708512927477937, 0.2613299714582581,
       0.10579574034673805, 0.10264220723822902, 0.18458330550249133, 0.2336085434325723, 0.09634994100823419,
       0.1843462260772796, 0.2540942223678578, 0.3061330912089214, 0.1300802950989623, 0.20391892599274672,
       0.10794245840432874, 0.03372305456658646, 0.17641916700809226, 0.3416489513124915, 0.2706022490852644,
       0.04225498742268097, 0.20698443983272466, 0.14077375435386305, 0.27042252357416086, 0.08980638546920114,
       0.23935955357813463, 0.10647275306041029, 0.032869900957063755, 0.06271860825895303, 0.11677856689931079,
       0.039660692145260724, 0.2802367498453548]

# pl_list = run_random_coreeff(core_eff_min=0.1, core_eff_max=0.9999, Tp=Tp, Mp=Mp, n_sample=150, n=n, tol=tol,
#                              restart=None, oxide_list=oxide_list, oxides=oxide_list,
#                              use_local_composition=True,
#                              existing_output_parent='/home/claire/Works/perple_x/output/stellar_comps/',
#                              existing_dir='stellar_comps/',
#                              local_composition_cuttoidx=5,
#                              suppress_output=True,
#                              perplex_path='/home/claire/Works/perple_x/'
#                              # perplex_path='/raid1/cmg76/perple_x/'
#                              )

pl_list = load_random_coreeff('/home/claire/Works/perple_x/output/random_coreeff/U10-99/')

rgb = hex_to_rgb('#60ce92')
fig, ax = subplot_density_mgfe(pl_list, oxide_list=oxide_list, save=True, fig_path=fig_path,
                               rho_scale=p.rho_E,
                               show_earth=True, cmf=cmf, c_a21=rgb,#'xkcd:cool green',
                               vmin=0.0, vmax=0.35, cmap='copper', err_c='xkcd:dark green blue',
                               legsize=12, labelsize=16, ticksize=12, xlabelpad=12,
                               # typical_uncert='u23',
                               ylim0=(0.91, 1.03),  # for no adibekyan data
                               typical_uncert='u23', ylim1=(0.4, 1.6),  # for adibekyan data
                               # typical_uncert='a21', ylim=(0.7, 1.25),
                               f_iron_type='stoichiometric', xlim=(23, 39.5),
                               # f_iron_type='abundances', xlim=(20, 50),
                               # ylabel=r'$f_{\rm iron}^{\rm star, refractory}$ (wt.%)',
                               xlabel=r'$f_{\rm iron}^{\rm star}$',
                               ylabel=r'$\rho/\rho_{\rm Earthlike}$',
                               fname='starplanet_stoich_subplot', fformat='.pdf'
                               )

# test mass versus f_iron
# fig, ax = plot_density_mgfe(pl_list, oxide_list=oxide_list, save=True, fig_path=fig_path,
#                              show_earth=False, show_exo=True, show_typical_uncert=True,
#                             vmin=0.01, vmax=0.35, cmap='copper', err_c='xkcd:dark green blue', c_a21=rgb,
#                             # typical_uncert='u23', ylim=(0.87, 1.08),  # for no adibekyan data
#                             typical_uncert=None, #ylim=(0.49, 1.55),  # for adibekyan data
#                             # typical_uncert='a21', ylim=(0.7, 1.25),
#                             f_iron_type='stoichiometric', xlim=(24, 36),
#                             # f_iron_type='abundances', xlim=(20, 50),
#                             # ylabel=r'$f_{\rm iron}^{\rm star, refractory}$ (wt.%)',
#                             ylabel=r'$M_p$ ($M_E$)',
#                             xlabel=r'$f_{\rm iron}^{\rm star}$',
#                             y_var='mass', rho_scale=1, ylim=(0, 14),
#                             fname='starplanet_stoich_mass', fformat='.pdf'
#                             )
# not a real relationship

plt.show()
