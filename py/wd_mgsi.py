# compare pwd and star mg/si

from pathlib import Path
import sys
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
import pandas as pd
import warnings

duplicates = ['GaiaJ0347+1624Photo',
              'GaiaJ0347+1624Spec',
              'GaiaJ0510+2315Photo',
              'GaiaJ0510+2315Spec',
              'GaiaJ0644-0352Photo',
              'GaiaJ0644-0352Spec',
              'GaiaJ2100+2122Photo',
              'GaiaJ2100+2122Spec',
              'WD0611-6931Spec',
              'WD1622+587Photo',
              'WD1622+587Spec',
              'SDSSJ0006+2858Photo',
              'SDSSJ0006+2858Spec',
              'GALEX1931+0117G12',  # UV
              'PG0843+516G12',  # UV
              'PG1015+161G12',  # UV
              'SDSSJ1228+1040'  # UV
              ]


# Gaensicke 2012 is UV


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def ratio_from_nH(mgh, sih):
    return 10 ** (mgh - sih)


# load hypatia
# oxide_list_default = ['SiO2', 'MgO', 'CaO', 'Al2O3', 'FeO']  # 'NaO
# subfolders = [f.path for f in os.scandir('/home/claire/Works/perple_x/output/stellar_comps/') if f.is_dir()]
# mgsi_hyp = []
# bad = []
# for path in subfolders:
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")
#         nH_star = np.loadtxt(path + '/nH_star.txt')
#     try:
#         mgsi = mgsi_from_nH(nH_star[0], nH_star[1])
#         mgsi_hyp.append(mgsi)
#     except IndexError as e:
#         # print(path, e, 'nH_star =', nH_star)
#         bad.append(path)
# print('n bad:', len(bad))


def filter_wds(df, exclude=None):
    """ exclude is a list of names to drop (duplicates? UV?) """
    # print(df_wd.columns)
    df.set_index('Name', inplace=True)

    # drop dupicates etc?
    if exclude:
        for name in exclude:
            df.drop(index=name, inplace=True)

    df.rename(columns={"log(Mg/Hx)": "Mg", "Mg (SS corrected)": "Mg_corr", "error log(Mg/Hx)": "Mg_err",
                       "log(Si/Hx)": "Si", "Si (SS corrected)": "Si_corr", "error log(Si/Hx)": "Si_err",
                       "log(Fe/Hx)": "Fe", "Fe (SS corrected)": "Fe_corr", "error log(Fe/Hx)": "Fe_err",
                       "log(Ca/Hx)": "Ca", "error log(Ca/Hx)": "Ca_err",
                       },
              inplace=True)

    # has some < upper limits, remove these
    df.Si = pd.to_numeric(df.Si, errors='coerce')
    df.Mg = pd.to_numeric(df.Mg, errors='coerce')
    df.Fe = pd.to_numeric(df.Fe, errors='coerce')
    df.Ca = pd.to_numeric(df.Ca, errors='coerce')

    # drop 0.0
    df[['Mg', 'Si', 'Fe', 'Ca']] = df[['Mg', 'Si', 'Fe', 'Ca']].replace(0.0000, np.nan)

    # add solar norm for comparison
    df['Mg_corr_norm'] = df.Mg_corr - p.mg_sol
    df['Si_corr_norm'] = df.Si_corr - p.si_sol
    df['Fe_corr_norm'] = df.Fe_corr - p.fe_sol
    # print('mean Mg pwd', np.nanmean(df_wd.Mg_corr_norm))

    df['Mg_norm'] = df.Mg - p.mg_sol
    df['Si_norm'] = df.Si - p.si_sol
    df['Fe_norm'] = df.Fe - p.si_sol
    df['Ca_norm'] = df.Ca - p.si_sol



    return df


""" plot scatter with errors """


def plot_scatter_MgX(df_hyp, df_wd, denominator='Si', c_wd=None, c_hyp=None, fig=None, ax=None, lims=(-4, 1), sun=True,
                     delta_lines=None, delta_pos=1, lws=None, pos_offset=0.1, pos_n=7,
                     labelsize=14, ticksize=10, legsize=10, save=False, fig_path=plotpx.fig_path):
    if c_wd is None:
        c_wd = [i / 255 for i in list(hex_to_rgb('#34013f'))]  # xkcd:dark violet
    if c_hyp is None:
        c_hyp = 'xkcd:kelly green'
    if delta_lines is None:
        delta_lines = [-0.7, -0.2, 0, 0.2, 0.7]
    if lws is None:
        lws = [0.3, 0.5, 0.9, 0.5, 0.3]

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(9, 9))
    ax.set_aspect('equal')

    y_hyp = df_hyp.Mg
    yerr_hyp = df_hyp.Mg_err
    y_wd = df_wd.Mg
    y_wd_corr = df_wd.Mg_corr
    yerr_wd = df_wd.Mg_err
    x_sol = eval('p.' + denominator.lower() + '_sol')
    x_hyp = eval('df_hyp.' + denominator)
    xerr_hyp = eval('df_hyp.' + denominator + '_err')
    x_wd = eval('df_wd.' + denominator)
    xerr_wd = eval('df_wd.' + denominator + '_err')
    xlabel = denominator + r'$_\star$/H$_\star$'
    ratiostr = 'Mg/' + denominator
    if denominator == 'Si':
        xcorr = df_wd.Si_corr
        delta_lines = [-0.7, -0.18, 0, 0.18, 0.7]  # corresponds to sio2 lower mantle (mgsi = 0.7) and opx-free (1.5)
    elif denominator == 'Fe':
        xcorr = df_wd.Fe_corr
        lims = (-10, -3.2)
    elif denominator == 'Ca':
        xcorr = [np.nan] * len(df_wd.Ca_err)
        lims = (-9, -3.2)
        delta_lines = [-0.7, 0, 0.7, 1, 1.3]
        lws = [0.5, 0.9, 0.5, 0.3, 0.3]
        delta_pos = 4
        pos_offset = 0.2
        pos_n = 6

    idx_core_rich = df_wd['Fragment Core Mass Fraction'] > 0.4
    # print('core_rich', idx_core_rich[idx_core_rich == True])
    # c_core = [q/255 for q in [0,191,255]] + [1]
    c_core = [q / 255 for q in [200, 64, 19]] + [1]

    # Hypatia stars (unnormalised from solar)
    ax.errorbar(x=x_hyp + x_sol, y=y_hyp + p.mg_sol, xerr=xerr_hyp, yerr=yerr_hyp, marker='.', lw=0,
                ecolor=c_hyp, c=c_hyp,
                elinewidth=0.5,  # markersize=10,
                capsize=None, alpha=0.2, zorder=10, label='Hypatia FGKM stars')

    # White Dwarfs (never normalised)
    if not np.isnan(xcorr).all():
        ax.errorbar(x=xcorr[~idx_core_rich], y=y_wd_corr[~idx_core_rich],
                    xerr=xerr_wd[~idx_core_rich], yerr=yerr_wd[~idx_core_rich],
                    marker='.', lw=0, c=c_wd,
                    elinewidth=1, markersize=15,
                    capsize=None, alpha=1, zorder=60, label='Polluted WDs, corrected')
        ax.errorbar(x=xcorr[idx_core_rich], y=y_wd_corr[idx_core_rich],
                    xerr=xerr_wd[idx_core_rich], yerr=yerr_wd[idx_core_rich],
                    lw=0, elinewidth=1, markersize=15, markeredgecolor=c_core, markerfacecolor=c_wd, alpha=1,
                    markeredgewidth=1.5, ecolor=c_wd, marker='.',
                    capsize=None, zorder=61, )

    ax.errorbar(x=x_wd[~idx_core_rich], y=y_wd[~idx_core_rich],
                xerr=xerr_wd[~idx_core_rich], yerr=yerr_wd[~idx_core_rich],
                marker='s', lw=0, c='k',
                elinewidth=1, markersize=8,
                capsize=None, alpha=0.3, zorder=50, label='Polluted WDs, uncorrected')

    ax.errorbar(x=x_wd[idx_core_rich], y=y_wd[idx_core_rich],
                xerr=xerr_wd[idx_core_rich], yerr=yerr_wd[idx_core_rich],
                lw=0, elinewidth=1, markersize=8, markeredgecolor=c_core, markerfacecolor=(0, 0, 0, 0.4),
                markeredgewidth=1.5, ecolor=(0, 0, 0, 0.4), marker='s',
                capsize=None, zorder=51, label='Polluted WDs, possibly core-rich')

    if sun:
        print('sun', eval('p.' + denominator.lower() + '_sol'), p.mg_sol)
        ax.scatter(eval('p.' + denominator.lower() + '_sol'), p.mg_sol, marker='*', s=80, c='xkcd:light yellow',
                   zorder=100, label='Sun')

    print('core rich:\n', df_wd[idx_core_rich].loc[:, ['Mg_corr', 'Si_corr', 'Mg', 'Si', 'Ca']], '\n')

    # print(df_wd.Si_corr_norm[~df_wd.Si_corr.isnull()])
    # print(df_wd.Mg_corr_norm[~df_wd.Mg_corr.isnull()])
    print('pwd corrected n =', df_wd.Si_corr.count(), df_wd.Mg_corr.count())

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel(r'Mg$_\star$/H$_\star$', fontsize=labelsize)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.legend(frameon=False, fontsize=legsize)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # 1:1 line
    for ii, (dl, lw) in enumerate(zip(delta_lines, lws)):
        x = np.linspace(lims[0], lims[1], endpoint=True, num=pos_n)
        y = x + dl
        ax.plot(x, y, c='k', ls=(1, (5, 10)), lw=lw, alpha=0.9, zorder=20)
        # print('x0, y0', x[0], y[0], 'log(mg/si)', y[0] - x[0])
        ms = ratio_from_nH(y[0], x[0])  # get mg/si of this ref. line
        ax.text(x[delta_pos] - (ii * pos_offset), y[delta_pos] - (ii * pos_offset), ratiostr + ' = {:0.1f}'.format(ms),
                fontsize=8, rotation=45, zorder=21,
                bbox=dict(boxstyle='square,pad=-0.1', fc='w', ec='none'))

    if save:
        fig.savefig(fig_path + 'mgsi_hypatia_error.pdf', bbox_inches='tight')
    return fig, ax


""" plot hist """


def plot_hist_MgX(df_hyp, df_wd, denominator='Si', c_wd=None, c_hyp=None, fig=None, ax=None, histrange=(0, 4),
                  labelsize=14, ticksize=10, legsize=10,
                  bins_wd=40, bins_hyp=50, ylim=None, lw_median=2, ls_median='--',
                  save=False, fig_path=plotpx.fig_path):
    # calculate hypatia ratio
    mgsi_hyp = ratio_from_nH(df_hyp.Mg + p.mg_sol, df_hyp.Si + p.si_sol)
    mgfe_hyp = ratio_from_nH(df_hyp.Mg + p.mg_sol, df_hyp.Fe + p.fe_sol)
    mgca_hyp = ratio_from_nH(df_hyp.Mg + p.mg_sol, df_hyp.Ca + p.ca_sol)
    # print('mean Mg hypatia', np.nanmean(df_hyp.Mg))

    # calculate wd ratio
    # molar ratius ignoring optical/UV problems
    # mgsi_wd = ratio_from_nH(df_wd.Mg_corr, df_wd.Si_corr)
    mgsi_wd = 10 ** df_wd['Mg/Si (SS corrected)']
    mgfe_wd = ratio_from_nH(df_wd.Mg_corr, df_wd.Fe_corr)

    # mgsi_wd_uncorr = 10 ** (df_wd['Mg/Si'][(df_wd['log(Mg/Hx)'] != 0.0) & (df_wd['log(Si/Hx)'] != 0.0)])
    # mgfe_wd_uncorr = 10 ** ((df_wd['Fe/Mg'][(df_wd['log(Mg/Hx)'] != 0.0) & (df_wd['log(Fe/Hx)'] != 0.0)]) ** -1)
    mgsi_wd_uncorr = ratio_from_nH(df_wd.Mg, df_wd.Si)
    mgfe_wd_uncorr = ratio_from_nH(df_wd.Mg, df_wd.Fe)
    mgca_wd_uncorr = ratio_from_nH(df_wd.Mg, df_wd.Ca)

    mgsi_wd_uncorr.replace(1.0000, np.nan, inplace=True)
    mgfe_wd_uncorr.replace(1.0000, np.nan, inplace=True)
    mgca_wd_uncorr.replace(1.0000, np.nan, inplace=True)
    mgca_wd_uncorr.replace(1.0000, np.nan, inplace=True)
    # print('mgsi wd uncorr\n', mgsi_wd_uncorr[~mgsi_wd_uncorr.isnull()])

    if c_wd is None:
        c_wd = [i / 255 for i in list(hex_to_rgb('#34013f'))]  # xkcd:dark violet
    if c_hyp is None:
        c_hyp = 'xkcd:kelly green'

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    if denominator == 'Si':
        ratio_hyp = mgsi_hyp
        ratio_wd = mgsi_wd
        ratio_uncorr = mgsi_wd_uncorr
        ratiostr = 'Mg/Si'
    elif denominator == 'Fe':
        ratio_hyp = mgfe_hyp
        ratio_wd = mgfe_wd
        ratio_uncorr = mgfe_wd_uncorr
        histrange = (0, 10)
        bins_wd = 40
        ratiostr = 'Mg/Fe'
    elif denominator == 'Ca':
        ratio_hyp = mgca_hyp
        ratio_wd = pd.Series([np.nan])
        ratio_uncorr = mgca_wd_uncorr
        histrange = (0, 60)
        ylim = (0, 0.25)
        bins_wd = 30
        ratiostr = 'Mg/Ca'

    # remove 0 values (undefined)

    ratio_hyp = ratio_hyp[ratio_hyp > 0.0]
    ratio_wd = ratio_wd[ratio_wd > 0.0]
    ratio_uncorr = ratio_uncorr[ratio_uncorr > 0.0]

    # plot hists and medians
    ax.axvline(np.nanmedian(ratio_hyp), lw=lw_median, c=c_hyp, alpha=1, ls=(4, (7, 3)))
    ax.hist(ratio_hyp, histtype='stepfilled', lw=1, fc=c_hyp, ec='k', alpha=0.6, range=histrange, density=True,
            bins=bins_hyp, label='Hypatia FGKM stars ($N$ = ' + '{:d})'.format(ratio_hyp.count()))
    print('Hypatia FGKM stars, 1 sigma percentiles:', np.nanpercentile(ratio_hyp, [16, 50, 84]))

    if not np.isnan(ratio_wd).all():
        ax.axvline(np.nanmedian(ratio_wd), lw=lw_median, c=c_wd, alpha=1, ls=(1, (7, 3)))
        ax.hist(ratio_wd, histtype='step', lw=1.5, ec=c_wd, hatch='\\\\', alpha=1, range=histrange, density=True,
                bins=12, label='Polluted WDs, corrected ($N$ = ' + '{:d})'.format(ratio_wd.count()))
        print('Polluted WDs, corrected, 1 sigma percentiles:', np.nanpercentile(ratio_wd, [16, 50, 84]))

    ax.axvline(np.nanmedian(ratio_uncorr), lw=lw_median, c='0.5', alpha=1, ls=(7, (7, 3)))
    ax.hist(ratio_uncorr, histtype='stepfilled', lw=1.5, fc='k', ec=c_wd, zorder=9, alpha=0.2, range=histrange,
            density=True,
            bins=bins_wd, label='Polluted WDs, uncorrected ($N$ = ' + '{:d})'.format(ratio_uncorr.count()))
    print('Polluted WDs, uncorrected, 1 sigma percentiles:', np.nanpercentile(ratio_uncorr, [16, 50, 84]))

    ax.legend(frameon=False, fontsize=legsize, loc='upper right')
    ax.set_xlabel(ratiostr, fontsize=labelsize)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.set_xlim(histrange)
    ax.set_ylim(ylim)
    ax.set_yticks([])
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    with pd.option_context('display.max_rows', None,
                           'display.max_columns', None,
                           'display.precision', 1,
                           ):
        # print(mgsi_wd_uncorr[~mgsi_wd_uncorr.isnull()])
        # print(mgca_wd_uncorr[~mgca_wd_uncorr.isnull()])
        print(mgfe_wd_uncorr[~mgfe_wd_uncorr.isnull()])
        print(df_wd[['Fe', 'Mg']][~mgfe_wd_uncorr.isnull()])

    if save:
        fig.savefig(fig_path + 'mgsi_hist.pdf', bbox_inches='tight')
    return fig, ax


""" load hypatia from csv """

df_hyp = pd.read_csv('/home/claire/Works/hypatia-25042023.tsv', delimiter='\t', header=0, index_col=None)
# print(df_hyp.columns)


""" load andy pwd """

df_wd = pd.read_csv('/home/claire/Works/andy_WD_master_spreadsheet.csv', delimiter=',', header=0, index_col=None,
                    na_values=['#VALUE!'], keep_default_na=True)
df_wd = filter_wds(df_wd, exclude=duplicates)

""" set up gridspec """

denominator = 'Si'

labelsize = 14
ticksize = 10
legsize = 10
fig = plt.figure(figsize=(9, 12))
gs = fig.add_gridspec(2, 1, height_ratios=(1, 4.5),
                      # bottom=0.1, top=0.9,
                      hspace=0.1)
ax_scat = fig.add_subplot(gs[1])
ax_hist = fig.add_subplot(gs[0])

fig, ax_scat = plot_scatter_MgX(df_hyp, df_wd, denominator=denominator, fig=fig, ax=ax_scat,
                                lims=(-8, -3.4), labelsize=labelsize, ticksize=ticksize, legsize=legsize)
fig, ax_hist = plot_hist_MgX(df_hyp, df_wd, denominator=denominator, fig=fig, ax=ax_hist,
                             histrange=(0, 3), labelsize=labelsize, ticksize=ticksize, legsize=legsize,
                             bins_wd=20, bins_hyp=70, ylim=(0, 3), lw_median=1, ls_median=(1, (5, 10)),
                             )

if denominator == 'Si':
    fname = 'wd_hypatia_mgsi'
elif denominator == 'Fe':
    fname = 'wd_hypatia_mgfe'
elif denominator == 'Ca':
    fname = 'wd_hypatia_mgca'
fig.savefig('/home/claire/Works/thesis/v0/Introduction/Figs/' + fname + '.pdf', bbox_inches='tight')

plt.show()

# todo: check if there is any correlation with core mass fraction (e.g. Si entering core)
