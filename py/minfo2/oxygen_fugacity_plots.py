import numpy as np
import pandas as pd
import os
import pathlib
import subprocess
import py.perplexdata as px
import py.bulk_composition as bulk
import py.main as rw
from scipy import interpolate
import py.parameters as p
import perplexfugacitydata as pf
import matplotlib.pyplot as plt
from py.useful_and_bespoke import colorize, colourbar, iterable_not_string
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
import matplotlib.cm as cm

figpath = '/home/claire/Works/min-fo2/figs_scratch/'
output_parent_px = '/home/claire/Works/min-fo2/perplex_output/'
plot_kwargs = {'labelsize': 16}
melt_phases = ['ctjL', 'dijL', 'enL'] + ['liquid']  # perplex + melts
c_phase_dict_stolper = {'Ol': 'tab:green', 'Opx': 'k', 'Cpx': 'tab:gray', 'Sp': 'tab:orange', 'Gt': 'tab:purple',
                        'q': 'tab:red', 'coe': 'tab:pink'}

c_phase_dict_stolper = {'olivine': 'tab:green', 'orthopyroxene': 'k', 'clinopyroxene': 'tab:gray',
                        'spinel': 'tab:orange', 'garnet': 'tab:purple',
                        'quartz': 'tab:red', 'coesite': 'tab:pink'}


def filter_silica_sat(df):
    """ filter out cases saturated in quartz - these get to constant fo2 """
    if 'X_q' in df.columns:
        # print('found q')
        df.iloc[:, 2:] = np.nan  # everything except p, T
        # print(df.head())
    return df

def remove_no_solution(df):
    """ remove rows with no stable perple_x solution """
    px_phase_cols = [col for col in df.columns if col.startswith('X_')]
    if len(px_phase_cols) > 0:
        s = df[px_phase_cols].sum(axis=1)
        idx = s == 0
        return df.loc[~idx]
    else:
        return df


def check_subsolidus(df):
    for phase in melt_phases:
        if (phase in df.columns) and ((df[phase] != 0).any()):
            return False  # False if there are bad ones
    return True


def check_nonmonotonic(df):
    """ find non-monotonically increasing fo2 """
    return ~df['logfo2'].is_monotonic_increasing


def apply_filters(df, name, p_min=None, p_max=None, output_parent_path=None):
    """ do all checks, p in GPa """
    try:
        tmp = df['P(bar)']
    except KeyError:
        # accidentally saved with external editor as ,-delimited
        df = pd.read_csv(output_parent_path + name + '/' + name + '_results.csv', sep=',')

    # checks for bad data
    if not check_subsolidus(df):
        raise Exception('ERROR: melt phases present in', name)
    df = remove_no_solution(df)  # drop solutionless rows

    if not check_nonmonotonic(df):
        print(name, 'non-monotonic')

    # extract and subset pressure
    flag = False
    if p_min is None:
        p_min = df['P(bar)'].min() * 1e-4
    else:
        flag = True
    if p_max is None:
        p_max = df['P(bar)'].max() * 1e-4
    else:
        flag = True
    if flag:
        df = df[(df['P(bar)'] * 1e-4 >= p_min) & (df['P(bar)'] * 1e-4 <= p_max)]

    # # check high-fo2 scenarios
    # if (df['delta_qfm'] > 0).any():
    #     print(name, ': fo2 > QFM')
    #     dat = pf.init_from_build(name, output_parent_path=output_parent_path)
    #     print('      Mg/Si =', dat.mgsi)

    return df


def fo2_xsection(name, output_parent_path=output_parent_px, fig=None, ax=None, xlabel=None, ylabel=None,
                 show_buffer=True, linec='k', lw=1, alpha=1, labelsize=16, save=True, fname=None, make_legend=True,
                 p_min=None, p_max=None, ymin=None, ymax=None, verbose=False, exclude_silica=True, **kwargs):
    """ plot isothermal cross section, p range in GPa """

    if fname is None:
        fname = name + '_fo2_xsection'
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    # read in results
    df = pd.read_csv(output_parent_path + name + '/' + name + '_results.csv', sep='\t', index_col=0)
    if verbose:
        print('\nloaded:\n', df.head(), '\nfrom', output_parent_path + name + '/' + name + '_results.csv')
    df = apply_filters(df, name, p_min, p_max)
    if exclude_silica:
        df = filter_silica_sat(df)

    try:
        pressure = df['P(bar)'] * 1e-4
        T = df['T(K)'].unique()[0]  # for 2D
    except Exception as e:
        print('\nloaded:\n', df.head(), '\nfrom', output_parent_path + name + '/' + name + '_results.csv')
        raise e

    # plot fo2 columns
    fo2 = df['logfo2']
    delta_qfm = df['delta_qfm']
    ax.plot(pressure, fo2, c=linec, lw=lw, alpha=alpha, label='log($f$O$_2$)')
    if show_buffer:
        fo2_qfm = df['logfo2_qfm']
        ax.plot(pressure, fo2_qfm, c='k', lw=1, ls=(0, (10, 4)), label='FMQ')

    # set legends and labels, clean up axes
    if make_legend:
        ax.legend()
    if xlabel is None:
        xlabel = 'Pressure (GPa)\nat ' + str(int(T)) + ' K'
    if ylabel is None:
        ylabel = 'log($f$O$_2$)'
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_xlim(p_min, p_max)
    ax.set_ylim(ymin, ymax)

    if save:
        fig.savefig(figpath + fname + '.png')
    return fig, ax


def phases_xsection(name, output_parent_path=output_parent_px, fig=None, ax=None, xlabel=None, ylabel=None,
                    linec='k', c_dict=c_phase_dict_stolper, lw=1, labelsize=16, save=True, fname=None, y_annot=45,
                    show_in_out=True, axes_all=None, make_legend=True, p_min=None, p_max=None,
                    model='perplex', **kwargs):
    if fname is None:
        fname = name + '_phase_xsection'
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    df = pd.read_csv(output_parent_path + name + '/' + name + '_results.csv', sep='\t')
    df = apply_filters(df, name, p_min, p_max)
    pressure = df['P(bar)'] * 1e-4
    T = df['T(K)'].unique()[0]

    # plot all composition columns
    for col in df.columns:
        if model == 'perplex':
            phase_name = col[2:]
            is_phase = col.startswith('X_') and (phase_name not in melt_phases) and ((df[col] != 0).any())
        elif model == 'melts':
            phase_name = col[:-2]
            is_phase = col.endswith('_0') and (phase_name not in melt_phases) and ((df[col] != 0).any())
        if is_phase:
            if c_dict is not None:
                c = c_dict[phase_name]
            y = df[col]
            ax.plot(pressure, y, c=c, lw=lw, label=phase_name)

            if show_in_out:
                y.replace(0, np.nan, inplace=True)
                idx_in = y.first_valid_index()
                idx_out = y.last_valid_index()

                if axes_all is None:
                    axes_all = [ax]

                for axx in axes_all:  # i.e. if this is part of a subplot, extend these vlines to full col of axes
                    if idx_in > pressure.index.min():  # don't plot if stable from start
                        axx.axvline(x=pressure[idx_in], c=c, lw=1, ls='--')
                        ax.text(pressure[idx_in], y_annot, phase_name + '-in', c=c, va='center', ha='left',
                                rotation=90)  # label (but only comp subplot)
                    if idx_out < pressure.index.max():  # don't plot if stable to end
                        axx.axvline(x=pressure[idx_out], c=c, lw=1, ls='--')
                        ax.text(pressure[idx_out], y_annot, phase_name + '-out', c=c, va='center',
                                ha='right', rotation=90)  # label (but only comp subplot)

    # make legends and labels, clean up axes
    if xlabel is None:
        xlabel = 'Pressure (GPa)\nat ' + str(int(T)) + ' K'
    if ylabel is None:
        ylabel = 'Phase fraction (wt%)'
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    if make_legend:
        ax.legend()
    ax.set_xlim(p_min, p_max)
    ax.set_ylim(0, ax.get_ylim()[1])

    if save:
        fig.savefig(figpath + fname + '.png')
    return fig, ax


def multicomp_xsection(output_parent_path=output_parent_px, fig=None, ax=None, save=True, hist_y=False,
                       ax_histy=None,
                       fname=None, cmap='summer', cmap_var=None, vmin=None, vmax=None, verbose=False, has_cbar=False,
                       cbar_label=None, exclude_names=[], bins=15, **kwargs):
    if fname is None:
        fname = 'fo2_variation'

    # setup fig
    if fig is None:
        fig = plt.figure(figsize=(6, 4))
    if hist_y:
        gs = fig.add_gridspec(1, 2, width_ratios=[10, 1],
                              left=0.1, right=0.9,
                              wspace=0.05)
        if ax is None:
            ax = fig.add_subplot(gs[0, 0])
        if ax_histy is None:
            ax_histy = fig.add_subplot(gs[0, 1], sharey=ax)
            ax_histy.spines.bottom.set_visible(False)
            ax_histy.spines.right.set_visible(False)
            ax_histy.spines.top.set_visible(False)
            ax_histy.tick_params(axis="both", labelbottom=False, labelleft=False)
            ax_histy.set_xticks([])
    elif ax is None:
        ax = fig.add_subplot(111)

    # setup colours
    if cmap_var is not None:
        if ((vmin is None) or (vmax is None)):
            raise NotImplementedError('ERROR: cannot colourise from variable without vmin, vmax')

    # get directory names in folder
    subfolders = rw.get_run_dirs(output_path=output_parent_path)
    once = True
    if subfolders:  # nonzero
        for ii, sub in enumerate(subfolders):
            if len(os.listdir(sub)) > 1:  # 1 if contains nH_star e.g.
                name = os.path.basename(sub)
                if name not in exclude_names:
                    # get colouring
                    if cmap_var is not None:
                        dat = pf.init_from_results(name, output_parent_path=output_parent_path)
                        try:
                            z = eval('dat.' + cmap_var)
                        except AttributeError as e:
                            print('cmap_var not valid')
                            raise e
                        linec = colorize(z, cmap=cmap, vmin=vmin, vmax=vmax)[0]  # override input linec
                        # print(name, cmap_var, z)

                    fig, ax = fo2_xsection(name=name, output_parent_path=output_parent_path, fig=fig, ax=ax,
                                           save=False, make_legend=False, verbose=verbose, linec=linec,
                                           show_buffer=once, **kwargs)
                    once = False  # only draw once
            elif verbose:
                print(sub, 'is empty')

    if hist_y:
        y = []

        # get line data
        lines = ax.get_lines()
        for line in lines:
            y.append(line.get_ydata()[-1])  # get rightmost value for hist plotting
        ax_histy.hist(y,  # range=ylim,
                      density=True, orientation='horizontal', color='w', edgecolor='k', bins=bins)

    if has_cbar:
        if cbar_label is None:
            cbar_label = cmap_var

        axins = inset_axes(ax,  # here using axis of the lowest plot
                           width="100%",  # width = 5% of parent_bbox width
                           height="5%",  # height : 340% good for a (4x4) Grid
                           loc='lower left',
                           bbox_to_anchor=(0, 1.05, 1, 1),
                           bbox_transform=ax.transAxes,
                           borderpad=0,
                           )
        cbar = colourbar(vector=np.linspace(vmin, vmax), ax=ax, vmin=vmin, vmax=vmax, label=cbar_label, labelsize=14,
                         ticksize=12,
                         labelpad=10, loc='top', cax=axins,
                         rot=None, cmap=cmap, pad=0.05)

    if save:
        fig.savefig(figpath + fname + '.png', bbox_inches='tight')


def stolper_subplot(name=None, fname=None, save=True, fig=None, axes=None, **kwargs):
    if fname is None:
        fname = name + '_fo2_subplot'

    if fig is None:
        fig, axes = plt.subplots(2, 1, figsize=(6, 8))
    fig, axes[0] = fo2_xsection(fig=fig, ax=axes[0], save=False, name=name, **kwargs)
    fig, axes[1] = phases_xsection(name=name, fig=fig, ax=axes[1], save=False, show_in_out=True, axes_all=axes,
                                   **kwargs)
    axes[0].set_xlabel('')
    axes[0].set_xticks([])

    axes[0].set_title(name)
    plt.subplots_adjust(hspace=0)
    if save:
        fig.savefig(figpath + fname + '.png')


def element_xplot(p_of_interest=1, components=[], output_parent_path=output_parent_px, fig=None, axes=None,
                  xlabel=None, ylabel=None, ylim=(-11.5, -7),
                  linec='k', c_dict=c_phase_dict_stolper, lw=1, labelsize=16, save=True, fname=None,
                  make_legend=True, verbose=False, exclude_names=[], exclude_silica=True, **kwargs):
    """ plot fo2 vs. wt% of some component at pressure of interest (in GPa)
    components can be bulk oxides or mineral phase proportion """

    if fname is None:
        fname = 'crossplot'

    # setup gridspec
    ncols = len(components)
    if fig is None:
        fig = plt.figure(figsize=(ncols * 4, 4))
    gs = fig.add_gridspec(1, ncols + 1, width_ratios=[10] * ncols + [1],
                          left=0.1, right=0.9,
                          wspace=0.05)
    axes = []
    for ii in range(ncols):
        axes.append(fig.add_subplot(gs[0, ii]))
        axes[ii].set_xlabel(components[ii], fontsize=labelsize)
        axes[ii].set_ylim(ylim)
        if ii > 0:
            axes[ii].set_yticks([])
    axes[0].set_ylabel('log($fO_2$)', fontsize=labelsize)

    # get directory names in folder
    subfolders = rw.get_run_dirs(output_path=output_parent_path)
    if subfolders:  # nonzero
        once = True
        # loop over runs
        ys = []
        for ii, sub in enumerate(subfolders):
            if len(os.listdir(sub)) > 1:  # 1 if contains nH_star e.g.
                name = os.path.basename(sub)
                if name not in exclude_names:
                    df = pd.read_csv(output_parent_path + name + '/' + name + '_results.csv', sep='\t')
                    df = apply_filters(df, name)
                    if exclude_silica:
                        df = filter_silica_sat(df)
                    idx = df['P(bar)'].sub(p_of_interest * 1e4).abs().idxmin()
                    row = df.iloc[idx]
                    d = pf.read_dict_from_build(name=name, output_parent_path=output_parent_path)

                    # loop over components to check (subplots)
                    for ii, (ax, component) in enumerate(zip(axes, components)):

                        # search for component in oxides
                        if component in d['wt_oxides'].keys():
                            x = d['wt_oxides'][component]
                        elif 'X_' + component in row.index:
                            x = row['X_' + component]
                        else:
                            if verbose:
                                print(name, ':', component, 'not found in', d['wt_oxides'].keys(), 'or', row.index)
                            x = np.nan
                        y = row['logfo2']
                        ax.scatter(x, y, c=linec, s=5)
                        if once:
                            ys.append(y)
                            T = df['T(K)'].unique()[0]  # isothermal

            elif verbose:
                print(sub, 'is empty')
        once = False

    # make hist
    ax_histy = fig.add_subplot(gs[0, -1], sharey=axes[-1])
    ax_histy.hist(ys,  # range=ylim,
                  density=True, orientation='horizontal', color='w', edgecolor='k', bins=15)
    ax_histy.set_xticks([])
    ax_histy.set_yticks([])
    ax_histy.spines.bottom.set_visible(False)
    ax_histy.spines.right.set_visible(False)
    ax_histy.spines.top.set_visible(False)
    ax_histy.tick_params(axis="both", labelbottom=False, labelleft=False)

    plt.suptitle(str(p_of_interest) + ' GPa, ' + str(T) + ' K', fontsize=labelsize)

    if save:
        fig.savefig(figpath + fname + '.png')
    return fig, axes


def compare_pop_hist(dirs, x_var, z_var=None, x_scale=1, z_scale=1, fname=None, figpath=figpath, save=True,
                     exclude_names=[], exclude_silica=True,
                     ls='-', cmap=None, c='k', vmin=None, vmax=None, xlabel=None, title=None, labelsize=16, legsize=12, bins=20,
                     hist_kwargs={}, legtitle=None, **kwargs):
    """ z_var is the colour/linestyle coding, must be consistent across all runs in dir[n] """
    if cmap and (not vmin or not vmax):
        raise NotImplementedError('Cannot colour-code without explicit vmin and vmax')
    if cmap:
        # get normalised cmap
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        colours = cm.ScalarMappable(norm=norm, cmap=cmap)
    if fname is None:
        fname = 'compare_pop_hist'
    if xlabel is None:
        xlabel = x_var
    if legtitle is None:
        legtitle = z_var

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    for jj, opp in enumerate(dirs):
        spl = os.path.dirname(opp).split('/')[-1].split('_')
        X_ferric = None
        for sp in spl:
            if 'ferric' in sp:
                X_ferric = int(''.join(filter(str.isdigit, sp)))

        # get directory names in folder
        subfolders = rw.get_run_dirs(output_path=opp)
        once = True
        x = []
        if subfolders:  # nonzero
            for ii, sub in enumerate(subfolders):
                if len(os.listdir(sub)) > 1:  # 1 if contains nH_star e.g.
                    name = os.path.basename(sub)
                    if name not in exclude_names:
                        dat = pf.init_from_results(name, X_ferric=X_ferric, output_parent_path=opp)
                        if exclude_silica:
                            dat.data = filter_silica_sat(dat.data)
                        dat.read_fo2_results()
                        x.append(eval('dat.' + x_var) * x_scale)
            print(opp, 'n_runs =', len(x))

            z = eval('dat.' + z_var) * z_scale  # only need 1 (last in this case)
            if cmap:
                c = colours.to_rgba(z)  # i.e., override input c
            else:
                c = c
            if iterable_not_string(ls):
                linestyle = ls[jj]
            else:
                linestyle = ls

            if exclude_silica:
                sd = np.nanstd(x)  # because otherwise some fo2 results set to nan
            else:
                sd = np.std(x)
            ax.hist(x, color=c, ls=linestyle, histtype='step', label=str(z) + r' ($\sigma =$ ' + '{:.2f})'.format(sd),
                    bins=bins, density=True, **hist_kwargs)

    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_title(title, fontsize=labelsize)
    leg = ax.legend(title=legtitle, fontsize=legsize, frameon=False)
    leg.get_title().set_fontsize(legsize)  # legend 'Title' fontsize

    if save:
        fig.savefig(figpath + fname + '.png')
    return fig, ax



""" checking weird cases? """
# names_nonmono = ['1M_88Ceff_2MASS07324421+3350061_999K', '1M_88Ceff_2MASS18484459+4426041_999K',
# '1M_88Ceff_HIP60644_999K', '1M_88Ceff_HIP12048_999K', '1M_88Ceff_HIP110813_999K', '1M_88Ceff_HIP30114_999K',
# '1M_88Ceff_HIP13291_999K', '1M_88Ceff_HIP99825_999K', '1M_88Ceff_HIP31246_999K', '1M_88Ceff_HIP90593_999K',
# '1M_88Ceff_HIP18387_999K', '1M_88Ceff_2MASS19171110+4627286_999K', '1M_88Ceff_HIP64792_999K',
# '1M_88Ceff_2MASS19535587+4647370_999K', '1M_88Ceff_HIP106353_999K', '1M_88Ceff_HIP68469_999K',
# '1M_88Ceff_2MASS19304273+4643361_999K', '1M_88Ceff_HIP2611_999K', '1M_88Ceff_HIP79219_999K',
# '1M_88Ceff_HIP10085_999K', '1M_88Ceff_HIP1499_999K', '1M_88Ceff_2MASS18452371+4417428_999K',
# '1M_88Ceff_HIP114855_999K', '1M_88Ceff_HIP34801_999K', '1M_88Ceff_HIP58237_999K']
# for name in names_nonmono:
#     stolper_subplot(name, output_parent_path=output_parent_path, show_buffer=True, save=True, p_min=1.1, lw=2, **plot_kwargs)


plt.show()
