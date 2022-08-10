import numpy as np
import pandas as pd
import os
import pathlib
import subprocess
import perplexdata as px
import bulk_composition as bulk
import main as rw
from scipy import interpolate
import parameters as p
import perplexfugacity as pf
import matplotlib.pyplot as plt
from useful_and_bespoke import colorize, colourbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

figpath = '/home/claire/Works/min-fo2/figs_scratch/'
output_parent_default = '/home/claire/Works/min-fo2/perplex_output/'
plot_kwargs = {'labelsize': 16}
melt_phases = ['ctjL', 'dijL', 'enL']
c_phase_dict_stolper = {'Ol': 'tab:green', 'Opx': 'k', 'Cpx': 'tab:gray', 'Sp': 'tab:orange', 'Gt': 'tab:purple',
                        'q': 'tab:red', 'coe': 'tab:pink'}

def remove_bad(df):
    """ remove rows with no stable perple_x solution """
    s = df[[col for col in df.columns if col.startswith('X_')]].sum(axis=1)
    idx = s == 0
    return df.loc[~idx]


def check_subsolidus(df):
    for phase in melt_phases:
        if (phase in df.columns) and ((df[phase] != 0).any()):
            return False  # False if there are bad ones
    return True


def check_nonmonotonic(df):
    """ find non-monotonically increasing fo2 """
    idx = [True] + [df['logfo2'].iloc[ii] < df['logfo2'].iloc[ii - 1] for ii in range(1, len(df))]
    tmp = df[idx]
    if len(tmp) > 0:
        print('bad idx', idx)
        return False
    return True


def fo2_xsection(name, output_parent_path=output_parent_default, fig=None, ax=None, xlabel=None, ylabel=None,
                 show_buffer=True, linec='k', lw=1, alpha=1, labelsize=16, save=True, fname=None, make_legend=True,
                 p_min=None, p_max=None, ymin=None, ymax=None, **kwargs):
    """ plot isothermal cross section, p range in GPa """

    if fname is None:
        fname = name + '_fo2_xsection'
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    # read in results
    df = pd.read_csv(output_parent_path + name + '/' + name + '_results.csv', sep='\t')

    # checks for bad data
    if not check_subsolidus(df):
        raise Exception('ERROR: melt phases present in', name)
    df = remove_bad(df)  # drop solutionless rows

    if not check_nonmonotonic(df):
        print(name, 'non-monotonic!!!!!!!!!!!')

    # extract and subset pressure
    if p_min is None:
        p_min = df['P(bar)'].min() * 1e-4
    if p_max is None:
        p_max = df['P(bar)'].max() * 1e-4
    df = df[(df['P(bar)'] * 1e-4 >= p_min) & (df['P(bar)'] * 1e-4 <= p_max)]
    pressure = df['P(bar)'] * 1e-4
    T = df['T(K)'].unique()[0]

    # plot fo2 columns
    fo2 = df['logfo2']
    delta_qfm = df['delta_qfm']
    ax.plot(pressure, fo2, c=linec, lw=lw, alpha=alpha, label='log($f$O$_2$)')
    if show_buffer:
        fo2_qfm = df['logfo2_qfm']
        ax.plot(pressure, fo2_qfm, c='k', lw=1, ls=(0, (10, 4)), label='FMQ')

    # check high-fo2 scenarios
    if (delta_qfm > 0).any():
        print(name, ': fo2 > QFM')
        dat = pf.init_from_build(name, output_parent_path=output_parent_path)
        print('      Mg/Si =', dat.mgsi)

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


def phases_xsection(name, output_parent_path=output_parent_default, fig=None, ax=None, xlabel=None, ylabel=None,
                    linec='k', c_dict=c_phase_dict_stolper, lw=1, labelsize=16, save=True, fname=None, y_annot=45,
                    show_in_out=True, axes_all=None, make_legend=True, p_min=None, p_max=None, **kwargs):

    if fname is None:
        fname = name + '_phase_xsection'
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    df = pd.read_csv(output_parent_path + name + '/' + name + '_results.csv', sep='\t')

    if not check_subsolidus(df):
        raise Exception('ERROR: melt phases present in', name)
    df = remove_bad(df)  # drop solutionless rows

    # extract and subset pressure
    if p_min is None:
        p_min = df['P(bar)'].min() * 1e-4
    if p_max is None:
        p_max = df['P(bar)'].max() * 1e-4
    df = df[(df['P(bar)'] * 1e-4 >= p_min) & (df['P(bar)'] * 1e-4 <= p_max)]
    pressure = df['P(bar)'] * 1e-4
    T = df['T(K)'].unique()[0]

    # plot all composition columns
    for col in df.columns:
        if col.startswith('X_') and (col[2:] not in melt_phases) and ((df[col] != 0).any()):
            if c_dict is not None:
                c = c_dict[col[2:]]
            y = df[col]
            ax.plot(pressure, y, c=c, lw=lw, label=col[2:])

            if show_in_out:
                y.replace(0, np.nan, inplace=True)
                idx_in = y.first_valid_index()
                idx_out = y.last_valid_index()

                if axes_all is None:
                    axes_all = [ax]

                for axx in axes_all:  # i.e. if this is part of a subplot, extend these vlines to full col of axes
                    if idx_in > pressure.index.min():  # don't plot if stable from start
                        axx.axvline(x=pressure[idx_in], c=c, lw=1, ls='--')
                        ax.text(pressure[idx_in], y_annot, col[2:] + '-in', c=c, va='center', ha='left', rotation=90)  # label (but only comp subplot)
                    if idx_out < pressure.index.max():  # don't plot if stable to end
                        axx.axvline(x=pressure[idx_out], c=c, lw=1, ls='--')
                        ax.text(pressure[idx_out], y_annot, col[2:] + '-out', c=c, va='center',
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


def multicomp_xsection(output_parent_path=output_parent_default, fig=None, ax=None, save=True, hist_y=False, ax_histy=None,
                       fname=None, cmap='summer', cmap_var=None, vmin=None, vmax=None, verbose=False, has_cbar=False,
                       cbar_label=None, exclude_names=None, **kwargs):

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
                if exclude_names and (name not in exclude_names):
                    # get colouring
                    if cmap_var is not None:
                        dat = pf.init_from_build(name, output_parent_path=output_parent_path)
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
        ax_histy.hist(y, #range=ylim,
                      density=True, orientation='horizontal', color='w', edgecolor='k')

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
        cbar = colourbar(vector=np.linspace(vmin, vmax), ax=ax, vmin=vmin, vmax=vmax, label=cbar_label, labelsize=14, ticksize=12,
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


def element_xplot():
    """ plot fo2 vs. wt% of some element """


output_parent_path = output_parent_default + 'hypatia88Fe/'
names = [
    'Stolper'
          # '1M_88Ceff_HIP58237_999K',
    # 'dmm'
]
exclude_names = ['1M_88Ceff_HIP58237_999K']

# # like figure 3 from Stolper+2020
# for name in names:
#     stolper_subplot(name, output_parent_path=output_parent_path, show_buffer=True, save=True, p_min=1.1, lw=2, **plot_kwargs)

# cmap_var, vmin, vmax = 'mg_number', 86, 94
cmap_var, vmin, vmax, cbar_label = 'mgsi', 0.7, 1.3, 'Mg/Si'
multicomp_xsection(output_parent_path=output_parent_path, cmap='spring', cmap_var=cmap_var, vmin=vmin, vmax=vmax, has_cbar=True, cbar_label=cbar_label,
                   save=True, fname=None, alpha=0.9, lw=0.5, p_min=1.1, ymin=-13, ymax=-6,
                   hist_y=True, exclude_names=exclude_names, **plot_kwargs)

plt.show()
