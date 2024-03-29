import numpy as np
import pandas as pd
import os
import pathlib
import subprocess
# import py.perplexdata as px
import py.bulk_composition as bulk
from py import main as rw
from scipy import interpolate
import py.parameters as p
import perplexfugacitydata as pfug
import meltsfugacitydata as mfug
import matplotlib.pyplot as plt
from py.useful_and_bespoke import colorize, colourbar, iterable_not_string
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
import matplotlib.cm as cm

figpath = '/home/claire/Works/min-fo2/figs_scratch/'
output_parent_px = '/home/claire/Works/min-fo2/perplex_output/hypatia/'
output_parent_mlt_earth = '/home/claire/Works/min-fo2/alphamelts_output/earth-tea23/'
output_parent_mlt = '/home/claire/Works/min-fo2/alphamelts_output/'
plot_kwargs = {'labelsize': 16}
melt_phases = ['ctjL', 'dijL', 'enL'] + ['liquid']  # perplex + melts
c_phase_dict_stolper = {'Ol': 'tab:green', 'Opx': 'k', 'Cpx': 'tab:gray', 'Sp': 'tab:orange', 'Gt': 'tab:purple',
                        'q': 'tab:red', 'coe': 'tab:pink'}


# c_phase_dict_stolper = {'olivine': 'tab:green', 'orthopyroxene': 'k', 'clinopyroxene': 'tab:gray',
#                         'spinel': 'tab:orange', 'garnet': 'tab:purple',
#                         'quartz': 'tab:red', 'coesite': 'tab:pink'}


def filter_silica_sat(df):
    """ filter out cases saturated in quartz - these get to constant fo2 """
    if ('X_q' in df.columns) or ('X_coe' in df.columns) or ('X_stv' in df.columns):
        # print('   found q or coe')
        df.iloc[:, 2:] = np.nan  # everything except p, T
        print(df.head())
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
    #     dat = pfug.init_from_build(name, output_parent_path=output_parent_path)
    #     print('      Mg/Si =', dat.mgsi)

    return df


def fo2_xsection(name=None, df=None, output_parent_path=output_parent_px, fig=None, ax=None, xlabel=None, ylabel=None,
                 c_buf='k', legsize=12,
                 show_buffer=True, linec='k', lw=1, alpha=1, labelsize=16, save=True, fname=None, make_legend=True,
                 p_min=None, p_max=None, ymin=None, ymax=None, verbose=False, exclude_silica=True, **kwargs):
    """ plot isothermal cross section, p range in GPa """

    if fname is None:
        fname = name + '_fo2_xsection'
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    # read in results
    if df is None:
        df = pd.read_csv(output_parent_path + name + '/' + name + '_results.csv', sep='\t', index_col=0)
        print(output_parent_path + name + '/' + name + '_results.csv', 'exists')
        if verbose:
            print('\nloaded:\n', df.head(), '\nfrom', output_parent_path + name + '/' + name + '_results.csv')
        df = apply_filters(df, name, p_min, p_max)
        if exclude_silica:
            df = filter_silica_sat(df)

    if not df['P(bar)'].isnull().all():
        try:
            pressure = df['P(bar)'] * 1e-4
            T = df['T(K)'].unique()[0]  # for 2D
        except Exception as e:
            print('\nWARNING: loaded:\n', df.head(), '\nfrom', output_parent_path + name + '/' + name + '_results.csv')
            return fig, ax

        # plot fo2 columns
        fo2 = df['logfo2']
        # delta_qfm = df['delta_qfm']
        ax.plot(pressure, fo2, c=linec, lw=lw, alpha=alpha)
        # print('plotted', pressure[0], fo2[0])
        if show_buffer:
            try:
                fo2_qfm = df['logfo2_qfm']
                ax.plot(pressure, fo2_qfm, c=c_buf, lw=1, ls=(0, (10, 4)), label='FMQ')
            except KeyError:
                print(name, 'logfo2_qfm not calculated')

        # set legends and labels, clean up axes
        if make_legend:
            ax.legend(frameon=False, fontsize=legsize)
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
                    show_in_out=True, axes_all=None, make_legend=True, p_min=None, p_max=None, mode='line',
                    model='perplex', **kwargs):
    if fname is None:
        fname = name + '_phase_xsection'
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    df = pd.read_csv(output_parent_path + name + '/' + name + '_results.csv', sep='\t')
    df = apply_filters(df, name, p_min, p_max)
    pressure = df['P(bar)'] * 1e-4
    T = df['T(K)'].unique()[0]
    print(name, '\n', df.head())
    # plot all composition columns
    for col in df.columns:
        phase_name = col[2:]
        is_phase = col.startswith('X_') and (phase_name not in melt_phases) and ((df[col] != 0).any())
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
        ax.legend(frameon=False)
    ax.set_xlim(p_min, p_max)
    ax.set_ylim(0, ax.get_ylim()[1])

    if save:
        fig.savefig(figpath + fname + '.png')
    return fig, ax


def multicomp_xsection(output_parent_path=output_parent_px, fig=None, ax=None, save=True, hist_y=False,
                       ax_histy=None, model='perplex',
                       fname=None, cmap='summer', cmap_var=None, vmin=None, vmax=None, verbose=False, has_cbar=False,
                       cbar_label=None, exclude_names=[], bins=15, **kwargs):
    if fname is None:
        fname = model + '_fo2_variation'

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
            name = os.path.basename(sub)
            # print(sub + name + '_results.csv')
            # print(os.path.exists(sub + name + '_results.csv'))

            if (len(os.listdir(sub)) > 1) and (
            os.path.exists(sub + '/' + name + '_results.csv')):  # 1 if contains nH_star e.g.

                if name not in exclude_names:
                    if model == 'perplex':
                        dat = pfug.init_from_results(name, output_parent_path=output_parent_path)
                    elif model == 'melts':
                        dat = mfug.init_from_results(name, output_parent_path=output_parent_path, verbose=False,
                                                     **kwargs)
                    if dat is not None:
                        # get colouring
                        if cmap_var is not None:
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
                        once = False  # only draw buffer once
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


def stolper_subplot(name=None, output_parent_path=output_parent_px, p_min=None, p_max=None, title=None, fname=None,
                    save=True, fig=None, axes=None, exclude_silica=True, phase_mode='stacked', verbose=True,
                    T_of_interest=1373, **kwargs):
    if fname is None:
        fname = name + '_fo2_subplot'

    df = pd.read_csv(output_parent_path + name + '/' + name + '_results' + str(int(T_of_interest)) + '.csv', sep='\t', index_col=0)
    if not df['logfo2'].isnull().all():
        if verbose:
            print('\nloaded:\n', df.head(), '\nfrom', output_parent_path + name + '/' + name + '_results.csv')
        df = apply_filters(df, name, p_min, p_max)
        if exclude_silica:
            df = filter_silica_sat(df)

        if fig is None:
            fig, axes = plt.subplots(2, 1, figsize=(6, 8))
        fig, axes[0] = fo2_xsection(df=df, fig=fig, ax=axes[0], save=False, name=name,
                                    output_parent_path=output_parent_path, **kwargs)

        if phase_mode == 'stacked':
            fig, axes[1] = single_composition(df, phases=['Gt', 'Cpx', 'Opx', 'Ol', 'Plag'], fig=fig, ax=axes[1],
                                              comp_stacked=True, fig_path=figpath, save=False,
                                              show=False, legtitle=None, p_max=p_max,
                                              ylabel='Phase fraction (wt%)', xlabel='Pressure (GPa)',
                                              verbose=verbose, make_legend=True,
                                              leg_bbox_to_anchor=(1, 1), orientation='horizontal', scale=1,
                                              output_parent_path=output_parent_path,
                                              **kwargs)
        else:
            fig, axes[1] = phases_xsection(name=name, df=df, fig=fig, ax=axes[1], save=False, show_in_out=True,
                                           axes_all=axes, output_parent_path=output_parent_path,
                                           **kwargs)
        axes[0].set_xlabel('')
        # axes[0].set_xticks([])

    fig.suptitle(title)
    plt.subplots_adjust(hspace=0)
    if save:
        fig.savefig(figpath + fname + '.png')
    return fig, axes


def single_composition(df, phases=['Plag', 'Gt', 'Cpx', 'Opx', 'Ol'],
                       comp_stacked=True, fig_path=figpath, save=True, filename='min-mode',
                       show=False, legtitle=None, override_ax_arrow=False, ylabelpad=None, xlabelpad=None,
                       xlabel=None, ylabel=None, cmap='tab20', labelsize=16, p_max=None,
                       verbose=False, fig=None, ax=None, title=None, extension='.png', make_legend=True,
                       leg_bbox_to_anchor=(1, 1), legsize=14, orientation='horizontal', scale=100, **plot_kwargs):
    import matplotlib
    """ comp_stacked: plot cumulative modes
    p_max: maximum pressure to plot in GPa (e.g., to see UM and MTZ better"""

    # setup plot
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    ax.set_title(title, fontsize=labelsize)

    pressures = df['P(bar)'].to_numpy() * 1e-4
    x = pressures
    # p_max_orig = p_max
    if p_max is None:
        p_max = np.max(x)

    if orientation == 'horizontal':
        ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
        ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
        ax.set_xlim(1000e-4, p_max)
    elif orientation == 'vertical':
        ax.set_ylabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
        ax.set_xlabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
        ax.set_ylim(1000e-4, p_max)
        ax.invert_yaxis()

    # get colours from cmap
    cmap_vals = matplotlib.cm.get_cmap(cmap)
    n_phases = len(phases)
    colours = [cmap_vals(j) for j in np.arange(0, n_phases) / n_phases]

    for ii, phase in enumerate(phases):
        # loop over phases and either plot instantly (not stacked) or create stackplot matrix
        col = 'X_' + phase  # weight fraction of phase
        try:
            # dat.data[col].replace(0, np.nan, inplace=True)
            y = df[col].to_numpy(dtype=float) * scale
        except KeyError:
            print(col, 'not found in data')
            y = np.zeros(len(df), dtype=float)  # assign null value so phase order is retained
        if comp_stacked:
            # add to array for stackplot
            if ii == 0:
                y_stacked = y  # initialise
            else:
                y_stacked = np.vstack((y_stacked, y))
        else:
            # or, if not a stackplot, plot immediately
            if orientation == 'horizontal':
                ax.plot(x, y, c=colours[ii], label=phase, **plot_kwargs)
            elif orientation == 'vertical':
                ax.plot(y, x, c=colours[ii], label=phase, **plot_kwargs)
    # make stackplot at this point
    print('y stacked\n', y_stacked)
    y_stacked = np.nan_to_num(y_stacked)
    if comp_stacked:
        # print('y_stacked', np.shape(y_stacked))
        if orientation == 'horizontal':
            ax.stackplot(x, y_stacked, labels=phases, colors=colours)
        elif orientation == 'vertical':
            y_trans = np.cumsum(y_stacked, axis=0)
            # print('y_trans', np.shape(y_trans))
            for s, col in enumerate(phases):
                ax.fill_betweenx(x, y_trans[s, :], label=phase, fc=colours[s], zorder=-s)
            ax.set_xlim(0, 100)
        # ax.scatter(x, [50]*len(x), marker='.', s=20, c='k', zorder=100)  # show perple_x resolution

    if make_legend:
        leg = ax.legend(loc='upper left', bbox_to_anchor=leg_bbox_to_anchor, frameon=False, fontsize=legsize,
                        title=legtitle)
        if legtitle is not None:
            leg.get_title().set_fontsize(legsize)  # legend 'Title' fontsize
    # print('p_max', p_max)
    # if p_max_orig is None:
    #     p_max = np.max(x)
    # else:
    #     # assume not all pressures are given, annotate axis to indicate there's more and not an error
    #     ax = annotate_ax_continuation(ax, 'x')
    # ax.set_xlim(0, p_max)

    if save:
        plt.tight_layout()
        fig.savefig(fig_path + filename + extension, bbox_inches='tight')
    if show:
        plt.show()
    return fig, ax


def pop_hist(dirs=None, x_var=None, z_var=None, z_val=None, x_scale=1, z_scale=1, fname=None, figpath=figpath,
             save=True, make_legend=False, T_of_interest=1373, z_labels=None,
             exclude_names=[], exclude_silica=True, model='perplex', hilite=None, figsize=(7, 5), fig=None,
             ax=None, show_sd='corner', xlim=None, ylim=None,
             ls='-', cmap=None, c='k', vmin=None, vmax=None, xlabel=None, title=None, labelsize=16, legsize=12,
             nbins=20, ticksize=12, labelpad=None, x_shift=0, legloc=None, alpha=0.5,
             hist_kwargs={}, legtitle=None, verbose=True, **kwargs):
    """ z_var is the colour/linestyle coding, must be consistent across all runs in dir[n] """
    if cmap and (not vmin or not vmax):
        raise NotImplementedError('Cannot colour-code without explicit vmin and vmax')
    if cmap:
        # get normalised cmap
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        colours = cm.ScalarMappable(norm=norm, cmap=cmap)
        # print('vmin, vmax', vmin, vmax)
    if fname is None:
        fname = 'pop_hist'
    if xlabel is None:
        xlabel = x_var
    if legtitle is None:
        legtitle = z_var

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    for jj, opp in enumerate(dirs):
        spl = os.path.dirname(opp).split('/')[-1].split('_')
        X_ferric = None
        core_eff = None
        for sp in spl:
            if 'ferric' in sp:
                X_ferric = int(''.join(filter(str.isdigit, sp)))
            if 'coreeff' in sp:
                core_eff = int(''.join(filter(str.isdigit, sp))) / 100

        # get directory names in folder
        subfolders = rw.get_run_dirs(output_path=opp)
        once = True
        x = []
        n_qtz = 0
        if subfolders:  # nonzero
            # print('subfolders', subfolders)
            for ii, sub in enumerate(subfolders):  # loop over stars
                if len(os.listdir(sub)) > 1:  # 1 if contains nH_star e.g.
                    name = os.path.basename(sub)
                    if name not in exclude_names:
                        if model == 'perplex':
                            dat = pfug.init_from_results(name, X_ferric=X_ferric,  # core_eff=core_eff,
                                                         output_parent_path=opp, T_iso=T_of_interest,
                                                         load_results_csv=True, verbose=False, **kwargs)
                            if exclude_silica and (dat is not None) and (T_of_interest != 1373) and \
                                    ~set(['X_Ol', 'X_Opx']).issubset(dat.data.columns):
                                # didn't load composition in some temps
                                dat1373 = pfug.init_from_results(name, X_ferric=X_ferric,  # core_eff=core_eff,
                                                                 output_parent_path=opp, T_iso=1373,
                                                                 load_results_csv=True, verbose=False, **kwargs)
                                if 'X_q' in dat1373.data.columns:  # don't print quartz saturation cases
                                    n_qtz += 1
                                    if verbose:
                                        print('dropping case with quartz:', dat.name)
                                    continue
                        elif model == 'melts':
                            # print('T_of_interest', T_of_interest)
                            dat = mfug.init_from_results(name, X_ferric=X_ferric, core_eff=core_eff,
                                                         output_parent_path=opp, T_final=T_of_interest,
                                                         load_results_csv=True, verbose=False, **kwargs)
                        if dat is not None:
                            if exclude_silica:
                                if 'X_q' in dat.data.columns:  # don't print quartz saturation cases
                                    n_qtz += 1
                                    if verbose:
                                        print('dropping case with quartz:', dat.name)
                                    continue
                            try:
                                ans = eval('dat.' + x_var) * x_scale
                                # if ans > 1.5:
                                #     print(dat.name)
                                #     print(dat.data.columns)
                                #     print(dat.data.tail())
                                x.append(eval('dat.' + x_var) * x_scale - x_shift)
                                # x.append(10**(eval('dat.' + x_var) * x_scale - x_shift))
                            except TypeError as e:
                                print(e)
                                print('excluding', name, eval('dat.' + x_var))
                            except AttributeError as e:
                                print(e)
                                print('excluding', name, ': missing attribute', x_var)

                            # get colour (do this here because final dat could be None
                            if not z_val:
                                z = eval('dat.' + z_var) * z_scale  # only need 1 (last in this case)
                            else:
                                z = z_val
            print(os.path.basename(os.path.normpath(opp)), 'n_runs =', len(x), '......................(dropped', n_qtz, 'cases with silica phase)')
            # n_qtz : 30 - 164 in melts from 70-99 ceff (i.e. iron stabilises
            # " "   : 31 - 162
            # 75

            if cmap:
                c = colours.to_rgba(z)  # i.e., override input c
                # print('z', z)
                try:
                    c_alpha = list(c)
                    c_alpha[3] = alpha
                except IndexError as e:
                    print('c', c)
                    print('z', z)
                    print('failed to get list of colours, no runs..?')
                    raise e
            else:
                c_alpha = c
            if iterable_not_string(ls):
                linestyle = ls[jj]
            else:
                linestyle = ls

            if exclude_silica:
                sd = np.nanstd(x)  # because otherwise some fo2 results set to nan
            else:
                sd = np.std(x)

            if z_labels is None:
                z_label = '{:.0f}%'.format(z)  # + r' ($\sigma =$ ' + '{:.2f})'.format(sd),
            else:
                z_label = z_labels[jj]
            if (hilite is not None and (jj == hilite)) or hilite is None:
                print('z_label', z_label, 'c_alpha', c_alpha)
                ax.hist(x, ec=c, fc=c_alpha, ls=linestyle, histtype='stepfilled',
                        label=z_label, bins=nbins, density=True, **hist_kwargs)
            else:
                ax.hist(x, color=c, alpha=0.5, ls=linestyle, histtype='step',
                        label=z_label,  bins=nbins, density=True, **hist_kwargs)

            if show_sd and ((hilite is not None and (jj == hilite)) or hilite is None):  # turn off for non-highlighted if multi
                if show_sd == 'corner':
                    xt, yt = 0.99, 0.9  # put text in upper right corner
                    transform = ax.transAxes
                    ha, va = 'right', 'top'
                elif show_sd == 'leftcorner':
                    xt, yt = 0.01, 0.9  # put text in upper right corner
                    transform = ax.transAxes
                    ha, va = 'left', 'top'
                elif show_sd == 'above':
                    xt, yt = np.median(x), 0.9 * ylim[1]  # put text above hist centre, ylim0 = 0 for hist
                    # xt, yt = np.median(x), 0.85 * ylim[1]   # put text above hist centre, ylim0 = 0 for hist
                    transform = ax.transData
                    ha, va = 'left', 'top'
                ax.text(xt, yt, r'$\sigma =$ ' + '{:.2f}'.format(sd), c=c_alpha[:3], alpha=1, fontsize=legsize,
                        ha=ha, va=va, transform=transform)
                print(os.path.basename(os.path.normpath(opp)), 'sd =', sd)


    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=labelpad)
    ax.set_title(title, fontsize=labelsize)
    ax.tick_params(axis='both', labelsize=ticksize)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if make_legend:
        leg = ax.legend(title=legtitle, fontsize=legsize, frameon=False, loc=legloc)
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
