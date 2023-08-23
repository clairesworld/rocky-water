import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
from py.useful_and_bespoke import cornertext, dark_background
import py.main as rw
import os
import pandas as pd
import numpy as np
import meltsfugacitydata as mfug
import perplexfugacitydata as pfug
import py.bulk_composition as bulk
from matplotlib import rc
import datetime
from matplotlib.ticker import FormatStrFormatter, PercentFormatter, FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmcrameri

today = str(datetime.date.today())
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

""" make cross plot of oxides / elements"""


def element_xplot(p_of_interest=1, T_of_interest=None, components=[], y_name='logfo2',
                  output_parent_path=fo2plt.output_parent_px, output_parent_path_melts=fo2plt.output_parent_mlt_earth,
                  fig=None, axes=[], fformat='.png', ax_histy=None, nbins=None, x_str_format=None,
                  ylim=(-11.5, -7), model='melts', ylabel=r'log($f_{{\rm O}_2}$)', xlim=None, xlabels=None,
                  labelsize=16, ticksize=14, save=True, fname=None, z_name=None, fig_path=fo2plt.figpath,
                  make_legend=True, verbose=False, exclude_names=[], exclude_silica=True, make_hist=False, **sc_kwargs):
    """ plot fo2 vs. wt% of some component at pressure of interest (in GPa)
    components can be bulk oxides or mineral phase proportion """

    if fname is None:
        fname = 'crossplot'
    if xlabels is None:
        xlabels = components
    if x_str_format is None:
        x_str_format = ['%.1f']*len(components)

    # setup gridspec
    ncols = len(components)
    if fig is None:
        fig = plt.figure(figsize=(ncols * 4, 4))
        gs = fig.add_gridspec(1, ncols + 1, width_ratios=[10] * ncols + [1],
                              left=0.1, right=0.9,
                              wspace=0.05)
        [axes.append(fig.add_subplot(gs[0, ii])) for ii in range(ncols)]

    for ii in range(ncols):
        axes[ii].set_xlabel(xlabels[ii], fontsize=labelsize)
        axes[ii].tick_params(axis='both', labelsize=ticksize)
        axes[ii].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        axes[ii].xaxis.set_major_formatter(FormatStrFormatter(x_str_format[ii]))
        axes[ii].set_ylim(ylim)
        if xlim:
            axes[ii].set_xlim(xlim[ii])
        if ii > 0:
            axes[ii].set_yticks([])
    axes[0].set_ylabel(ylabel, fontsize=labelsize)

    # get directory names in folder
    subfolders = rw.get_run_dirs(output_path=output_parent_path)
    if subfolders:  # nonzero
        once = True
        # loop over runs
        ys = []
        for ii, sub in enumerate(subfolders):
            name = os.path.basename(sub)
            results_file = output_parent_path + name + '/' + name + '_results' + str(int(T_of_interest)) + '.csv'
            if (len(os.listdir(sub)) > 0) and os.path.exists(results_file):  # > 1 if contains nH_star e.g.
                if name not in exclude_names:

                    if model == 'melts':
                        dat = mfug.init_from_results(name=name, output_parent_path=output_parent_path,
                                                     load_results_csv=True, verbose=False)
                    elif model == 'perplex':
                        dat = pfug.init_from_results(name=name, output_parent_path=output_parent_path,
                                                     load_results_csv=True, verbose=False)

                    if dat is not None:
                        if exclude_silica:
                            if 'X_q' in dat.data.columns:  # don't print quartz saturation cases
                                if verbose:
                                    print('dropping case with quartz:', dat.name)
                                continue

                        # get pressure row
                        try:
                            idx = dat.data['P(bar)'].sub(np.round(p_of_interest) * 1e4).abs().idxmin()
                            row = dat.data.iloc[idx]
                        except (TypeError, IndexError) as e:
                            # pressure not found
                            print('idx', idx, name)
                            print(row)
                            continue

                        # loop over components to check (subplots)
                        for ii, (ax, component) in enumerate(zip(axes, components)):
                            # search for component in oxides
                            if component in dat.wt_oxides.keys():
                                x = dat.wt_oxides[component]
                            # search for component in phase comp
                            elif 'X_' + component in row.index:
                                x = row['X_' + component]
                            # maybe ratio with H
                            elif '/H' in component:
                                try:
                                    x = dat.nH_star[[ox[:2] for ox in dat.oxide_list].index(component[:2])]
                                except AttributeError as e:
                                    # nH_star not loaded yet
                                    try:
                                        nH_file = '/home/claire/Works/min-fo2/alphamelts_output/hypatia_local2/hypatia_88coreeff_3ferric_ext/' + name + '/nH_star.txt'
                                        nH_star = np.loadtxt(nH_file)
                                        x = nH_star[[ox[:2] for ox in dat.oxide_list].index(component[:2])]
                                        # print('loaded nH_star, x', x)
                                    except FileNotFoundError:
                                        print('missing', nH_file)
                                        x = np.nan
                                    except IndexError as e:
                                        x = np.nan  # nHstar file is empty

                            elif component == 'O_total':
                                # total refractory O in wt%, add up
                                x = bulk.total_refractory_O(dat.wt_oxides)

                            # maybe it's an element ratio
                            elif '/' in component:
                                x = bulk.get_element_ratio(component, dat.wt_oxides)
                            else:
                                if verbose:
                                    print(name, ':', component, 'not found in', dat.wt_oxides.keys(), 'or', row.index)
                                x = np.nan
                            try:
                                y = row[y_name]
                                # print(y)
                            except KeyError as e:
                                print(name, 'keyerror', e)
                                print(row, '\n\n')
                                y = np.nan
                            if z_name:
                                if z_name in dat.wt_oxides.keys():
                                    c = dat.wt_oxides[component]
                                # search for component in phase comp
                                elif z_name in row.index:
                                    c = row[z_name]
                                elif 'X_Fe3' in z_name:
                                    phase = z_name.split('_')[2]
                                    # perplex haven't loaded these in results
                                    if (p_of_interest == 4) and model == 'perplex':
                                        pprime = 3.9  # correct for fact u didnt run px at 4 gpa lol
                                    else:
                                        pprime = p_of_interest
                                    d = dat.get_phase_composition_dict(p_of_interest=pprime,
                                                                       T_of_interest=T_of_interest,
                                                                       component='Fe2O3',
                                                                       phases=[phase],
                                                                       to_absolute_abundance=False,
                                                                       verbose=verbose)
                                    # print(z_name, 'd', d)
                                    c = d[phase]
                                else:
                                    raise NotImplementedError(z_name, 'not implemented for scatter colouring')
                                # print('c', c, name)
                                sc_kwargs.update({'c': c})
                                # todo: need to check for vmin vmac cmap etc
                            ax.scatter(x, y, **sc_kwargs)
                            # print(dat.name, '(x, y) =', x, y)
                            # if once:  # for histogram, don't repeat run for all columns
                            if ii == 0:
                                ys.append(y)
                                T = dat.data['T(K)'].unique()[0]  # isothermal
                    else:
                        print('  ...problem initialising data object for', name, '(returned None)')
            elif verbose:
                print(sub, 'is empty')
        once = False

    # make hist
    if make_hist:
        if ax_histy is None:
            ax_histy = fig.add_subplot(gs[0, -1], sharey=axes[-1])
        ax_histy.hist(ys,  range=ylim,
                      density=True, orientation='horizontal', color='0.5', edgecolor='k', linewidth=0.4, bins=nbins)
        ax_histy.set_xticks([])
        ax_histy.set_yticks([])
        ax_histy.spines.bottom.set_visible(False)
        ax_histy.spines.right.set_visible(False)
        ax_histy.spines.top.set_visible(False)
        ax_histy.tick_params(axis="both", labelbottom=False, labelleft=False)

    print('num points:', len(ys))

    # plt.suptitle(str(p_of_interest) + ' GPa, ' + str(T) + ' K', fontsize=labelsize)

    if save:
        fig.savefig(fig_path + fname + fformat)
    return fig, axes





""" 
plot fo2 vs. element ratios (cols) for two pressures (rows) 
colour by Fe3+ content in opx
"""

core_eff = 88
Xf = 3
T_of_interest = 1373
exclude_silica = True
y_name = 'delta_qfm'
ylabel = r'$f_{{\rm O}_2}$ ($\Delta$FMQ)'
model = 'melts'
# model = 'perplex'
# components = ['Mg/Si', 'Fe/Si', 'Al/Si', 'Ca/Si']
# xlims = [(0.8, 1.65), (0.05, 0.2), (0, 0.18), (0.025, 0.13)]
components = ['Mg/Si', 'Al/Si', 'Fe/H', 'O_total']
pressures = [4]
# xlims = [(0.75, 1.75), (0, 0.19), (-5.4, -3.9), (41.7, 47.3)]
xlims = [(0.75, 1.75), (0, 0.19), (-5.4, -3.9), (43.1, 47.3)]
x_labels = ['Mg/Si', 'Al/Si', r'[Fe/H]$_{\star}$', r'Total O$_{\rm rock}$ (wt.\%)']
x_str_format = ['%.1f', '%.2f', '%.1f', '%.0f']
markersize = 40
labelsize = 16
legsize = 12
ticksize = 10
alpha = 0.3
vmin, vmax = 0, 0.8
mec = 'xkcd:midnight blue'
cmap = cmcrameri.cm.hawaii
fformat = '.png'
make_hist = True

ncols = len(components)
fig = plt.figure(figsize=(ncols * 2.5, len(pressures) * 2.5))
gs = fig.add_gridspec(len(pressures), ncols + 1, width_ratios=[10] * ncols + [1], height_ratios=[10] * len(pressures),
                      left=0.1, right=0.9,
                      wspace=0.05, hspace=0.05)

for jj, p_of_interest in enumerate(pressures):
    if model == 'melts':
        source = fo2plt.output_parent_mlt_earth
        # ylims = [(-11, -8.5), (-8, -5.3)]  # 1 GPa, 4 GPa - for logfO2 axis
        ylims = [(-1.5, 2.5)]  # delta QFM
        # yticks = [-1, 0, 1, 2]
        mec = 'xkcd:midnight blue'
        title = r'$f_{{\rm O}_2}$ distribution with pMELTS'
    elif model == 'perplex':
        source = fo2plt.output_parent_px
        # ylims = [(-13.5, -8.7), (-10, -6.5)]  # 1 GPa, 4 GPa - for logfO2 axis
        ylims = [(-4.4, 1.1)]  # delta QFM
        mec = 'xkcd:midnight blue'
        title = r'$f_{{\rm O}_2}$ distribution with Perple_X'

    output_parent_path = source + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/'

    axes = []  # create new ax list just for this row
    [axes.append(fig.add_subplot(gs[jj, ii])) for ii in range(ncols)]
    ax_histy = None
    if make_hist:
        ax_histy = fig.add_subplot(gs[jj, -1])  # for y histogram if making one

    fig, axs = element_xplot(p_of_interest=p_of_interest, T_of_interest=T_of_interest,
                             components=components, y_name=y_name, ylabel=ylabel,
                             output_parent_path=output_parent_path,
                             xlim=xlims, xlabels=x_labels, x_str_format=x_str_format,
                             ylim=ylims[jj],
                             edgecolors=mec, s=markersize, alpha=alpha,   # sc_kwargs
                             labelsize=labelsize, ticksize=ticksize,
                             save=False, fig=fig, axes=axes,
                             model=model, verbose=False,
                             z_name='X_Fe3_Opx', vmin=vmin, vmax=vmax, cmap=cmap,
                             ax_histy=ax_histy, make_hist=True, nbins=35,
                             exclude_silica=exclude_silica, fformat=fformat)

    # print('ylim', axs[0].get_ylim())
    if model == 'melts':
        plt.setp(axs[0].get_yticklabels()[0], visible=False)
        plt.setp(axs[0].get_yticklabels()[-1], visible=False)

    axs[0] = cornertext(axs[0], str(p_of_interest) + ' GPa', size=labelsize, pos='top left', c='xkcd:off white')
    axs[1].set_xticks([0.05, 0.1, 0.15])  # Al/Si
    # axs[1].set_xticks([0.075, 0.125, 0.175])  # Fe/Si
    if jj == 0:
    #     for ax in axs:
    #         ax.set_xticks([])
    #         ax.set_xlabel(None)

        # make colorbar top row
        dum = np.linspace(vmin, vmax)
        im = axs[-2].scatter(dum, dum, c=dum, cmap=cmap, s=0, vmin=vmin, vmax=vmax)
        cax = inset_axes(
            axs[-2],
            width="50%",  # width: 50% of parent_bbox width
            height="5%",  # height: 5%
            loc="lower right",
            bbox_to_anchor=(1, 1.05, 1, 1),
            bbox_transform=axs[-2].transAxes,
            borderpad=0,
        )


        def percent_tick_fmt(x, pos):
            return "{:.1f}".format(x) + r"\%"

        cbar = fig.colorbar(im, cax=cax, orientation="horizontal", format=FuncFormatter(percent_tick_fmt),
                            ticks=[vmin, vmax])

        cbar.ax.set_xlabel(r'[Fe$_2$O$_3$]$^{\rm opx}$', fontsize=ticksize+1, labelpad=-1)
        cbar.ax.xaxis.set_label_position('top')
        # cax.text(-0.5, 0.5, r'Fe$_2$O$_3$ in opx', fontsize=ticksize, transform=cax.transAxes, )
        cax.xaxis.set_ticks_position("top")
        cbar.ax.tick_params(axis="x", labelsize=ticksize-2, bottom=False, labelbottom=False, top=True, labeltop=True)
        # ax.yaxis.set_major_formatter(PercentFormatter())



fig.suptitle(title, fontsize=labelsize, y=0.99, c='xkcd:off white')

fig, *axs = dark_background(fig, list(axs) + [cbar.ax, ax_histy])
fig.savefig(fo2plt.figpath + 'crossplot_xtra_dark_' + model + today + fformat, bbox_inches='tight', dpi=300)
# plt.show()






