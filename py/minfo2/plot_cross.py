import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
from py.useful_and_bespoke import cornertext
import py.main as rw
import os
import pandas as pd
import numpy as np
import meltsfugacitydata as mfug
import perplexfugacitydata as pfug
import py.bulk_composition as bulk
from matplotlib import rc
import datetime
from matplotlib.ticker import FormatStrFormatter

today = str(datetime.date.today())
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

""" make cross plot of oxides / elements"""


def element_xplot(p_of_interest=1, T_of_interest=None, components=[], y_name='logfo2',
                  output_parent_path=fo2plt.output_parent_px, output_parent_path_melts=fo2plt.output_parent_mlt_earth,
                  fig=None, axes=[], fformat='.png', ax_histy=None,
                  ylim=(-11.5, -7), model='melts', ylabel=r'log($f_{{\rm O}_2}$)', xlim=None, xlabels=None,
                  labelsize=16, ticksize=14, save=True, fname=None, z_name=None, fig_path=fo2plt.figpath,
                  make_legend=True, verbose=False, exclude_names=[], exclude_silica=True, make_hist=False, **sc_kwargs):
    """ plot fo2 vs. wt% of some component at pressure of interest (in GPa)
    components can be bulk oxides or mineral phase proportion """

    if fname is None:
        fname = 'crossplot'
    if xlabels is None:
        xlabels = components

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
        axes[ii].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
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
                                except AttributeError:
                                    try:
                                        nH_star = np.loadtxt(output_parent_path_melts + name + '/nH_star.txt')
                                        x = nH_star[[ox[:2] for ox in dat.oxide_list].index(component[:2])]
                                    except FileNotFoundError:
                                        x = np.nan
                                    except IndexError as e:
                                        x = np.nan  # nHstar file is empty

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
                                    d = dat.get_phase_composition_dict(p_of_interest=p_of_interest,
                                                                       T_of_interest=T_of_interest,
                                                                       component='Fe2O3',
                                                                       phases=[phase],
                                                                       to_absolute_abundance=False,
                                                                       verbose=verbose)
                                    print(z_name, 'd', d)
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
                      density=True, orientation='horizontal', color='w', edgecolor='k', bins=15)
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
ylabel = r'$\Delta$QFM'
# model = 'melts'
model = 'perplex'
components = ['Mg/Si', 'Fe/Si', 'Al/Si', 'Ca/Si']
markersize = 40
labelsize = 16
legsize = 12
ticksize = 10
alpha = 0.3
vmin, vmax = 0, 0.6
mec = 'xkcd:midnight blue'
cmap = 'viridis'
fformat = '.pdf'
make_hist = True

ncols = len(components)
fig = plt.figure(figsize=(ncols * 3, 6))
gs = fig.add_gridspec(2, ncols + 1, width_ratios=[10] * ncols + [1], height_ratios=[10, 10],
                      left=0.1, right=0.9,
                      wspace=0.05, hspace=0.05)

for jj, p_of_interest in enumerate((1, ncols)):
    if model == 'melts':
        source = fo2plt.output_parent_mlt_earth
        # ylims = [(-11, -8.5), (-8, -5.3)]  # 1 GPa, 4 GPa - for logfO2 axis
        ylims = [(-2.5, 1.5), (-1.5, 2.5)]  # delta QFM
        mec = 'xkcd:midnight blue'
        title = 'pMELTS'
    elif model == 'perplex':
        source = fo2plt.output_parent_px
        # ylims = [(-13.5, -8.7), (-10, -6.5)]  # 1 GPa, 4 GPa - for logfO2 axis
        ylims = [(-2.5, 1.5), (-1.5, 2.5)]  # delta QFM
        mec = 'xkcd:brick red'
        title = 'Perple_X'
    output_parent_path = source + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/'

    axes = []  # create new ax list just for this row
    [axes.append(fig.add_subplot(gs[jj, ii])) for ii in range(ncols)]
    ax_histy = None
    if make_hist:
        ax_histy = fig.add_subplot(gs[jj, -1])  # for y histogram if making one

    fig, axs = element_xplot(p_of_interest=p_of_interest, T_of_interest=T_of_interest,
                             components=components, y_name=y_name, ylabel=ylabel,
                             output_parent_path=output_parent_path,
                             xlim=[(0.8, 1.65), (0.05, 0.2), (0, 0.18), (0.025, 0.13)],
                             ylim=ylims[jj],
                             edgecolors=mec, s=markersize, alpha=alpha,   # sc_kwargs
                             labelsize=labelsize, ticksize=ticksize,
                             save=False, fig=fig, axes=axes,
                             model=model, verbose=False,
                             z_name='X_Fe3_Opx', vmin=vmin, vmax=vmax, cmap=cmap,
                             ax_histy=ax_histy, make_hist=True,
                             exclude_silica=exclude_silica, fformat=fformat)

    # print('ylim', axs[0].get_ylim())

    axs[0] = cornertext(axs[0], str(p_of_interest) + ' GPa', size=labelsize, pos='top left')
    axs[1].set_xticks([0.075, 0.125, 0.175])  # Fe/Si
    if jj == 0:
        for ax in axs:
            ax.set_xticks([])
            ax.set_xlabel(None)

fig.suptitle(title, fontsize=labelsize, y=0.92)
fig.savefig(fo2plt.figpath + 'crossplot_elements_' + model + today + fformat, bbox_inches='tight')
# plt.show()














""" other phase abundances """

# if p_of_interest == 1:
#     ylim = (-14, -9)
#     phcomps = ['Ol', 'Opx', 'Cpx', 'Sp']
# elif p_of_interest == 4:
#     ylim = (-11.5, -7)
#     phcomps = ['Ol', 'Opx', 'Cpx', 'Gt']
# if not exclude_silica:
#     phcomps.append('q')
# fo2plt.element_xplot(p_of_interest=p_of_interest, components=phcomps,
#                      output_parent_path=output_parent_path,
#                      ylim=ylim, linec='k', labelsize=16, save=True, fname='crossplot_phases_' + str(p_of_interest),
#                      model=model, verbose=True,
#                      exclude_silica=exclude_silica)


# fo2plt.element_xplot(p_of_interest=p_of_interest, components=['MgO', 'SiO2', 'Al2O3', 'FeO', 'CaO'],
#                      xlim=[(30, 45), (43, 57), (1, 5), (3,10), (1,5)],
#                      y_name='X_Gt', ylabel='Gt abundance (wt%)',
#                      output_parent_path=output_parent_path,
#                      ylim=(0, 16),
#                      linec='k', labelsize=16, save=True, fname='crossplot_gt_' + str(p_of_interest),
#                      model=model, verbose=False,
#                      exclude_silica=exclude_silica)
#
# fo2plt.element_xplot(p_of_interest=p_of_interest, components=['MgO', 'SiO2', 'Al2O3', 'FeO', 'CaO'],
# xlim=[(30, 45), (43, 57), (1, 5), (3,10), (1,5)],
#                      y_name='X_Opx', ylabel='Opx abundance (wt%)',
#                      output_parent_path=output_parent_path,
#                      ylim=(0,80),
#                      linec='k', labelsize=16, save=True, fname='crossplot_opx_' + str(p_of_interest),
#                      model=model, verbose=False,
#                      exclude_silica=exclude_silica)


