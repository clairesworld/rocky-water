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


def fit_with_outliers(x, y):
    from sklearn.linear_model import HuberRegressor
    from sklearn.preprocessing import StandardScaler

    # # Huber Regressor
    # # standardize
    # x_scaler, y_scaler = StandardScaler(), StandardScaler()
    # x_train = x_scaler.fit_transform(np.array(x)[..., None])
    # y_train = y_scaler.fit_transform(np.array(y)[..., None])
    #
    # # fit model
    # try:
    #     mdl = HuberRegressor(epsilon=1)
    #     m = mdl.fit(x_train, y_train.ravel())
    #     coef = [m.coef_[0], m.intercept_]
    # except ValueError:
    #     # HuberRegregressor failed, probably bad fit anyways
    #     coef = np.polyfit(x, np.array(y), deg=1, cov=False)

    coef = np.polyfit(np.array(x), np.array(y), deg=1, cov=False)

    o = np.poly1d(coef)
    xp = np.linspace(np.min(x), np.max(x), 100)
    yp = o(xp)

    print('fitted to x, y:', np.shape(x), np.shape(y))

    # do some predictions
    # yp = y_scaler.inverse_transform(
    #     mdl.predict(x_scaler.transform(xp[..., None]))
    # )

    return coef, xp, yp


def element_xplot(p_of_interest=1, T_of_interest=None, components=[], y_name='logfo2',
                  output_parent_path=fo2plt.output_parent_px, output_parent_path_melts=fo2plt.output_parent_mlt_earth,
                  fig=None, axes=[], fformat='.png', ax_histy=None, nbins=None, x_str_format=None,
                  ylim=(-11.5, -7), model='melts', ylabel=r'log($f_{{\rm O}_2}$)', xlim=None, xlabels=None,
                  labelsize=16, ticksize=14, save=True, fname=None, z_name=None, fig_path=fo2plt.figpath,
                  make_legend=True, verbose=False, exclude_names=[], exclude_silica=True, make_hist=False,
                  fit=True, **sc_kwargs):
    """ plot fo2 vs. wt% of some component at pressure of interest (in GPa)
    components can be bulk oxides or mineral phase proportion """

    if fname is None:
        fname = 'crossplot'
    if xlabels is None:
        xlabels = components
    if x_str_format is None:
        x_str_format = ['%.1f'] * len(components)

    # setup gridspec
    ncols = len(components)
    if fig is None:
        fig = plt.figure(figsize=(ncols * 4, 4))
        gs = fig.add_gridspec(1, ncols + 1, width_ratios=[10] * ncols + [1],
                              left=0.1, right=0.9,
                              wspace=0.05)
        [axes.append(fig.add_subplot(gs[0, ii])) for ii in range(ncols)]

    all_x = []
    all_y = []
    all_names = []
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
        all_x.append([])
        all_y.append([])
        all_names.append([])
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
                            elif component in row.index:
                                x = row[component]
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

                            elif 'X_Fe3' in component:
                                # maybe not loaded
                                if (p_of_interest == 4) and model == 'perplex':
                                    pprime = 3.9  # correct for fact u didnt run px at 4 gpa lol
                                else:
                                    pprime = p_of_interest
                                x = dat.get_phase_composition_dict(p_of_interest=pprime,
                                                                   T_of_interest=T_of_interest,
                                                                   component='Fe2O3',
                                                                   phases=[component.split('_')[-1]],
                                                                   to_absolute_abundance=False,
                                                                   verbose=verbose)[component.split('_')[-1]]
                            else:
                                if verbose:
                                    print(name, ':', component, 'not found in', dat.wt_oxides.keys(), 'or', row.index)
                                x = np.nan
                            try:
                                y = row[y_name]
                                # print(y)
                            except KeyError as e:
                                if '/' in y_name:
                                    y = bulk.get_element_ratio(y_name, dat.wt_oxides)
                                else:
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
                                    try:
                                        c = d[phase]
                                    except TypeError as e:
                                        print(dat.name, model)
                                        raise e
                                else:
                                    raise NotImplementedError(z_name, 'not implemented for scatter colouring')
                                # print('c', c, name)
                                sc_kwargs.update({'c': c})
                                # todo: need to check for vmin vmac cmap etc
                            ax.scatter(x, y, **sc_kwargs)
                            all_x[ii].append(x)
                            all_y[ii].append(y)
                            all_names[ii].append(name)

                            # # sanity
                            # if (p_of_interest == 4) and (model == 'melts'):
                            #     if name in (fit_exclude_4GPa_melts):
                            #         ax.scatter(x, y, c='g', s=markersize)

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

    if fit:
        fit_c = 'xkcd:light green blue'  # 'xkcd:purple/pink' #'xkcd:pistachio'
        for ii, (xii, yii, namesii) in enumerate(zip(all_x, all_y, all_names)):
            if components[ii] == 'Mg/Si':
                print(components[ii], 'p =', p_of_interest)
                print('xii =')
                print(repr(xii))
                print('yii =')
                print(repr(yii))
                print('names')
                print(repr(namesii))

                # filter
                if model == 'melts' and p_of_interest == 1:
                    exclude = fit_exclude_1GPa_melts
                elif model == 'melts' and p_of_interest == 4:
                    exclude = fit_exclude_4GPa_melts
                elif model == 'perplex' and p_of_interest == 1:
                    exclude = fit_exclude_1GPa_perplex
                elif model == 'perplex' and p_of_interest == 4:
                    exclude = fit_exclude_4GPa_perplex

                xfit, yfit, nfit = [], [], []
                for x0, y0, n0 in zip(xii, yii, namesii):
                    if (n0 not in exclude) and (x0 <= xmax_fit):
                        xfit.append(x0)
                        yfit.append(y0)
                        nfit.append(n0)
                    else:
                        print('excluding', n0, '(x={:.2f}, y={:.2f})'.format(x0, y0))
                        # axes[ii].scatter(x0, y0, c=(0, 0, 0, 0), edgecolors='b', s=markersize)
                        axes[ii].scatter(x0, y0, c=(0, 0, 0, 0.5), marker='x', s=markersize)

                coef, xp, yp = fit_with_outliers(np.array(xfit), np.array(yfit))
                print(components[ii], 'fit:', coef)

                # sanity
                axes[ii].plot(xp, yp, ls='--', c=fit_c, lw=1.5)
                axes[ii] = cornertext(axes[ii], 'm = {:.2f}'.format(coef[0]), pos='bottom right', c='k')

    # make hist
    if make_hist:
        if ax_histy is None:
            ax_histy = fig.add_subplot(gs[0, -1], sharey=axes[-1])
        ax_histy.hist(ys, range=ylim,
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
        fig.savefig(fig_path + fname + fformat, dpi=300)
    return fig, axes


def xplot_p_grid(model, components, pressures):
    ncols = len(components)
    nrows = len(pressures)
    fig = plt.figure(figsize=(ncols * 2.5, nrows * 2.5))
    gs = fig.add_gridspec(nrows, ncols + 1, width_ratios=[10] * ncols + [1], height_ratios=[10] * nrows,
                          left=0.1, right=0.9,
                          wspace=0.05, hspace=0.05)

    for jj, p_of_interest in enumerate(pressures):
        if model == 'melts':
            source = fo2plt.output_parent_mlt_earth
            # ylims = [(-11, -8.5), (-8, -5.3)]  # 1 GPa, 4 GPa - for logfO2 axis
            ylims = [(-2.5, 1.5), (-1.5, 2.5)]  # delta QFM
            # yticks = [-1, 0, 1, 2]
            mec = 'xkcd:midnight blue'
            title = r'pMELTS'
        elif model == 'perplex':
            source = fo2plt.output_parent_px
            # ylims = [(-13.5, -8.7), (-10, -6.5)]  # 1 GPa, 4 GPa - for logfO2 axis
            ylims = [(-5.1, 0.5), (-4.4, 1.1)]  # delta QFM
            mec = 'xkcd:midnight blue'
            title = r'Perple_X'

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
                                 edgecolors=mec, s=markersize, alpha=alpha,  # sc_kwargs
                                 labelsize=labelsize, ticksize=ticksize,
                                 save=False, fig=fig, axes=axes,
                                 model=model, verbose=False,
                                 z_name='X_Fe3_Opx', vmin=vmin, vmax=vmax, cmap=cmap,
                                 ax_histy=ax_histy, make_hist=make_hist, nbins=35,
                                 exclude_silica=exclude_silica, fformat=fformat)

        # print('ylim', axs[0].get_ylim())
        if model == 'melts':
            plt.setp(axs[0].get_yticklabels()[0], visible=False)
            plt.setp(axs[0].get_yticklabels()[-1], visible=False)

        axs[0] = cornertext(axs[0], str(p_of_interest) + ' GPa', size=labelsize, pos='top left')

        if jj < nrows - 1:
            for ax in axs:
                ax.set_xticks([])
                ax.set_xlabel(None)

        if jj == 0:
            # make colorbar top row
            dum = np.linspace(vmin, vmax)
            im = axs[0].scatter(dum, dum, c=dum, cmap=cmap, s=0, vmin=vmin, vmax=vmax)
            cax = inset_axes(
                axs[0],
                width="50%",  # width: 50% of parent_bbox width
                height="5%",  # height: 5%
                loc="lower right",
                bbox_to_anchor=(1, 1.05, 1, 1),
                bbox_transform=axs[0].transAxes,
                borderpad=0,
            )

            def percent_tick_fmt(x, pos):
                return "{:.1f}".format(x) + r"\%"

            cbar = fig.colorbar(im, cax=cax, orientation="horizontal", format=FuncFormatter(percent_tick_fmt),
                                ticks=[vmin, vmax])

            cbar.ax.set_xlabel(r'[Fe$_2$O$_3$]$^{\rm opx}$', fontsize=ticksize + 1, labelpad=-1)
            cbar.ax.xaxis.set_label_position('top')
            # cax.text(-0.5, 0.5, r'Fe$_2$O$_3$ in opx', fontsize=ticksize, transform=cax.transAxes, )
            cax.xaxis.set_ticks_position("top")
            cbar.ax.tick_params(axis="x", labelsize=ticksize - 2, bottom=False, labelbottom=False, top=True,
                                labeltop=True)
            # ax.yaxis.set_major_formatter(PercentFormatter())

    fig.suptitle(title, fontsize=labelsize, y=0.93)
    fig.savefig(fo2plt.figpath + 'crossplot_mgsifit_' + model + today + fformat, bbox_inches='tight')
    print('saved figure')

    # plt.show()
    # fig, *axs = dark_background(fig, axs)
    # fig.savefig(fo2plt.figpath + 'crossplot_xtra_dark_' + model + today + fformat, bbox_inches='tight')

    return fig


""" 
plot fo2 vs. element ratios (cols) for two pressures (rows) 
colour by Fe3+ content in opx
"""

core_eff = 88
Xf = 3
T_of_interest = 1373
exclude_silica = True
y_name = 'delta_qfm'
ylabel = r'$\Delta$FMQ'
components = ['Mg/Si']
pressures = (1, 4)
xlims = [(0.75, 1.75)]
x_labels = ['Mg/Si']
x_str_format = ['%.1f']
markersize = 40
labelsize = 16
legsize = 12
ticksize = 10
alpha = 0.3
vmin, vmax = 0, 0.8
mec = 'xkcd:midnight blue'
cmap = cmcrameri.cm.hawaii
fformat = '.png'
make_hist = False

fit = True
xmax_fit = 10  # 1.5  # Mg/Si max

# 0.05 chisqr diff
# fit_exclude_1GPa_melts = np.array(['1M_88Ceff_2MASS19394601-2544539_999K_3,0fer'], dtype='<U48')

# 0.01 chisqr diff
fit_exclude_1GPa_melts = np.array(['1M_88Ceff_2MASS19451246+5040203_999K_3,0fer',
                                   '1M_88Ceff_2MASS04215269+5749018_999K_3,0fer',
                                   '1M_88Ceff_2MASS19064452+4705535_999K_3,0fer',
                                   '1M_88Ceff_2MASS19004386+4349519_999K_3,0fer',
                                   '1M_88Ceff_HIP63510_999K_3,0fer',
                                   '1M_88Ceff_2MASS19403794+4053558_999K_3,0fer',
                                   '1M_88Ceff_2MASS19283611+3938152_999K_3,0fer',
                                   '1M_88Ceff_2MASS19250148+4915322_999K_3,0fer',
                                   '1M_88Ceff_2MASS19422610+3844086_999K_3,0fer',
                                   '1M_88Ceff_2MASS19394601-2544539_999K_3,0fer',
                                   '1M_88Ceff_2MASS19470525+4145199_999K_3,0fer',
                                   '1M_88Ceff_HIP48331_999K_3,0fer', '1M_88Ceff_HIP15578_999K_3,0fer',
                                   '1M_88Ceff_HIP77838_999K_3,0fer'], dtype='<U48')

# 0.05 chisqr diff
fit_exclude_4GPa_melts = np.array(['1M_88Ceff_2MASS19394601-2544539_999K_3,0fer',
                                   '1M_88Ceff_2MASS19191922+5035104_999K_3,0fer',
                                   '1M_88Ceff_HIP28941_999K_3,0fer', '1M_88Ceff_HIP61028_999K_3,0fer'],
                                  dtype='<U48')

# 0.01 chisqr diff
fit_exclude_1GPa_perplex = np.array(['1M_88Ceff_2MASS04215269+5749018_999K_3,0fer',
                                     '1M_88Ceff_2MASS19064452+4705535_999K_3,0fer',
                                     '1M_88Ceff_2MASS19411202+4033237_999K_3,0fer',
                                     '1M_88Ceff_2MASS19422610+3844086_999K_3,0fer',
                                     '1M_88Ceff_2MASS19231995+3811036_999K_3,0fer',
                                     '1M_88Ceff_2MASS19443633+5005449_999K_3,0fer',
                                     '1M_88Ceff_2MASS19191922+5035104_999K_3,0fer',
                                     '1M_88Ceff_HIP77838_999K_3,0fer', '1M_88Ceff_HIP17054_999K_3,0fer',
                                     '1M_88Ceff_HIP28941_999K_3,0fer', '1M_88Ceff_HIP61028_999K_3,0fer',
                                     '1M_88Ceff_2MASS19490218+4650354_999K_3,0fer',
                                     '1M_88Ceff_2MASS19211746+4402089_999K_3,0fer',
                                     '1M_88Ceff_2MASS19213437+4321500_999K_3,0fer',
                                     '1M_88Ceff_2MASS14474655+0103538_999K_3,0fer'], dtype='<U48')

fit_exclude_4GPa_perplex = np.array(['1M_88Ceff_2MASS04215269+5749018_999K_3,0fer',
                                     '1M_88Ceff_HIP58374_999K_3,0fer',
                                     '1M_88Ceff_2MASS19064452+4705535_999K_3,0fer',
                                     '1M_88Ceff_2MASS19125618+4031152_999K_3,0fer',
                                     '1M_88Ceff_2MASS19490218+4650354_999K_3,0fer',
                                     '1M_88Ceff_2MASS19211746+4402089_999K_3,0fer',
                                     '1M_88Ceff_2MASS19213437+4321500_999K_3,0fer',
                                     '1M_88Ceff_2MASS14474655+0103538_999K_3,0fer'], dtype='<U48')

fig1 = xplot_p_grid('melts', components, pressures)
fig2 = xplot_p_grid('perplex', components, pressures)

# plt.show()
