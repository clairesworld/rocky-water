import matplotlib.pyplot as plt
import matplotlib
from matplotlib.gridspec import GridSpec
import numpy as np
from useful_and_bespoke import colorize, cornertext, dark_background, iterable_not_string, colourbar, not_string
import perplexdata as px
from mpl_toolkits.axisartist.axislines import SubplotZero
from saturation import TO
import parameters as p
import main as rw

fig_path = '/home/claire/Works/rocky-water/figs_scratch/'  # path to save figures
test_phases_order = ['cpx', 'gt', 'opx', 'hpcpx', 'ol', 'pl', 'wad', 'ring', 'st', 'capv', 'wus', 'aki', 'cf',
                     'pv', 'ppv']  # phases to include in coloured composition plots

c_Earth = 'xkcd:olive green'


def annotate_ax_continuation(ax, which):
    # end axis with arrow to suggest it continues
    if which == 'y':
        ax.plot((0), (1), ls="", marker="^", ms=10, color="k",
                transform=ax.get_xaxis_transform(), clip_on=False)  # for vertical pressure
    elif which == 'x':
        ax.plot(1, (0), ls="", marker=">", ms=10, color="k",
                transform=ax.get_yaxis_transform(), clip_on=False)
    return ax


def single_structure(dat, fig_path=fig_path, save=True, **kwargs):
    fig, ax = plt.subplots(4, 2, sharex=True, figsize=[12, 8])

    ax[0, 0].plot(dat.radius / 1000, dat.density, color='black', alpha=1)
    ax[0, 0].axvline(dat.radius[dat.i_cmb] / 1000, ls='--')
    ax[0, 0].set_ylabel(r"$\rho$ (kg m$^{-3}$)")

    ax[1, 0].plot(dat.radius / 1000, dat.gravity, color='black', alpha=1)
    ax[1, 0].set_ylabel("$g$ (m s$^{-2}$)")

    ax[0, 1].plot(dat.radius / 1000, dat.pressure / 10 ** 9, color='black', alpha=1)
    ax[0, 1].set_ylabel("$P$ (GPa)")

    ax[1, 1].plot(dat.radius / 1000, dat.temperature, color='black', alpha=1)
    ax[1, 1].set_ylabel("$T$ (K)")

    ax[2, 0].plot(dat.radius / 1000, dat.alpha * 10 ** 5, color='black', alpha=1)
    ax[2, 0].set_ylabel(r"$\alpha$ ($10^{-5}$ K$^{-1}$)")

    ax[2, 1].plot(dat.radius / 1000, dat.cp, color='black', alpha=1)
    ax[2, 1].set_ylabel("$C_p$ (J kg$^{-1}$ K$^{-1}$)")

    # ax[3, 0].plot(dat.radius[1:] / 1000, dat.mass[1:], color='black', alpha=1)
    # ax[3, 0].set_ylabel("$\delta m$ (kg)")
    # ax[3, 0].set_xlabel("radius (km)")

    ax[3, 1].plot(dat.radius / 1000, dat.cum_mass, color='black', alpha=1)
    ax[3, 1].set_ylabel("$m$ (kg)")
    ax[3, 1].set_xlabel("radius (km)")

    plt.suptitle(dat.name)
    plt.tight_layout()
    if save:
        fig.savefig(fig_path + dat.name + '_structure.png', bbox_inches='tight')
    else:
        plt.show()


def single_composition(dat, which='pressure', modality_type='phase', comp_stacked=True, fig_path=fig_path, save=True,
                       show=False,
                       xlabel=None, ylabel=None, cmap='tab20', labelsize=16, plot_phases_order=None, p_max=None,
                       verbose=False, fig=None, ax=None, title=None, extension='.png', make_legend=True, **kwargs):
    """ comp_stacked: plot cumulative modes
    p_max: maximum pressure to plot in GPa (e.g., to see UM and MTZ better"""
    if which != 'pressure':
        raise NotImplementedError('independent variables other than pressure not implemented')

    # setup plot
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    x, xlabelauto = set_independent_var(dat, which)
    if modality_type == 'phase':
        ylabelauto = 'Phase abundance (wt %)'
        filename = dat.name + '_composition'
    elif modality_type == 'water':
        ylabelauto = 'Water modality in NAMs (wt %)'
        filename = dat.name + '_watermodality'
    if xlabel is None:
        xlabel = xlabelauto
    if ylabel is None:
        ylabel = ylabelauto
    if title is None:
        title = dat.name
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    fig.suptitle(title, fontsize=labelsize)

    if p_max is None:
        p_max = np.max(x)
    else:
        ax = annotate_ax_continuation(ax, 'x')
    ax.set_xlim(0, p_max)

    # get modes of each phase from dataframe
    if plot_phases_order is None:  # use all phases from dataframe (default)
        phases = rw.rename_phases(dat.phases_px)
    else:  # use a fixed order of phases - neater but may miss some
        phases = plot_phases_order
        if verbose:
            print('using input phase order', phases, '\noriginal phases are:', dat.phases_px)

    # get colours from cmap
    cmap_vals = matplotlib.cm.get_cmap(cmap)
    n_phases = len(phases)
    colours = [cmap_vals(j) for j in np.arange(0, n_phases) / n_phases]

    for ii, phase in enumerate(phases):
        # loop over phases and either plot instantly (not stacked) or create stackplot matrix
        if modality_type == 'phase':
            col = 'X_' + phase.lower()  # weight fraction of phase
        elif modality_type == 'water':
            col = 'frac_h2o_' + phase.lower()  # fraction of layer's water in this phase
        try:
            y = dat.df_all[col].to_numpy(dtype=float) * 100
        except KeyError:
            print(col, 'not found in df_all')
            if plot_phases_order is not None:  # this should always be true...
                y = np.zeros(len(dat.df_all), dtype=float)  # assign null value so phase order is retained
            else:
                raise Exception('phase name mismatch:', phase)
        if comp_stacked:
            # add to array for stackplot
            if ii == 0:
                y_stacked = y  # initialise
            else:
                y_stacked = np.vstack((y_stacked, y))
        else:
            # or, if not a stackplot, plot immediately
            ax.plot(x, y, c=colours[ii], label=phase)

            # make stackplot at this point
    if comp_stacked:
        ax.stackplot(x, y_stacked, labels=phases, colors=colours)
        ax.set_ylim(0, 100)
        # ax.scatter(x, [50]*len(x), marker='.', s=20, c='k', zorder=100)  # show perple_x resolution

    if make_legend:
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
    if p_max is None:
        p_max = np.max(x)
    else:
        # assume not all pressures are given, annotate axis to indicate there's more and not an error
        ax = annotate_ax_continuation(ax, 'x')
    ax.set_xlim(0, p_max)

    if save:
        plt.tight_layout()
        fig.savefig(fig_path + filename + extension, bbox_inches='tight')
    if show:
        plt.show()
    return fig, ax


def set_independent_var(dat, var_name):
    if var_name in ('p', 'pressure', 'P'):
        y = dat.df_all['P(bar)'].to_numpy(dtype=float)
        label = 'Pressure (GPa)'
        scale = 1e-4  # bar to GPa
    elif var_name in ('T', 'temperature'):
        y = dat.df_all['T(K)'].to_numpy(dtype=float)
        label = 'T (K)'
        scale = 1
    elif var_name in ('r', 'radius'):
        y = dat.df_all['r'].to_numpy(dtype=float)
        label = 'Radius (km)'
        scale = 1e-3
    elif var_name in ('z', 'depth'):
        y = dat.R_p - dat.df_all['r'].to_numpy(dtype=float)
        label = 'Depth (km)'
        scale = 1e-3
    return y * scale, label


def profile(dat, parameter, independent_ax='p', reverse_y=True, ax_label=None, scale=1, c='k', lw=2, alpha=1,
            xmin=None, xmax=None, ymin=None, ymax=None, ax_ticks=None, label_x=True, label_y=True,
            fig=None, ax=None, log=False, leg_label=None, labelsize=14, orientation='vertical',
            y2var=None, y2label=None, y2scale=None, save=True, fig_path=fig_path, **kwargs):
    if fig is None and ax is None:
        fig, ax = plt.subplots(1, 1)

    # get data and axis labels
    x = dat.df_all[parameter].to_numpy() * scale
    y, independent_ax_label = set_independent_var(dat, independent_ax)
    if ax_label is None:
        ax_label = parameter

    if orientation == 'vertical':
        ax.plot(x, y, c=c, lw=lw, label=leg_label, alpha=alpha)
        if log:
            ax.set_xscale('log')
        if ax_ticks is not None:
            ax.set_xticks(ax_ticks)
        if reverse_y:
            ax.invert_yaxis()
        xlabel = ax_label
        ylabel = independent_ax_label
    elif orientation == 'horizontal':
        ax.plot(y, x, c=c, lw=lw, label=leg_label, alpha=alpha)
        if log:
            ax.set_yscale('log')
        if ax_ticks is not None:
            ax.set_yticks(ax_ticks)
        xlabel = independent_ax_label
        ylabel = ax_label

    # set plot limits and labels
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin, xmax)
    if label_x:
        ax.set_xlabel(xlabel, fontsize=labelsize)
    if label_y:
        ax.set_ylabel(ylabel, fontsize=labelsize)

    if y2var is not None:
        raise NotImplementedError('second independent axis not implemented')
        # y2, y2label = set_independent_var(dat, y2var, y2scale)
        # ax2 = ax.twinx()
        # ax2.plot(x, y2, c=c, lw=lw, label=label, **kwargs)
        # ax2.set_ylabel(y2label, fontsize=labelsize)
        # ax2.set_ylim(min(y2), max(y2))

    if save:
        plt.tight_layout()
        fig.savefig(fig_path + dat.name + '_' + parameter + '.png', bbox_inches='tight')
    return fig, ax


def melting_curve(dat, labelsize=14, fig_path=fig_path):
    p = dat.pressure[dat.i_cmb:0] * 1e-9  # GPa
    T_a = dat.temperature[dat.i_cmb:0]  # K
    # Fe_num = dat.wt_oxides['FeO'] / dat.wt_oxides['MgO']
    # T_sol = 1409.15 + 134.2 * p - 6.581 * p ** 2 + 0.1054 * p ** 3 + (102.0 + 64.1 * p - 3.62 * p ** 2) * (0.1 - Fe_num)
    # T_m_stx = 5400 * (p / 140) ** 0.48 / (1 - np.log(1 - Fe_num))
    T_melt = dat.melting_temperature()

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    ax.plot(T_a, p, label='Mantle adiabat', c='w', ls='-')
    ax.plot(T_melt, p, label='Peridotite solidus (Fiquet+ 2010)', c='w', ls='--')

    ax.invert_yaxis()
    ax.set_xlabel('T (K)', fontsize=labelsize)
    ax.set_ylabel('p (GPa)', fontsize=labelsize)

    # # compare to Dorn parameterisation
    # i_189 = np.argmax(p >= 189.75)
    # x_FeO = dat.wt_oxides['FeO'] * 1e-2  # wt fraction FeO in mantle, lowers melting tempGt
    # a1 = 1831
    # a2 = 4.6
    # a3 = 0.33    # compare to Dorn parameterisation
    # i_189 = np.argmax(p >= 189.75)
    # x_FeO = dat.wt_oxides['FeO'] * 1e-2  # wt fraction FeO in mantle, lowers melting tempGt
    # a1 = 1831
    # a2 = 4.6
    # a3 = 0.33
    # b1 = 5400
    # b2 = 140
    # b3 = 0.48
    # c1 = 360
    # c2 = 0.0818
    # c3 = 102
    # c4 = 64.1
    # c5 = 3.62
    # if i_189 == 0:  # no pressures above this
    #     T_melt_Dorn = a1 * (1 + p / a2) ** a3 + c1 * (c2 - x_FeO)
    # else:
    #     T_melt_Dorn[:i_189] = a1 * (1 + p / a2) ** a3 + c1 * (c2 - x_FeO)
    #     T_melt_Dorn[i_189:] = b1 * (p / b2) ** b3 + (c3 + c4 * p - c5 * p ** 2) * (c2 - x_FeO)
    # b1 = 5400
    # b2 = 140
    # b3 = 0.48
    # c1 = 360
    # c2 = 0.0818
    # c3 = 102
    # c4 = 64.1
    # c5 = 3.62
    # if i_189 == 0:  # no pressures above this
    #     T_melt_Dorn = a1 * (1 + p / a2) ** a3 + c1 * (c2 - x_FeO)
    # else:
    #     T_melt_Dorn[:i_189] = a1 * (1 + p / a2) ** a3 + c1 * (c2 - x_FeO)
    #     T_melt_Dorn[i_189:] = b1 * (p / b2) ** b3 + (c3 + c4 * p - c5 * p ** 2) * (c2 - x_FeO)
    # ax.plot(T_melt_Dorn, p, label='solidus Dorn eq (dry)')

    leg = ax.legend(frameon=False, fontsize=12, bbox_to_anchor=(-0.1, 1.02, 1., .102), loc='lower left', )
    # fig, *ax = dark_background(fig, ax, )
    fig.savefig(fig_path + 'melt.png', bbox_inches='tight', bbox_extra_artists=(leg,),
                facecolor=fig.get_facecolor())
    plt.show()


def pop_hist1D(dats, x_name, scale=1, earth=None, xlabel=None, title=None, c_hist='k', ls='-', fig_path=fig_path,
               filename=None, extension='.png', save=True, show=True, data_label=None, fig=None, ax=None,
               xlim=None, labelsize=12, legsize=12, **kwargs):
    """ histogram of 1 or more variables on same x axis. x_name, data_label, and ls can be list, to differentiate """
    if xlabel is None:
        xlabel = x_name
    if filename is None:
        filename = x_name + '_hist'

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    def do_plot(attr, c, label, linestyle):
        x = []
        for pl in dats:
            try:
                x.append(eval('pl.' + attr))
            except AttributeError:
                pass  # blank data, didn't work because star not measured maybe
        ax.hist([a * scale for a in x], color=c, label=label, ls=linestyle, **kwargs)
        if earth:
            ax.axvline(eval('earth.' + attr) * scale, label='Earth', c=c_Earth, ls=linestyle)

    if iterable_not_string(x_name):
        for ii, attr in enumerate(x_name):
            do_plot(attr, c_hist, data_label[ii], ls[ii])
    else:
        do_plot(x_name, c_hist, data_label, ls)

    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.legend(frameon=False, fontsize=legsize)
    cornertext(ax, text='n = ' + str(len(dats)), pos='top left', size=legsize, pad=0.02)
    plt.title(title)
    if xlim is not None:
        ax.set_xlim(xlim)

    if save:
        plt.savefig(fig_path + filename + extension, bbox_inches='tight')
    if show:
        plt.show()
    return fig, ax


def pop_scatter(dats, x_name, y_name, x_scale=1, y_scale=1, title=None, xlabel=None, ylabel=None, filename=None,
                data_label=None, earth=False, fig_path=fig_path, save=True, show=True, Earth_kwargs={},
                extension='.png', **kwargs):
    # for list of planet objects, make scatter plot of 2 properties
    x = []
    y = []
    if xlabel is None:
        xlabel = x_name
    if ylabel is None:
        ylabel = y_name
    if filename is None:
        filename = x_name + '_' + y_name + '_scatter'
    for pl in dats:
        try:
            xi = eval('pl.' + x_name)
            yi = eval('pl.' + y_name)
            x.append(xi)
            y.append(yi)
        except KeyError:
            print(pl.name, 'does not have attribute', x_name, 'or', y_name)

    fig, ax = plt.subplots(1, 1)
    ax.scatter([a * x_scale for a in x], [b * y_scale for b in y], label=data_label, **kwargs)
    if earth:
        ax.scatter(eval('earth.' + x_name) * x_scale, eval('earth.' + y_name) * y_scale, label='Earth', c=c_Earth,
                   marker='$\oplus$', s=100)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False)
    cornertext(ax, text='n = ' + str(len(x)), pos='top left', size=10, pad=0.02)
    plt.title(title)

    if save:
        plt.savefig(fig_path + filename + extension, bbox_inches='tight')
    if show:
        plt.show()
    return fig, ax


def dict_across_pops(dirs, x_name, y_name, v_name=None, exclude_params=None, exclude_min=-1e50, exclude_max=1e50,
                     **kwargs):
    """ build dict for same planets across multiple datasets (e.g., changing mass) for use with next two fns """
    # get star names from first dir
    datsdict = {}  # will save each star as key
    # initialise dictionary
    dats0 = rw.read_dir(dirs[0])
    for dat in dats0:
        if (exclude_params is not None and exclude_points(dat, exclude_params, exclude_min, exclude_max)) or (
                exclude_params is None):
            if v_name is None:  # no colorbar
                datsdict[dat.star] = ([], [])  # placeholder for x, y
            else:
                datsdict[dat.star] = ([], [], [])  # placeholder for x, y, v
            print('added key', dat.star)

    # load x, y values
    for dir_ in dirs:
        datsj = rw.read_dir(dir_)  # [:head]
        for dat in datsj:
            # print('attempting', dat.star, 'in', dir_)
            try:  # make sure data is there
                xj = eval('dat.' + x_name)
                yj = eval('dat.' + y_name)
                if v_name is not None:
                    vj = eval('dat.' + v_name)
                datsdict[dat.star][0].append(xj)  # datsdict is dict of 2-tuples of lists
                datsdict[dat.star][1].append(yj)
                if v_name is not None:
                    datsdict[dat.star][2].append(vj)
            except AttributeError:
                # if data is missing, remove this star everywhere
                print(dat.star, 'does not have', x_name, 'or', y_name)
                datsdict.pop(dat.star, None)
            except KeyError as e:
                print(e)
                pass  # move to next star folder
    return datsdict


def compare_pop_fillbetween(dirs, x_name, y_name, x_scale=1, y_scale=1, xlog=False, ylog=False, ylabel=None,
                            xlabel=None,
                            title=None, filename=None, fig_path=fig_path,
                            save=True, show=True, extension='.png', labelsize=12, show_med=True, c='xkcd:peach',
                            legsize=10, ticksize=12, earth=None, head=-1, xlim=None, ylim=None, patch_kwargs=None,
                            line_kwargs=None,
                            dark=False, fig=None, ax=None, sigma=1, **kwargs):
    """ for a list of directories containing runs with different x_name values, plot y_name vs. x_name as fillbetween
        dats to plot is taken from first entry in directory list (dirs) so make sure that dir is complete """
    if xlabel is None:
        xlabel = x_name
    if ylabel is None:
        ylabel = y_name
    if filename is None:
        filename = 'pops_dist_' + x_name + '_' + y_name

    datsdict = dict_across_pops(dirs, x_name, y_name, **kwargs)

    # plot
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    # get into plotting format: y should be nested list of different runs' y-value at each x
    x = datsdict[list(datsdict.keys())[0]][0]  # should be the same for all
    y = [[] for _ in range(len(x))]
    for jj in range(len(x)):
        for star in list(datsdict.keys()):
            try:
                y[jj].append(datsdict[star][1][jj])
            except IndexError as e:
                print(star, 'missing index', jj, ':', e)

    # get 1-sigma distribution of y
    if patch_kwargs is None:
        patch_kwargs = {'color': c, 'alpha': 0.5}
    if line_kwargs is None:
        line_kwargs = {'color': c, 'lw': 3}

    if sigma == 1:
        q = [0.16, 0.50, 0.84]  # percentiles of 1 standard deviation above and below mean
    elif sigma == 2:
        q = [0.0227, 0.50, 0.9773]  # percentiles of 1 standard deviation above and below mean
    else:
        raise NotImplementedError('sigma level ' + str(sigma) + ' not implemented')

    y_min, y_mid, y_max = [], [], []
    for ys in y:
        mini, midi, maxi = np.quantile(ys, q)
        y_min.append(mini)
        y_mid.append(midi)
        y_max.append(maxi)

    ax.fill_between(np.array(x)*x_scale, np.array(y_min)*y_scale, np.array(y_max)*y_scale,
                    label=str(sigma) + '-sigma distribution',**patch_kwargs)
    if show_med:
        ax.plot(np.array(x)*x_scale, np.array(y_mid)*y_scale, **line_kwargs)

    # define min and max
    ax.plot(np.array(x)*x_scale, np.array(y_min)*y_scale, lw=0.5, c=patch_kwargs['color'])
    ax.plot(np.array(x)*x_scale, np.array(y_max)*y_scale, lw=0.5, c=patch_kwargs['color'])

    # add earth point if desired
    if earth:
        ax.scatter(eval('earth.' + x_name) * x_scale, eval('earth.' + y_name) * y_scale, label='Earth', c=c_Earth,
                   marker='$\oplus$', s=100, zorder=100)

    # finish plot labelling, axis limits etc
    cornertext(ax, 'n = ' + str(len(list(datsdict.keys())[:head])), pos='bottom right', size=legsize, c='k')
    ax.legend(frameon=False, fontsize=legsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_title(title)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if save:
        fig.savefig(fig_path + filename + extension, bbox_inches='tight')
    if show:
        plt.show()
    return fig, ax


def compare_pop_scatter(dirs, x_name, y_name, v_name=None, x_scale=1, y_scale=1, xlog=False, ylog=False, ylabel=None,
                        xlabel=None,
                        title=None, filename=None, fig_path=fig_path,
                        save=True, show=True, extension='.png', labelsize=12, lw=0.5, vmin=None, vmax=None,
                        legsize=10, ticksize=12, c='k', earth=None, head=-1, xlim=None, ylim=None, alpha=0.5,
                        dark=False, fig=None, ax=None, **kwargs):
    """ for a list of directories containing runs with different x_name values, plot y_name vs. x_name scatterplot
    dats to plot is taken from first entry in directory list (dirs) so make sure that dir is complete """
    if xlabel is None:
        xlabel = x_name
    if ylabel is None:
        ylabel = y_name
    if filename is None:
        filename = 'pops_scatter_' + x_name + '_' + y_name

    datsdict = dict_across_pops(dirs, x_name, y_name, v_name, **kwargs)

    # add colormap values
    if v_name is not None:
        v = [datsdict[star][2][0] for star in
             list(datsdict.keys())]  # 3-tuple of datsdict[star][2] should have identical values
    else:
        v = range(len(datsdict))
    try:
        colors = colorize(v, c, vmin=vmin, vmax=vmax)[0]
    except:
        # not a valid colormap
        colors = [c] * len(datsdict)

    # plot
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    for kk, star in enumerate(list(datsdict.keys())[:head]):
        x = datsdict[star][0]
        y = datsdict[star][1]
        ax.plot([xx * x_scale for xx in x], [yy * y_scale for yy in y], c=colors[kk], alpha=alpha, lw=lw, marker='o',
                markersize=5)

    # add colourbar
    if v_name is not None:
        # print('v')
        cbar = colourbar(mappable=None, vector=v, ax=ax, label='Mg/Si', labelsize=labelsize, ticksize=ticksize,
                         ticks=None, ticklabels=None, labelpad=17, loc='right', vmin=vmin, vmax=vmax,
                         rot=None, cmap=c, pad=0.05, log=False)

    # add earth point if desired
    if earth:
        ax.scatter(eval('earth.' + x_name) * x_scale, eval('earth.' + y_name) * y_scale, label='Earth', c=c_Earth,
                   marker='$\oplus$', s=100, zorder=100)

    # finish plot labelling, axis limits etc
    cornertext(ax, 'n = ' + str(len(list(datsdict.keys())[:head])), pos='bottom right', size=legsize, c='k')
    ax.legend(frameon=False, fontsize=legsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_title(title)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if save:
        fig.savefig(fig_path + filename + extension, bbox_inches='tight')
    if show:
        plt.show()
    return fig, ax


def single_phase_subfig(dat, var_name, var_scale=1, var_log=False, var_label=None, vertical_pressure=False,
                        title=None, annotation=None, p_max=None,
                        filename=None, fig_path=fig_path, save=True, show=True,
                        extension='.png', phase_order=None,
                        labelsize=12, legsize=10, ticksize=12, xpad=10, ypad=10,
                        linec='xkcd:navy', lw=2, dark=False, **kwargs):
    if filename is None:
        filename = dat.name + '_modes_' + var_name
    if var_label is None:
        var_label = var_name

    # setup figure
    fig = plt.figure()
    fig.suptitle(title, fontsize=labelsize)
    if not vertical_pressure:  # i.e., horizontal pressure axis, profile on top subplot
        gs = GridSpec(2, 1, height_ratios=[1, 4])
        ax_var = fig.add_subplot(gs[0])  # axis for plotting profile
        ax_mod = fig.add_subplot(gs[1])  # axis for plotting phase modes
        if var_log:
            ax_var.set_yscale('log')
    else:  # i.e., vertical pressure axis, profile on right subplot
        gs = GridSpec(1, 2, width_ratios=[2, 1])
        ax_mod = fig.add_subplot(gs[0])
        ax_var = fig.add_subplot(gs[1])

    # make composition subplot
    fig, ax_mod = single_composition(dat, which='pressure', comp_stacked=True, save=False,
                                     xlabel=None, ylabel=None, labelsize=labelsize,
                                     # use default labels because always in labelling position
                                     plot_phases_order=phase_order, p_max=p_max, fig=fig, ax=ax_mod, make_legend=False)

    leg = ax_mod.legend(bbox_to_anchor=(1.01, 0), loc='lower left', frameon=False, fontsize=legsize)

    # plot profile
    if p_max is None:
        p_max = dat.pressure[dat.i_cmb + 1] * 1e-9  # need to say something otherwise bug in label limits
    fig, ax_var = profile(dat, var_name, independent_ax='p', scale=var_scale, c=linec, lw=lw, alpha=1,
                          ax_label=var_label,
                          xmin=0, xmax=p_max, fig=fig, ax=ax_var, log=var_log, leg_label=None, orientation='horizontal',
                          save=False, labelsize=labelsize, **kwargs)
    if not vertical_pressure:
        ax_var.set_xlabel('')  # hide label
        # ax_var.set_xticks([])

    if annotation is not None:
        ax_var = cornertext(ax_var, annotation, size=labelsize, c=linec)

    axes = (ax_var, ax_mod)
    for ax in axes:
        # ax.set_xlim(0, p_max)
        ax.tick_params(axis='both', which='major', labelsize=ticksize)
    # plt.tight_layout()
    if dark:
        fig, *axes = dark_background(fig, axes, )
    if save:
        fig.savefig(fig_path + filename + extension, bbox_extra_artists=(leg,), bbox_inches='tight',
                    facecolor=fig.get_facecolor())
    if show:
        plt.show()
    return fig, axes


def exclude_points(dat, exclude_params, minval=-1e50, maxval=1e50):
    """ function to plot only planets inside a range for a given parameter(s)
    parameter must be scalar but can be any string that works with eval('dat.' + string)"""
    if iterable_not_string(exclude_params):
        # given a list of names, minval and maxval must also be a list
        for ii, param in enumerate(exclude_params):
            p = eval('dat.' + param)
            if p < minval[ii]:
                return None
            if p > maxval[ii]:
                return None
    else:
        p = eval('dat.' + exclude_params)
        if p < minval:
            return None
        if p > maxval:
            return None
    return dat


def find_rich_planets(dir_name, phase, n=1, output_base='output/', **kwargs):
    """ find n planets most rich in phase """
    dats = rw.read_dir(px.perplex_path_default + output_base + dir_name + '/')

    phase_props = []
    names = []
    for dat in dats:
        try:
            phase_props.append(np.sum(eval('dat.df_all.X_' + phase).to_numpy()) / len(dat.df_all))  # wt frac of mantle
            names.append(dat.name)
        except AttributeError:  # phase not found
            pass
    yx = list(zip(phase_props, names))
    yx.sort(reverse=True)
    names_sorted = [x for y, x in yx]
    # props_sorted = [y for y, x in yx]
    print('\n', n, 'planets with the most', phase)
    print('-------------------------------')
    for val, pl in yx[:n]:
        print(pl, '|', val)
    return names_sorted[:n]
