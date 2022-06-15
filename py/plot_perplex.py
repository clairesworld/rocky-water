import matplotlib.pyplot as plt
import matplotlib
from matplotlib.gridspec import GridSpec
import numpy as np
from useful_and_bespoke import colorize, cornertext, dark_background, iterable_not_string, colourbar, not_string
import perplexdata as px
from mpl_toolkits.axisartist.axislines import SubplotZero
from saturation import TO
from pprint import pprint
import parameters as p
import main as rw

fig_path = '/home/claire/Works/rocky-water/figs_scratch/'  # path to save figures
test_phases_order = ['cpx', 'gt', 'opx', 'hpcpx', 'ol', 'pl', 'wad', 'ring', 'st', 'capv', 'wus', 'aki', 'cf',
                     'pv', 'ppv']  # phases to include in coloured composition plots

c_Earth = 'xkcd:olive green'


def get_phase_display_names(phases_list):
    alias_dict = {}
    for s in phases_list:
        alias_dict[s] = s
        if s == 'an':
            alias_dict[s] = 'plg'
        elif s == 'st':
            alias_dict[s] = 'stv'
        elif s == 'wus':
            alias_dict[s] = 'fp'
        elif s == 'hpcpx':
            alias_dict[s] = 'hcpx'
    return alias_dict


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
                       show=False, legtitle=None, override_ax_arrow=False, ylabelpad=None, xlabelpad=None,
                       xlabel=None, ylabel=None, cmap='tab20', labelsize=16, plot_phases_order=None, p_max=None,
                       verbose=False, fig=None, ax=None, title=None, extension='.png', make_legend=True,
                       leg_bbox_to_anchor=(1, 1), legsize=10, **kwargs):
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
    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
    fig.suptitle(title, fontsize=labelsize)

    # p_max_orig = p_max
    if p_max is None:
        p_max = np.max(x)
    elif not override_ax_arrow:
        ax = annotate_ax_continuation(ax, 'x')
    ax.set_xlim(1000e-4, p_max)

    # get modes of each phase from dataframe
    if plot_phases_order is None:  # use all phases from dataframe (default)
        phases = rw.rename_phases(dat.phases_px)  # list of strings
    else:  # use a fixed order of phases - neater but may miss some
        phases = plot_phases_order
        if verbose:
            print('using input phase order', phases, '\noriginal phases are:', dat.phases_px)
    phase_map = get_phase_display_names(phases)

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
            ax.plot(x, y, c=colours[ii], label=phase_map[phase])

    # make stackplot at this point
    if comp_stacked:
        ax.stackplot(x, y_stacked, labels=[phase_map[s] for s in phases], colors=colours)
        ax.set_ylim(0, 100)
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


def profile(dat, parameter, independent_ax='pressure', reverse_y=True, ax_label=None, scale=1, c='k', lw=2, ls='-',
            alpha=1,
            xmin=None, xmax=None, ymin=None, ymax=None, ax_ticks=None, label_x=True, label_y=True, labelsize=14,
            fig=None, ax=None, log=False, leg_label=None, orientation='vertical', figsize=None, fname=None,
            legtitle=None,
            y2var=None, y2label=None, y2scale=None, save=True, fig_path=fig_path, legsize=12, leg_bbox_to_anchor=None,
            **kwargs):
    if fig is None and ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)

    # get data and axis labels
    try:
        x = dat.df_all[parameter].to_numpy() * scale
    except KeyError as e:
        print(dat.df_all.head())
        raise e
    if fname is None:
        fname = dat.name + '_' + parameter

    y, independent_ax_label = set_independent_var(dat, independent_ax)

    # set plot limits and labels
    if ymax == 'max':
        ymax = np.max(y)
    ax.set_ylim(ymin, ymax)
    if (xmin is not None) or (xmax is not None):
        ax.set_xlim(xmin, xmax)

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
        ax.plot(y, x, c=c, lw=lw, ls=ls, label=leg_label, alpha=alpha)
        if log:
            ax.set_yscale('log')
        if ax_ticks is not None:
            ax.set_yticks(ax_ticks)
        xlabel = independent_ax_label
        ylabel = ax_label

    if label_x:
        ax.set_xlabel(xlabel, fontsize=labelsize)
    if label_y:
        ax.set_ylabel(ylabel, fontsize=labelsize)

    leg = ax.legend(frameon=False, fontsize=legsize, title=legtitle, loc='upper left',
                    bbox_to_anchor=leg_bbox_to_anchor,
                    )
    if legtitle is not None:
        leg.get_title().set_fontsize(legsize)  # legend 'Title' fontsize

    if y2var is not None:
        raise NotImplementedError('second independent axis not implemented')
        # y2, y2label = set_independent_var(dat, y2var, y2scale)
        # ax2 = ax.twinx()
        # ax2.plot(x, y2, c=c, lw=lw, label=label, **kwargs)
        # ax2.set_ylabel(y2label, fontsize=labelsize)
        # ax2.set_ylim(min(y2), max(y2))

    if save:
        plt.tight_layout()
        fig.savefig(fig_path + fname + '.png', bbox_inches='tight')
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


def pop_hist1D(dats, x_name, scale=1, earth=None, xlabel=None, title=None, c_hist='k', ls_stats='-', fig_path=fig_path,
               filename=None, extension='.png', save=True, show=True, data_label=None, fig=None, ax=None,
               annotate_n=True, xmin=None, xmax=None,
               xlim=None, labelsize=12, legsize=12, bins=None, showmedian=False, showsigma=False, **kwargs):
    """ histogram of 1 or more variables on same x axis. x_name, data_label, and ls can be list, to differentiate """
    if xlabel is None:
        xlabel = x_name
    if filename is None:
        filename = x_name + '_hist'

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    def do_plot(attr, c, label, ls_stats):
        x = []
        for dat in dats:
            try:
                x.append(eval('dat.' + attr) * scale)
            except (AttributeError, KeyError):
                pass  # blank data, didn't work because star not measured maybe
            except TypeError:
                # attr is None?
                pass
        arr = np.array(x)
        # print('hist; y < 0.2', yarr[yarr < 0.2], 'n', len(yarr[yarr < 0.2]))
        if xmin is not None:
            arr = arr[arr > xmin]
        if xmax is not None:
            arr = arr[arr < xmax]
        ax.hist(arr, color=c, label=label, bins=bins, **kwargs)
        if earth:
            ax.axvline(eval('earth.' + attr) * scale, label='Earth', c=c_Earth)
        if showmedian:
            med = np.median(arr)
            print(attr, 'median', med)
            ax.axvline(med, label='median', c=c, ls=ls_stats)
        if showsigma:
            if showsigma == 1:
                q = [16, 84]
            elif showsigma == 2:
                q = [2.27, 97.93]
            qs = np.percentile(arr, q, method='nearest')
            print('qs', qs)
            siglabel = '$\pm$' + str(showsigma) + '$\sigma$'
            for q in qs:
                ax.axvline(q, c='r', ls=ls_stats, label=siglabel)
                siglabel = None

    if iterable_not_string(x_name):
        for ii, attr in enumerate(x_name):
            do_plot(attr, c_hist, data_label[ii], ls_stats[ii])
    else:
        do_plot(x_name, c_hist, data_label, ls_stats)

    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.legend(frameon=False, fontsize=legsize)
    if annotate_n:
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
                ticksize=10, legsize=12, xlabelpad=None, ylabelpad=None, range_min=None, range_max=None,
                data_label=None, earth=False, fig_path=fig_path,
                extension='.png', fig=None, ax=None, xlim=None, ylim=None, labelsize=14, save=True, show=True, x=None,
                y=None,
                return_data=False, annotate_n=True, ms=200, **scatter_kwargs):
    # for list of planet objects, make scatter plot of 2 properties
    print('range_min', range_min)
    if x is None and y is None:
        x = []
        y = []
        if xlabel is None:
            xlabel = x_name
        if ylabel is None:
            ylabel = y_name
        if filename is None:
            filename = x_name + '_' + y_name + '_scatter'
        for dat in dats:
            try:
                xi = eval('dat.' + x_name)
                yi = eval('dat.' + y_name)
                if (range_min is None) or (range_min is not None and xi >= range_min):
                    if (range_max is None) or (range_max is not None and xi <= range_max):
                        x.append(xi * x_scale)
                        y.append(yi * y_scale)

            except KeyError:
                print(dat.name, 'does not have attribute', x_name, 'or', y_name)

    # create figure
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    print('max', y_name, '=', np.max(y))
    ax.scatter(x, y, label=data_label, **scatter_kwargs)
    if earth:
        ax.scatter(eval('earth.' + x_name) * x_scale, eval('earth.' + y_name) * y_scale, label='Earth', c='k',
                   marker='$\oplus$', s=ms, zorder=200)

    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.legend(frameon=False, fontsize=legsize)
    if annotate_n:
        cornertext(ax, text='n = ' + str(len(x)), pos='top left', size=8, pad=0.02)
    plt.title(title)
    if ylim is not None:
        ax.set_ylim(ylim)
    if xlim is not None:
        ax.set_xlim(xlim)

    if save:
        plt.savefig(fig_path + filename + extension, bbox_inches='tight')
    if show:
        plt.show()
    if return_data:
        return fig, ax, x, y
    return fig, ax


def pop_scatterhist_subplotwrapper(dats, x_name, y_names, y_scales=None, ylabels=None, filename=None,
                                   fig_path=fig_path, extension='.png', ylims=None,
                                   save=True, show=True, ratios=[7, 1], lims_histy=None, fig=None, **kwargs):
    # assume both histx and histy
    ncols = len(y_names)
    if y_scales is None:
        y_scales = [1] * ncols
    if ylabels is None:
        ylabels = [None] * ncols
    if ylims is None:
        ylims = [None] * ncols
    if lims_histy is None:
        lims_histu = [None] * ncols

    if fig is None:
        fig = plt.figure()
    gs = fig.add_gridspec(1 + ncols, 2, width_ratios=[ratios[0] + 2, ratios[1]],
                          height_ratios=[ratios[1]] + [ratios[0]] * ncols,
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    ax_histx = fig.add_subplot(gs[0, 0])
    axes = []
    axesy = []
    show_histx = True
    for ii in range(ncols):
        axes.append(fig.add_subplot(gs[ii + 1, 0], sharex=ax_histx))
        axesy.append(fig.add_subplot(gs[ii + 1, 1], sharey=axes[ii]))

        if y_names[ii] is not None:
            fig, axes[ii], ax_histx, axesy[ii] = pop_scatterhist(dats, x_name, y_names[ii], x_scale=1, y_scale=y_scales[ii],
                                                                 ylabel=ylabels[ii], show_histx=show_histx, show_histy=True,
                                                                 fig=fig, ax=axes[ii], ax_histx=ax_histx,
                                                                 ax_histy=axesy[ii], ylim=ylims[ii],
                                                                 save=False, show=False, lim_histy=lims_histy[ii], **kwargs)
            show_histx = False  # only plot once
    if save:
        plt.savefig(fig_path + filename + extension, bbox_inches='tight')
    if show:
        plt.show()
    return fig, axes, ax_histx, axesy


def pop_scatterhist(dats, x_name, y_name, x_scale=1, y_scale=1, title=None, xlabel=None, ylabel=None, filename=None,
                    show_histx=False, show_histy=False, histx_kwargs={}, histy_kwargs={}, data_label=None, earth=False,
                    fig_path=fig_path, range_min=None, range_max=None,
                    extension='.png', fig=None, ax=None, ax_histx=None, ax_histy=None, xlim=None, ylim=None,
                    save=True, show=True, ratios=[7, 1], lim_histx=None, lim_histy=None, **kwargs):
    x, y = [], []
    for dat in dats:
        try:
            xi = eval('dat.' + x_name)
            yi = eval('dat.' + y_name)
            if (range_min is None) or (range_min is not None and xi >= range_min):
                if (range_max is None) or (range_max is not None and xi <= range_max):
                    x.append(xi * x_scale)
                    y.append(yi * y_scale)
        except KeyError:
            print(dat.name, 'does not have attribute', x_name, 'or', y_name)

    if fig is None:
        fig = plt.figure()
    if (ax is None) and (not show_histx) and (not show_histy):
        ax = fig.add_subplot(111)
    elif show_histx and show_histy:
        gs = fig.add_gridspec(2, 2, width_ratios=[ratios[0] + 2, ratios[1]],
                              height_ratios=[ratios[1], ratios[0]],
                              left=0.1, right=0.9, bottom=0.1, top=0.9,
                              wspace=0.05, hspace=0.05)
        if ax is None:
            ax = fig.add_subplot(gs[1, 0])
        if ax_histx is None:
            ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
        if ax_histy is None:
            ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    elif show_histx:
        gs = fig.add_gridspec(2, 1, height_ratios=list(reversed(ratios)),
                              bottom=0.1, top=0.9,
                              hspace=0.05)
        if ax is None:
            ax = fig.add_subplot(gs[1, 0])
        if ax_histx is None:
            ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    elif show_histy:
        gs = fig.add_gridspec(1, 2, width_ratios=ratios,
                              left=0.1, right=0.9,
                              wspace=0.05)
        if ax is None:
            ax = fig.add_subplot(gs[1, 0])
        if ax_histy is None:
            ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    fig, ax = pop_scatter(dats, x_name, y_name, x_scale, y_scale, title, xlabel, ylabel, filename='',
                          data_label=data_label, earth=earth, fig_path=fig_path, range_min=range_min, range_max=range_max,
                          extension=extension, fig=fig, ax=ax, xlim=xlim, ylim=ylim, save=False,
                          show=False, x=x, y=y, **kwargs)

    if show_histx:
        if 'range' not in histx_kwargs.keys():
            histxrange = np.array(ax.get_xlim())
            if range_min is not None:
                histxrange[0] = range_min
            if range_max is not None:
                histxrange[1] = range_max
            histx_kwargs['range'] = histxrange
        ax_histx.hist(x, density=True, **histx_kwargs)
        ax_histx.spines.right.set_visible(False)
        ax_histx.spines.left.set_visible(False)
        ax_histx.spines.top.set_visible(False)
        ax_histx.tick_params(axis="both", labelbottom=False, labelleft=False)
        ax_histx.set_yticks([])
        ax_histx.set_ylim(lim_histx)

    if show_histy:
        ax_histy.hist(y, range=ylim, density=True, orientation='horizontal', **histy_kwargs)
        ax_histy.spines.bottom.set_visible(False)
        ax_histy.spines.right.set_visible(False)
        ax_histy.spines.top.set_visible(False)
        ax_histy.tick_params(axis="both", labelbottom=False, labelleft=False)
        ax_histy.set_xticks([])
        ax_histy.set_xlim(lim_histy)

    if save:
        plt.savefig(fig_path + filename + extension, bbox_inches='tight')
    if show:
        plt.show()
    return fig, ax, ax_histx, ax_histy


def dict_across_pops(dirs, x_name, y_name, v_name=None, exclude_params=None, exclude_min=-1e50, exclude_max=1e50,
                     subsample=None, verbose=False, **kwargs):
    """ build dict for same planets across multiple datasets (e.g., changing mass) for use with next two fns """
    # get star names from first dir
    datsdict = {}  # will save each star as key, datsdict is dict of 2-tuples of lists
    dats0 = rw.read_dir(dirs[0], subsample=subsample)
    for dat in dats0:
        if (exclude_params is not None and exclude_points(dat, exclude_params, exclude_min, exclude_max)) or (
                exclude_params is None):
            if v_name is None:  # no colorbar
                datsdict[dat.star] = ([], [])  # placeholder for x, y
            else:
                datsdict[dat.star] = ([], [], [])  # placeholder for x, y, v
            if verbose:
                print('added key', dat.star)

    # load x, y values
    for dir_ in dirs:
        datsj = rw.read_dir(dir_)  # [:head]
        for dat in datsj:
            # print('attempting', dat.star, 'in', dir_)
            if dat.star in datsdict:
                try:  # make sure data is there
                    xj = eval('dat.' + x_name)
                    yj = eval('dat.' + y_name)
                    if v_name is not None:
                        vj = eval('dat.' + v_name)
                    datsdict[dat.star][0].append(xj)  # datsdict[key][0] is list of x values
                    datsdict[dat.star][1].append(yj)  # datsdict[key][1] is list of x values
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
                            title=None, filename=None, fig_path=fig_path, dpi=300,
                            save=True, show=True, extension='.png', labelsize=12, show_med=True, c='xkcd:peach',
                            legsize=10, ticksize=12, earth=None, sun=None, head=-1, xlim=None, ylim=None, patch_kwargs=None,
                            line_kwargs=None, show_scaling_fn=None, figsize=(4, 4), scalinglabel=None, earth_real=None,
                            dark=False, fig=None, ax=None, sigma=1, show_n=True, datalabel=None, legend=True, **kwargs):
    """ for a list of directories containing runs with different x_name values, plot y_name vs. x_name as fillbetween
        dats to plot is taken from first entry in directory list (dirs) so make sure that dir is complete """
    if xlabel is None:
        xlabel = x_name
    if ylabel is None:
        ylabel = y_name
    if filename is None:
        filename = 'pops_dist_' + x_name + '_' + y_name
    if datalabel is None:
        datalabel = str(sigma) + '$\sigma$ distribution'

    print('x_name', x_name)
    datsdict = dict_across_pops(dirs, x_name, y_name, **kwargs)

    # plot
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)

    # add earth point if desired

    if sun:
        ax.scatter(eval('sun.' + x_name) * x_scale, eval('sun.' + y_name) * y_scale, label='Solar', c='k', #alpha=0.4,
                   marker='*', s=100, zorder=100)
    if earth:
        ax.scatter(eval('earth.' + x_name) * x_scale, eval('earth.' + y_name) * y_scale, label='Earth', c='xkcd:dark blue',
                   marker='$\oplus$', s=200, zorder=100)
        print('earth:', y_name, eval('earth.' + y_name) * y_scale)


    if earth_real:
        ax.scatter(eval('earth.' + x_name) * x_scale, eval('earth.' + y_name) * y_scale, label='Earth (observed)',
                   c='xkcd:silver',
                   marker='$\oplus$', s=200, zorder=100)

    # get into plotting format: y should be nested list of different runs' y-value at each x
    x = datsdict[list(datsdict.keys())[0]][0]  # should be the same for all
    print('x', x)
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
        line_kwargs = {'color': c, 'lw': 2}

    if sigma == 1:
        q = [0.16, 0.50, 0.84]  # percentiles of 1 standard deviation above and below mean
    elif sigma == 2:
        q = [0.0227, 0.50, 0.9773]  # percentiles of 2 standard deviation above and below mean
    else:
        raise NotImplementedError('sigma level ' + str(sigma) + ' not implemented')

    y_min, y_mid, y_max = [], [], []
    for ys in y:
        mini, midi, maxi = np.quantile(ys, q)
        y_min.append(mini)
        y_mid.append(midi)
        y_max.append(maxi)

    ax.fill_between(np.array(x) * x_scale, np.array(y_min) * y_scale, np.array(y_max) * y_scale,
                    label=datalabel, **patch_kwargs)
    if show_med:
        ax.plot(np.array(x) * x_scale, np.array(y_mid) * y_scale, **line_kwargs)
        print('med', [ys*y_scale for ys in y_mid])

    # # define min and max
    # ax.plot(np.array(x) * x_scale, np.array(y_min) * y_scale, lw=0.5, c=patch_kwargs['color'])
    # ax.plot(np.array(x) * x_scale, np.array(y_max) * y_scale, lw=0.5, c=patch_kwargs['color'])
    # # print('x\n', np.array(x) * x_scale)
    # # print('\n\ny\n',  np.array(y_mid) * y_scale)

    if show_scaling_fn:
        # add naive scaling
        xhat = np.linspace(np.min(x), np.max(x))
        y_naive = show_scaling_fn(xhat)
        ax.plot(xhat * x_scale, y_naive * y_scale, c='k', lw=0.5, ls=':', label=scalinglabel)

    # finish plot labelling, axis limits etc
    if show_n:
        cornertext(ax, 'n = ' + str(len(list(datsdict.keys())[:head])), pos='top left', size=legsize, c='k')
    if legend:
        ax.legend(frameon=False, fontsize=legsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_title(title, fontsize=labelsize)
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
        fig.savefig(fig_path + filename + extension, bbox_inches='tight', dpi=dpi)
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
    dats to plot is taken from first entry in directory list (dirs) so make sure that dir is complete
    use lw=0 for no joining line"""
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


def composition_subfig(dat, var_name, var_scale=1, var_log=False, var_label=None, vertical_pressure=False,
                       title=None, annotation=None, p_max=None,
                       filename=None, fig_path=fig_path, save=True, show=True,
                       extension='.png', phase_order=None, show_UM=False,
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
                                     xlabel=None, ylabel=None, labelsize=labelsize, title=title,
                                     # use default labels because always in labelling position
                                     plot_phases_order=phase_order, p_max=p_max, fig=fig, ax=ax_mod, make_legend=False)

    leg = ax_mod.legend(bbox_to_anchor=(1.01, 0), loc='lower left', frameon=False, fontsize=legsize)

    # plot profile
    if p_max is None:
        p_max = dat.pressure[dat.i_cmb + 1] * 1e-9  # need to say something otherwise bug in label limits
    fig, ax_var = profile(dat, var_name, independent_ax='pressure', ax_label=var_label, scale=var_scale, c=linec, lw=lw,
                          alpha=1, xmin=0, xmax=p_max, fig=fig, ax=ax_var, log=var_log, leg_label=None,
                          labelsize=labelsize, orientation='horizontal', save=False, **kwargs)

    if not vertical_pressure:
        ax_var.set_xlabel('')  # hide label
        ax_var.tick_params(axis="x", labelbottom=False)
        # ax_var.set_xticks([])

    if annotation is not None:
        ax_var = cornertext(ax_var, annotation, size=labelsize, c=linec)

    i_lm = dat.find_lower_mantle()
    # print(dat.name, 'i_lm', i_lm, ', p =', dat.pressure[i_lm] * 1e-9, 'GPa')
    # print('   mass um', dat.mass_um, 'kg', '| mass_h2o um', dat.mass_h2o_um / p.TO, 'TO', '| c_h2o',
    #       dat.mass_h2o_um / dat.mass_um * 1e6, 'ppm')
    axes = (ax_var, ax_mod)
    for ax in axes:
        # ax.set_xlim(0, p_max)
        ax.tick_params(axis='both', which='major', labelsize=ticksize)
        if show_UM:
            ax.axvline(dat.pressure[i_lm] * 1e-9, c='k', lw=10)

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


def find_rich_planets(dir_name, phase, n=1, output_base='output/', get_min=False, **kwargs):
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
    if get_min:
        print('\n', n, 'planets with the lowest', phase)
        print('-------------------------------')
        for val, pl in reversed(yx[-n:]):
            print(pl, '|', val)
        return names_sorted[-n:]
    else:
        print('\n', n, 'planets with the most', phase)
        print('-------------------------------')
        for val, pl in yx[:n]:
            print(pl, '|', val)
    return names_sorted[:n]


def find_extrema(dir_name, attr, get_min=True, get_max=True, n=1, output_base='output/', scale=1, **kwargs):
    """ find n planets with most extreme values """
    if get_min and get_max:
        raise NotImplementedError('find_extrema() can only take min or max at once')

    dats = rw.read_dir(px.perplex_path_default + output_base + dir_name + '/')
    props = []
    names = []
    for dat in dats:
        try:
            props.append(eval('dat.' + attr) * scale)
            names.append(dat.name)  # property not found
        except AttributeError:  # attr not found
            pass
    yx = list(zip(props, names))

    yx.sort(reverse=True)
    names_sorted = [x for y, x in yx]
    # props_sorted = [y for y, x in yx]
    if get_max:
        print('\n', n, 'planets with the highest', attr)
        print('-------------------------------')
        for val, pl in yx[:n]:
            print(pl, '|', val)
        extrema = names_sorted[:n]
    if get_min:
        print('\n', n, 'planets with the lowest', attr)
        print('-------------------------------')
        for val, pl in reversed(yx[-n:]):
            print(pl, '|', val)
        extrema = names_sorted[-n:]
    return extrema
