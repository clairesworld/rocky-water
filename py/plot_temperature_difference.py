import numpy as np
import parameters as p
import main as rw
import perplexdata as px
import plot_perplex as plotpx
import parameters as p
import matplotlib.pyplot as plt
from useful_and_bespoke import cornertext, colorize, colourbar, colourised_legend
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
import matplotlib.transforms as transforms
import matplotlib.patches as patches
from matplotlib.patches import ConnectionPatch

""" difference in sat between 1600 and 1900 """


def plot_sat_difference(Mp=1, pmin=1000, pmax=140e4, fig=None, axes=None, c='k', figsize=(10, 4), save=False,
                        labelsize=16, ticksize=12, legsize=12, title=None, sigma=2, inset=[0.34, 0.2, 0.6, 0.2], hist_kwargs={},
                        dir1=None, dir2=None, show_layers=True):
    from scipy.interpolate import interp1d
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})

    if fig is None:
        fig, axes = plt.subplots(1, 2, figsize=(5, 5))
        # gs = fig.add_gridspec(1, 2,
        #                       width_ratios=(4, 1),
        #                       # left=0.1, right=0.9, bottom=0.1, top=0.9,
        #                       wspace=0.05, hspace=0.2)

        # ax = fig.add_subplot(111)
        ax0_hist = axes[0].inset_axes(inset)
        # ax1_hist = axes[1].inset_axes(inset)ectory_1900

    if dir1 is None and dir2 is None:
        dir1 = px.perplex_path_default + 'output/' + 'hypatia' + str(Mp) + 'M_1600K_80Fe/'
        dir2 = px.perplex_path_default + 'output/' + 'hypatia' + str(Mp) + 'M_1900K_80Fe/'
    datsdict = plotpx.dict_across_pops([dir1, dir2],  # subsample=2,
                                       x_name="df_all['mass_h2o(kg)']", y_name="df_all['P(bar)']")

    # get distribution across stars
    # interpolate pressure
    p_hat = np.linspace(pmin, pmax, 1000)
    x_all = []
    total_mass_diff_all = []
    for key, item in datsdict.items():
        # print(key)
        x = item[0]
        y = item[1]
        try:
            x0 = np.array(x[0])  # 1600 K, water mass profile
            p0 = np.array(y[0])  # 1600 K, pressure profile
            x1 = np.array(x[1])  # 1900 K, water mass profile
            p1 = np.array(y[1])  # 1900 K, pressure profile
            total_mass_diff = (np.sum(x0) - np.sum(x1)) / p.TO
            if total_mass_diff < 0:
                print(key, 'm diff:', total_mass_diff)

            fp0 = interp1d(p0, x0, kind='linear')  # 1600 K
            fp1 = interp1d(p1, x1, kind='linear')  # 1900 K
            try:
                x_hat = fp0(p_hat) / fp1(p_hat)
            except ValueError as e:
                print(e, '| p', np.min(p0), np.max(p0))
                # raise e
            x_all.append(x_hat)
            total_mass_diff_all.append(total_mass_diff)
        except IndexError as e:
            # runs not complete?
            print(e)
            pass

    # get percentiles
    x_all = np.array(x_all)
    # print('x_all', np.shape(x_all))
    if sigma == 1:
        q = [0.16, 0.50, 0.84]
    elif sigma == 2:
        q = [0.0227, 0.50, 0.9773]
    x_min, x_med, x_max = np.quantile(x_all, q, axis=0)

    p_GPa = p_hat * 1e-4

    # first of all, the base transformation of the data points is needed
    # base = plt.gca().transData
    # rot = transforms.Affine2D().rotate_deg(90)  # then add  transform=rot + base to plot function

    ax = axes[0]
    ax.fill_betweenx(p_GPa, x_min, x_max, color=c, alpha=0.2,)
    ax.plot(x_med, p_GPa, c=c, lw=1, )  # median
    ax.set_ylabel('Pressure (GPa)', fontsize=labelsize)
    ax.set_xlabel(r'$\left. w_{\rm cold} \middle/ w_{\rm hot} \right.$', fontsize=labelsize)
    ax.set_xlim(0, 20)
    ax.set_xticks([1, 10, 20])
    ax.set_ylim(np.min(p_GPa), np.max(p_GPa))
    ax.invert_yaxis()

    # add zoom subplot
    ax = axes[1]
    ax.fill_betweenx(p_GPa, x_min, x_max, color=c, alpha=0.2,  zorder=5)
    ax.plot(x_med, p_GPa, c=c, lw=1, zorder=6)  # median
    ax.set_xlabel(r'$\left. w_{\rm cold} \middle/ w_{\rm hot} \right.$', fontsize=labelsize)
    ax.set_xlim(0, 20)
    ax.set_xticks([1, 10, 20])
    ax.set_ylim(np.min(p_GPa), 30)
    ax.yaxis.tick_right()
    ax.invert_yaxis()

    # make hist
    print('mean difference:', np.mean(total_mass_diff_all))
    ax0_hist.hist(total_mass_diff_all, bins=25, range=(0, 2.5), **hist_kwargs)
    # some SiO2 rich planets increase in sat with T because of stish parameterisation increasing with T
    ax0_hist.set_xlabel(r'$\Delta w_{\rm m}$ (OM)', fontsize=labelsize)
    ax0_hist.set_ylim((0, 300))
    ax0_hist.set_yticks([])
    ax0_hist.set_xticks([0, 1, 2])

    # zoom connection
    axes[0].add_artist(patches.Rectangle((0, 00), 20, 30, linewidth=0, edgecolor=None, facecolor='k', alpha=0.1))
    for y in [30, np.min(p_GPa)]:
        axes[1].add_artist(ConnectionPatch(xyA=(0, y), xyB=(20, y), coordsA="data", coordsB="data",
                                           axesA=axes[1], axesB=axes[0], color="0.1", lw=0.1))

    if show_layers:
        # read in again (inefficient)
        dats = rw.read_dir(px.perplex_path_default + 'output/' + 'hypatia' + str(Mp) + 'M_1600K_88Fe/', subsample=1000)
        hline_kwargs = {'lw': 0.1, 'alpha': 0.3, 'c': 'xkcd:clay', 'zorder': 0}
        for dat in dats:
            # only show for ol-bearing
            if 'X_ol' in dat.df_all.columns:
                y_mtz = dat.p_mtz * 1e-9  # GPa
                axes[1].axhline(y_mtz, **hline_kwargs)
                y_lm = dat.p_lm * 1e-9  # GPa
                axes[1].axhline(y_lm, **hline_kwargs)
        axes[1].text(19, 13, 'OBM', fontsize=legsize, va='bottom', ha='right', color=hline_kwargs['c'])
        axes[1].text(19, 22.6, 'MTZ', fontsize=legsize, va='bottom', ha='right', color=hline_kwargs['c'])
        axes[1].text(19, 29.7, 'LM', fontsize=legsize, va='bottom', ha='right', color=hline_kwargs['c'])

    for axe in (axes[0], ax0_hist, axes[1]):
        axe.tick_params(axis='both', which='major', labelsize=ticksize)
    fig.suptitle(title, fontsize=labelsize)

    if save:
        fig.savefig(plotpx.fig_path + 'sat_T_diff' + '.png', bbox_inches='tight', dpi=300)


plot_sat_difference(Mp=1, pmin=1000, pmax=150e4, c='xkcd:green blue', title='',  # 'Mantle cooling',
                    hist_kwargs={'color': '0.9', 'edgecolor': 'k', 'linewidth': 0.5}, sigma=1,
                    # dir1=px.perplex_path_default + 'output/test_compositional_solidus_60Fe',
                    # dir2=px.perplex_path_default + 'output/hypatia1M_1900K_60Fe',
                    save=True
                    )

plt.show()
