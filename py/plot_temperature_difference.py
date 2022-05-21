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

""" difference in sat between 1600 and 1900 """
def plot_sat_difference(Mp=1, pmin=1000, pmax=140e4, fig=None, ax=None, c='k', figsize=(10, 4), save=True,
                        labelsize=14, title=None, hist_kwargs={}):
    from scipy.interpolate import interp1d

    if fig is None:
        fig = plt.figure(figsize=(2.5, 5))
        # gs = fig.add_gridspec(1, 2,
        #                       width_ratios=(4, 1),
        #                       # left=0.1, right=0.9, bottom=0.1, top=0.9,
        #                       wspace=0.05, hspace=0.2)

        ax = fig.add_subplot(111)
        ax_hist = ax.inset_axes([0.4, 0.15, 0.4, 0.2])  # fig.add_subplot(gs[1])
        # fig, (ax, ax_hist) = plt.subplots(1, 2, figsize=figsize)

    directory_1600 = px.perplex_path_default + 'output/' + 'hypatia' + str(Mp) + 'M/'
    directory_1900 = px.perplex_path_default + 'output/' + 'hypatia' + str(Mp) + 'M_1900K/'
    datsdict = plotpx.dict_across_pops([directory_1600, directory_1900], #subsample=10,
                                       x_name="df_all['mass_h2o(kg)']", y_name="df_all['P(bar)']")

    # print('datsdict\n', datsdict)

    # test
    # figs, axes2 = plt.subplots(3, 1)
    # dat_tup = datsdict[list(datsdict.keys())[0]]
    # dat_tupx = dat_tup[0]
    # dat_tupy = dat_tup[1]
    # axes2[0].plot(dat_tupy[0] * 1e-4, dat_tupx[0] / p.TO, c=c, lw=2)  #  test 1600 K
    # axes2[0].set_title('1600 K')
    # axes2[1].plot(dat_tupy[1] * 1e-4, dat_tupx[1] / p.TO, c=c, lw=2)  #  test 1600 K
    # axes2[1].set_title('1900 K')
    # axes2[2].plot(dat_tupy[1] * 1e-4, dat_tupx[0]/dat_tupx[1], c=c, lw=2)
    # axes2[0].set_ylabel('H2O mass (TO)')
    # axes2[1].set_ylabel('H2O mass (TO)')
    # axes2[2].set_ylabel('H2O mass (ratio)')
    # axes2[-1].set_xlabel('Pressure (GPa)')

    # get distribution across stars
    # interpolate pressure
    p_hat = np.linspace(pmin, pmax, 1000)
    x_all = []
    total_mass_diff_all = []
    for key, item in datsdict.items():
        x = item[0]
        y = item[1]
        try:
            x0 = np.array(x[0])  # 1600 K, water mass profile
            p0 = np.array(y[0])  # 1600 K, pressure profile
            x1 = np.array(x[1])  # 1900 K, water mass profile
            p1 = np.array(y[1])  # 1900 K, pressure profile
            total_mass_diff = (np.sum(x0) - np.sum(x1))/p.TO
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
    x_min, x_med, x_max = np.quantile(x_all, [0.0227, 0.50, 0.9773],  #[0.16, 0.50, 0.84],
                                      axis=0)

    p_GPa = p_hat * 1e-4

    # first of all, the base transformation of the data points is needed
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(90)  # then add  transform=rot + base to plot function

    ax.fill_betweenx(p_GPa, x_min, x_max, color=c, alpha=0.2,)
    ax.plot(x_med, p_GPa, c=c, lw=1)  # median
    ax.set_ylabel('Pressure (GPa)', fontsize=labelsize)
    ax.set_xlabel('H$_2$O mass ratio', fontsize=labelsize)
    ax.set_ylim(np.min(p_GPa), np.max(p_GPa))
    ax.invert_yaxis()

    # make hist
    ax_hist.hist(total_mass_diff_all, bins=25, range=(0, 2.5), **hist_kwargs)
    # some SiO2 rich planets increase in sat with T because of stish parameterisation increasing with T
    ax_hist.set_xlabel('Total difference (OM)', fontsize=labelsize-4)
    ax_hist.set_ylim((0, 300))
    ax_hist.set_yticks([])
    ax_hist.set_xticks([0, 1, 2])

    ax.set_title(title, fontsize=labelsize)

    if save:
        fig.savefig(plotpx.fig_path + 'sat_T_diff' + '.png', bbox_inches='tight')


plot_sat_difference(Mp=1, pmin=1000, pmax=150e4, c='xkcd:green blue', title='Mantle cooling',
                    hist_kwargs={'color': '0.9', 'edgecolor': 'k', 'linewidth': 0.5})

plt.show()