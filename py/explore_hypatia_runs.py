import numpy as np
import parameters as p
import ask_hypatia as hyp
import perplexdata as px
import plot_perplex as plotpx
import saturation as sat
import main as rw
import parameters as p
import matplotlib.pyplot as plt
from useful_and_bespoke import cornertext, colorize, colourbar, colourised_legend
from UM_mass_scaling import um_mass_scaling
import saturation as sat
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec

""" get earth benchmark """
# Tp = 1600
# earth = rw.build_planet(M_p=1 * p.M_E, test_CMF=0.325, test_oxides=px.wt_oxides_Earth,
#                         maxIter=30, tol=1e-4, n=800, Tp=Tp,  # core_efficiency=0.8,
#                         plot_all=False, get_saturation=True, verbose=True, clean=True,
#                         vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
#                         name='Earth_' + str(Tp) + 'K',
#                         )
# earth.get_garnet_composition()
# plotpx.single_composition(earth, which='pressure', modality_type='water', comp_stacked=True, save=True,
#                           show=True, cmap='tab20', labelsize=16, plot_phases_order=None, p_max=None, make_legend=True)
# plotpx.single_phase_subfig(earth, 'c_h2o', var_scale=1e6, var_log=True, vertical_pressure=False,
#                                            title=earth.name,
#                                            var_label='Water capacity\n(ppm)', save=False, show=True, phase_order=None,
#                                            labelsize=12, legsize=10, ticksize=12, xpad=10, ypad=10,
#                                            # ax_ticks=[1e1, 1e3, 1e5],
#                                            annotation='total = {0:3.1f} earth oceans'.format(earth.mass_h2o_total / p.TO),
#                                            cmap_phases='tab20', linec='xkcd:navy', linew=2, ymin=1e0, ymax=1e5)

""" check a few cases and plot instantly to verify """
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia2/', subsample=10)
# stars = ['2MASS 19141179+3833548', 'HIP 86087', 'HIP 79126']
# dats_1600 = [rw.read_name(output_path=px.output_parent_default + 'hypatia1M/', star=st, M_p=p.M_E,
#                      core_efficiency=0.8, Tp=1600) for st in stars]
# dats_1900 = [rw.read_name(output_path=px.output_parent_default + 'hypatia1M_1900K/', star=st, M_p=p.M_E,
#                      core_efficiency=0.8, Tp=1900) for st in stars]
# for dat in (dats_1600[0], dats_1900[0]):
#     m_w_tot = dat.mass_h2o_total / p.TO
#     fig, axes = plotpx.composition_subfig(dat, 'c_h2o', var_scale=1e6, var_log=True, vertical_pressure=False,
#                                            title='',
#                                           # phase_order=['gt', 'cpx', 'ol', 'opx', 'hpcpx', 'wad', 'capv', 'ring', 'st', 'wus', 'pv', 'ppv', 'fapv'],
#                                            var_label='Water capacity\n(ppm)', save=True, show=False,
#                                            labelsize=12, legsize=10, ticksize=12, xpad=10, ypad=10,
#                                            # ax_ticks=[1e1, 1e3, 1e5],
#                                            annotation='total = {0:3.1f} earth oceans'.format(m_w_tot),
#                                            cmap_phases='tab20', linec='xkcd:navy', linew=2, ymin=1e0, ymax=1e5)


""" try some crossplots """
# earth.get_phase_masses()

# UM water vs. mg/si
# dats = rw.read_dir(px.perplex_path_default + 'output/hypatia1M/')
# plotpx.pop_scatter(dats, 'mgsi', 'mass_h2o_um', x_scale=1, y_scale=p.TO**-1, c='k', alpha=0.4,
#                         xlabel='Mg/Si', ylabel='Upper mantle water capacity (Earth oceans)', filename=None,
#                         earth=earth, save=True)

# # mass of phase vs. water
# dats = rw.update_dir(px.perplex_path_default + 'output/hypatia2/', px.PerplexData.get_phase_masses, store=False)
# plotpx.pop_scatter(dats[:500], "phase_mass['O']", "mass_h2o_um", x_scale=1, y_scale=p.TO**-1, c='k', alpha=0.4,
#                         ylabel='Upper mantle water capacity (Earth oceans)', xlabel='mass Ol (kg)', filename='mass_h2o_ol',
#                         earth=earth, save=True)

# mgsi effect on upper mantle mass
# dats = rw.update_dir(px.perplex_path_default + 'output/hypatia2M/', px.PerplexData.get_phase_masses, store=False)
# plotpx.pop_scatter(dats, "mgsi", "mass_um", c='xkcd:bordeaux', alpha=0.4,
#                         ylabel='Upper mantle mass (kg)', xlabel='Mg/Si',
#                         earth=earth, save=True)


""" upper mantle water mass - histograms """

def hist_saturation_subfig(which='um', masses=None, xlim='default', bins=None, labelsize=14, cmap='rainbow', save=True,
                           show=True, earth=None, **kwargs):
    if masses is None:
        masses = [0.1, 0.5, 1, 2, 3, 4]
    if which == 'um':
        if xlim == 'default':
            xlim = (0.1, 4.0)
        key = 'mass_h2o_um'
        xlabel = 'Mantle water capacity (Earth oceans)'
        data_label = None
        ls = '-'
        fname = 'histsubplot-Mp_w_um'
    elif which == 'total':
        if xlim == 'default':
            xlim = (0.1, 8)
        key = 'mass_h2o_total'
        xlabel = 'UM water capacity (Earth oceans)'
        data_label = None
        ls = '-'
        fname = 'histsubplot-Mp_w_tot'
    elif which == 'both':
        if xlim == 'default':
            xlim = (0.1, 8)
        key = ['mass_h2o_um', 'mass_h2o_total']
        xlabel = 'Water capacity (Earth oceans)'
        data_label = ['Upper mantle', 'Whole mantle']
        ls = ['-', '--']
        fname = 'histsubplot-Mp_w_all'

    fig, axes = plt.subplots(len(masses), 1)
    c = colorize(masses, cmap)[0]
    for ii, Mp in enumerate(masses):
        if isinstance(Mp, float):
            mass_str = str(Mp).replace('.', ',')
        elif isinstance(Mp, int):
            mass_str = str(Mp)
        else:
            print('MP is', type(Mp))
        planets = rw.read_dir(px.perplex_path_default + 'output/hypatia' + mass_str + 'M/')
        try:
            ax = axes[ii]
        except:
            ax = axes
        fig, ax = plotpx.pop_hist1D(planets, key, scale=p.TO ** -1, earth=None,
                                    xlabel='', c_hist=c[ii], ls=ls,
                                    save=False, show=False, data_label=data_label,
                                    alpha=0.7, xlim=xlim, legsize=10, bins=bins,
                                    histtype='step', fig=fig, ax=ax, labelsize=labelsize, **kwargs)
        # add naive Earth scaling
        if earth is not None:
            w_um_earth = earth.mass_h2o_um / p.TO  # in earth oceans
            print('w UM earth', w_um_earth, 'scaled', w_um_earth * Mp)
        ax.axvline(w_um_earth * Mp, c='k', ls='--', lw=0.5, label=str(Mp) + r' $\times$ Earth value')
        ax.legend(fontsize=10, frameon=False)
        if ii == len(masses) - 1:
            ax.set_xlabel(xlabel)
        else:
            ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(str(Mp) + ' $M_\oplus$', fontsize=labelsize)
    # fig.suptitle('Mantle wantle water capacities, $T_p$ = 1600')
    if save:
        fig.savefig(plotpx.fig_path + fname + '.png')
    if show:
        plt.show()


# hist_saturation_subfig(which='um', masses=[0.1], xlim=(0, 1),
#                        bins=100, showmedian=True, show=False)


def find_percentile(q, dir_name, prop):
    if np.size(q) > 1:
        raise NotImplementedError('find_percentile() - multiple percentiles at once not implemented. use scalar q.')
    dats = rw.read_dir(px.perplex_path_default + 'output/' + dir_name + '/')
    x = []
    names = []
    test_vals = []
    for dat in dats:
        try:
            x.append(eval('dat.' + prop))
            names.append(dat.star)
        except AttributeError:
            pass
    pcen = np.percentile(x, q, method='nearest')
    i_near = abs(x - pcen).argmin()
    print(q, 'pcen', pcen, 'val')
    return names[i_near]


# print(find_percentile(10, 'hypatia1M', 'mgsi'))
# print(find_percentile(50, 'hypatia1M', 'mgsi'))
# print(find_percentile(90, 'hypatia1M', 'mgsi'))

""" c_h2o profiles with T dependence, for an average star"""

# star = '2MASS 19375133+4945541'
# dat = rw.read_name(px.perplex_path_default + 'output/hypatia1M/', star=star, M_p=1*p.M_E, core_efficiency=0.8, Tp=1600)
# fig, ax = plotpx.profile(dat, 'mass_h2o(kg)', independent_ax='pressure', reverse_y=True, ax_label='Water capacity (Earth oceans)', scale=p.TO**-1,
#                          c='xkcd:deep blue', lw=2, alpha=0.9,
#             xmin=None, xmax=None, ymin=0, ymax='max', ax_ticks=None, label_x=True, label_y=True, labelsize=14, legsize=12,
#             fig=None, ax=None, figsize=(4, 6), log=True, leg_label='1600 K', orientation='vertical',
#             y2var=None, y2label=None, y2scale=None, save=False)
#
# dat = rw.read_name(px.perplex_path_default + 'output/hypatia1M_1900K/', star=star, M_p=1*p.M_E, core_efficiency=0.8, Tp=1900)
# fig, ax = plotpx.profile(dat, 'mass_h2o(kg)', independent_ax='pressure', reverse_y=True, ax_label='Water capacity (Earth oceans)', scale=p.TO**-1,
#                          c='xkcd:peach', lw=2, alpha=0.9,
#             xmin=None, xmax=None, ymin=0, ymax='max', ax_ticks=None, label_x=True, label_y=True, labelsize=14, legsize=12,
#             fig=fig, ax=ax, log=True, leg_label='1900 K', orientation='vertical', save=True, fname='massh2o_profiles_hotcold')


""" find n richest planets """
# plotpx.find_rich_planets(dir_name='hypatia3M', phase='st', n=10)
# plotpx.find_rich_planets(dir_name='hypatia3M', phase='seif', n=10)


plt.show()
