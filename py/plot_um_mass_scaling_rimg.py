import numpy as np
import parameters as p
from perplexdata import wt_oxides_Earth, wt_oxides_MD95, perplex_path_default
import plot_perplex as plotpx
import main as rw
import parameters as p
import matplotlib.pyplot as plt
import perplexdata as px
import pickle as pkl
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FuncFormatter
from useful_and_bespoke import cornertext, colorize, colourbar, colourised_legend
from UM_mass_scaling import um_mass_scaling
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec

# plt.rc('text', usetex=True)
# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})

perplex_path_aopp = '/home/g/guimond/Work/perple_x/'
fig_path_aopp = '/home/g/guimond/Work/rocky-water/figs_scratch/'

""" get earth benchmark """
Tp = 1600
earth = rw.build_planet(M_p=1 * p.M_E, test_CMF=0.325, test_oxides=wt_oxides_MD95,
                        maxIter=30, tol=1e-4, n=1200,
                        Tp=Tp,  # core_efficiency=0.8,
                        plot_all=False, get_saturation=False, verbose=True, clean=True,
                        vertex_data='stx21ver', option_file='perplex_option_claire', excluded_phases=[],
                        name='Earth_' + str(Tp) + 'K',
                        perplex_path=perplex_path_aopp, output_parent_path=perplex_path_aopp + 'output/',
                        suppress_output=False,
                        )
# earth = rw.update_single(perplex_path_aopp + 'output/', 'Earth_' + str(Tp) + 'K', px.PerplexData.get_lm_mass, store=True)


""" plot variable across masses with composition continuity """
# x, y = [None, None, None], [None, None, None]
with open(fig_path_aopp + 'mass_scaling_figdata.pkl', "rb") as pfile:
    x, y = pkl.load(pfile)

cmap = 'Oranges'
colours = colorize(np.arange(3), cmap=cmap)[0]
fig, ax = None, None
for ii, fe in enumerate([70, 88, 99]):
    dirs = [perplex_path_aopp + 'output/hypatia_full/hypatia' + s + 'M_1600K_' + str(fe) + 'Fe_hires/' for s in
            ('0,1', '0,3', '0,5', '1', '1,5', '2', '2,5', '3', '4', '5')]

    # for folder in dirs:
    #     dats = rw.update_dir(folder, px.PerplexData.get_lm_mass, store=True)

    xi, yi = x[ii], y[ii]
    fig, ax, xi, yi = plotpx.compare_pop_fillbetween(dirs, 'M_p', 'mass_frac_lm',
                                                   x_scale=p.M_E ** -1, y_scale=100, xlog=False, ylog=False,
                                                   ylabel='Pv- or Ppv-bearing mantle fraction (wt.%)',
                                                   xlabel='Planet mass ($M_\oplus$)',
                                                   sigma=2,
                                                   save=False, show=False, extension='.pdf', xlim=(0, 5), ylim=(0, 100),
                                                   # (0.4, 1.8),
                                                   show_n=False, show_med=False,
                                                   labelsize=10, legsize=10, ticksize=8, #earth=earth,
                                                   legend=False, fig=fig, ax=ax,
                                                   earth_real=None,  # 1.06e24,
                                                   # show_scaling_fn=um_mass_scaling, scalinglabel=r'Constant-$\rho$ scaling',
                                                   c=colours[ii],  patch_kwargs = {'color': colours[ii], 'alpha': 1},
                                                     verbose=True, return_data=True, x=xi, y=yi)
    # x[ii] = xi
    # y[ii] = yi
# with open(fig_path_aopp + "mass_scaling_figdata.pkl", "wb") as pfile:
#     pkl.dump((x, y), pfile)



# Earth
ax.scatter(1, earth.mass_frac_lm * 100, label='Earth', c='#7ef4cc', # marker='$\oplus$', s=200,
           zorder=100)

# Mars too?
# ax.scatter(0.03, 0.03, label='Mars', c='xkcd:dark blue', marker=u'$\u2642$', s=200, zorder=100, transform=ax.transAxes)

# 1900K
x1900, y1900 = None, None
# with open(fig_path_aopp + 'mass_scaling_figdata2.pkl', "rb") as pfile:
#     x, y = pkl.load(pfile)
dirs = [perplex_path_aopp + 'output/Earth_masseffect_1900K/' + str(s) + 'M/' for s in
            ('0,1', '0,3', '0,5', '1', '1,5', '2', '2,5', '3', '4', '5')]
# for folder in dirs:
#     dats = rw.update_dir(folder, px.PerplexData.get_lm_mass, store=True)
fig, ax, x1900, y1900 = plotpx.compare_pop_fillbetween(dirs, 'M_p', 'mass_frac_lm',
                                                 x_scale=p.M_E ** -1, y_scale=100, xlog=False, ylog=False,
                                                 ylabel='Pv- or Ppv-bearing mantle fraction (wt.%)',
                                                 xlabel='Planet mass ($M_\oplus$)',
                                                 sigma=2,
                                                 save=False, show=False, xlim=(0, 5), ylim=(0, 100),
                                                 # (0.4, 1.8),
                                                 show_n=False, show_med=False,
                                                 labelsize=10, legsize=10, ticksize=8,  # earth=earth,
                                                 legend=False, fig=fig, ax=ax,
                                                 earth_real=None,  # 1.06e24,
                                                 # show_scaling_fn=um_mass_scaling, scalinglabel=r'Constant-$\rho$ scaling',
                                                 c='k', patch_kwargs={'color': 'k', 'alpha': 0.5},
                                                 verbose=True, return_data=True, x=x1900, y=y1900)

with open(fig_path_aopp + "mass_scaling_figdata2.pkl", "wb") as pfile:
    pkl.dump((x, y), pfile)

# cbar
dum = np.linspace(.70, .99, 3)
cbaxes = inset_axes(ax, width="30%", height="3%", loc='lower right', borderpad=2)
mappable = ax.scatter(dum, dum, c=dum, cmap=cmap, s=0)
fig.colorbar(mappable=mappable, cax=cbaxes, ticks=[.70,.99], orientation='horizontal', cmap=cmap, )
cbaxes.tick_params(labelsize=8)
cbaxes.set_title(r'Fe$_{core}$/Fe$_{bulk}$', fontsize=10)
cbaxes.xaxis.set_major_formatter(FuncFormatter(lambda q, _: '{:.0%}'.format(q)))

# fig.savefig(fig_path_aopp + 'LM_scaling.pdf', bbox_inches='tight')
plt.show()
