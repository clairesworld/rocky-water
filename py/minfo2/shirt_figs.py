import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter, PercentFormatter, FuncFormatter
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from py.useful_and_bespoke import dark_background

rc('text', usetex=True)

""" 
compare fo2 hist for multiple runs
 need to run these individually and make gif in gimp
 for slides - dark mode
"""

# # set these
z_var = 'X_ferric'
# z_var = 'core_eff'
model = 'melts'
p_of_interest = 4
T_of_interest = 1473

core_eff_constant = 88
X_ferric_constant = 3
core_eff_list = [99, 95, 88, 80, 70]
X_ferric_list = [1, 3, 5, 7]
labelsize = 24
ticksize = 14
legsize = 14
ylim = (0, 1.5)
xlim = (-3.5, 3)
xlabel = r'Upper mantle log$f{\rm O_2}$ relative to Earth' #$\Delta$FMQ'

if model == 'melts':
    source = fo2plt.output_parent_mlt_earth
    z_scale = 100
elif model == 'perplex':
    source = fo2plt.output_parent_px
    z_scale = 1
if z_var == 'X_ferric':
    z_list = X_ferric_list
    dirs = [source + 'hypatia_' + str(core_eff_constant) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
            X_ferric_list]
    vmin, vmax = min(X_ferric_list), max(X_ferric_list)
    legtitle = r'Fe$^{3+}/\Sigma$Fe'
    fname_base = 'hist_fer_dark'
    # xlim = (-3.5, 3)
    ylim = (0, 1.5)
    nbins = 35
    range = None #(-3, 3)
elif z_var == 'core_eff':
    z_list = core_eff_list
    dirs = [source + 'hypatia_' + str(ce) + 'coreeff_' + str(X_ferric_constant) + 'ferric_ext/' for ce in core_eff_list]
    vmin, vmax = min(core_eff_list), max(core_eff_list)
    legtitle = r'Fe$_{\rm core}$/Fe$_{\rm bulk}$'
    fname_base = 'hist_core_dark'
    # xlim = (-2, 2)
    ylim = (0, 1.8)
    nbins = 50
    range = (-2, 1.5)
if p_of_interest == 1:
    x_var = 'delta_qfm_1GPa'  #
    # xlabel = 'log(fO$_2$) at 1 GPa'
elif p_of_interest == 4:
    x_var = 'delta_qfm_4GPa'
    # xlabel = 'log(fO$_2$) at 4 GPa'

z_labels = ['{:.0f}%'.format(z) for z in z_list]  # + r' ($\sigma =$ ' + '{:.2f})'.format(sd),
z_labels[1] = z_labels[1] + ' (Earth)'
for ii, z in enumerate(z_list):
    if ii == 1:  # only 3% highlight
        fig, ax = fo2plt.pop_hist(dirs, x_var=x_var, z_var=z_var, z_scale=z_scale, model=model, xlabel=xlabel,
                                  hilite=ii, figsize=(11, 3), T_of_interest=T_of_interest,
                                  labelsize=labelsize, ticksize=ticksize, legsize=legsize, labelpad=12,
                                  z_labels=z_labels,
                                  legtitle='     ',  legloc='upper left', make_legend=True,
                                  nbins=35, lw=4, cmap='autumn', vmin=vmin, vmax=vmax, alpha=0.9,
                                  exclude_silica=True, save=False, verbose=False, show_sd=False,
                                  ylim=ylim, xlim=xlim, range=range,
                                  # x_shift=-14.431326430547
                                  )

        ax.set_yticks([])
        ax.xaxis.set_minor_locator(AutoMinorLocator())

        # ax.text(0.95, 0.9, '4 GPa', fontsize=labelsize, transform=ax.transAxes, c='xkcd:off white', ha='right', va='top')
        ax.text(0.02, 0.95, legtitle, fontsize=legsize, transform=ax.transAxes, c='xkcd:off white', ha='left', va='top',
                zorder=100)

        ax.set_title(r'$f{\rm O_2}$ modulation by mantle minerals, pMELTS calculations', fontsize=labelsize)

        # fig, ax = dark_background(fig, ax)
        fig.savefig('/home/claire/Works/min-fo2/shirt/images/histogram.png', bbox_inches='tight', facecolor=None,
                    dpi=300)

plt.show()

