import numpy as np
import main as rw
import perplexdata as px
import plot_perplex as plotpx
import matplotlib.pyplot as plt
import pandas as pd
import bulk_composition as bulk
import parameters as p
from useful_and_bespoke import colorize, cornertext

""" histograms of different mg/si distributions for various scenarios 
- Si content of core
- devolatilisation
"""

fig_path = '/home/g/guimond/Work/rocky-water/figs_scratch/'
perplex_path = '/home/g/guimond/Work/perple_x/'
oxide_list = ['MgO', 'SiO2', 'CaO', 'Al2O3', 'FeO']  # 'Na2O

def mgsi_mantle_post_core(core_eff=0.85, core_Si_wtpt=5, oxide_list=oxide_list,
                              path_to_tsv='/home/g/guimond/Work/hypatia-compositions/hypatia-04122023.tsv'):
    # read stellar abundances
    df = pd.read_csv(path_to_tsv, sep='\t', header=0)

    # remove rows with no measurements
    df.replace('', np.nan, inplace=True)  # first convert empty string to nan
    df.dropna(subset=[ox[:2] for ox in oxide_list], inplace=True)

    mgsi_partitioned_all = []
    for ii, star in enumerate([name.rstrip() for name in df.Name.to_numpy()]):
        nH_star = []
        for metal in oxide_list:
            el = metal[:2]
            sol_val = eval('p.' + el.lower() + '_sol')
            nH_star.append(df.iloc[ii][el] + sol_val)

        # get bulk mantle oxides
        wt_oxides_partitioned = bulk.stellar_mantle_corelightelements(oxide_list=oxide_list, nH_star=nH_star, core_eff=core_eff,
                                                     core_Si_wtpt=core_Si_wtpt, verbose=False)
        mgsi_partitioned = (wt_oxides_partitioned['MgO'] / p.M_MgO) / (wt_oxides_partitioned['SiO2'] / p.M_SiO2)

        mgsi_partitioned_all.append(mgsi_partitioned)
    return mgsi_partitioned_all


core_eff = 0.85
wt_Si_cores = [0, 5, 10]
colour = colorize(np.arange(len(wt_Si_cores)), cmap='copper')[0]
xmin, xmax = (0.2, 2.25)
fig, ax = plt.subplots(1,1, figsize=(6, 3))
labelsize = 14
ticksize = 9
legsize = 9
bins = 70

for ii, si in enumerate(wt_Si_cores):
    print('calculating', si, 'wt% Si')
    mgsi = mgsi_mantle_post_core(core_eff=core_eff, core_Si_wtpt=si)
    ax.hist(mgsi, color=colour[ii], label=str(si), lw=2, bins=bins, histtype='step', range=(xmin, xmax), density=True)
    ax.axvline(np.median(mgsi), c=colour[ii], lw=1, ls=':')
    ax.hist(mgsi, facecolor=colour[ii], alpha=0.5, lw=0, bins=bins, histtype='bar', range=(xmin, xmax), zorder=0, density=True)


ax.text(0.98, 0.05, 'n = ' + str(len(mgsi)), fontsize=legsize, transform=ax.transAxes, va='bottom', ha='right')
ax.set_xlabel('Mantle Mg/Si', fontsize=labelsize)
ax.set_ylabel('')
ax.set_xlim((xmin, xmax))
legend = ax.legend(frameon=False, fontsize=legsize, title='Core Si (wt.%)', loc='center left')
plt.setp(legend.get_title(),fontsize=legsize)
ax.tick_params(axis='both', which='major', labelsize=ticksize)
ax.set_yticks([])

# add pyrolitic range [0.67608298 0.83176377] - [1.44543977 1.77827941]
# ax.annotate('', xy=(0.83178, 2.6), xycoords='data',
#              xytext=(1.4454, 2.6), textcoords='data',
#              arrowprops=dict(facecolor='black', width=1, arrowstyle='<->'))
import matplotlib.patches as patches
y = 3.3
height = 0.2
c_arrow = '0.5'
ax.plot([0.83176377, 1.44543977], [y, y], c=c_arrow, zorder=10)
ax.add_patch(patches.Rectangle((1.44543977, y - (height/2)), (1.77827941 - 1.44543977), height, linewidth=1, hatch='//', edgecolor=c_arrow, facecolor='w'))
ax.add_patch(patches.Rectangle((0.67608298, y - (height/2)), (0.83176377 - 0.67608298), height, linewidth=1, hatch='//', edgecolor=c_arrow, facecolor='w'))
ax.text(1.09, 3.17, 'pyroxene-olivine mantles', fontsize=ticksize, va='top', ha='left', c=c_arrow)

# add sun/BSE
y_ref = 0.2
ax.scatter((10**p.mg_sol / 10**p.si_sol), y_ref, c='#feff00')
ax.scatter(1.25, y_ref, c='#7ef4cc')  # '#00abdc'

print('done')
plt.savefig(fig_path + 'mgsi_histograms_' + str(core_eff*100).replace('.', ',') + 'Fecore.pdf', bbox_inches='tight')
plt.show()
