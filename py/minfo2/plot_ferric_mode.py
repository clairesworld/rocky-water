import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(os.path.dirname(PARENT_DIR))

import ternary
import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
from py.useful_and_bespoke import dark_background
import numpy as np
import meltsfugacitydata as mfug
import perplexfugacitydata as pfug
from py import main as rw
from py.useful_and_bespoke import colourbar

# for running remotely
import matplotlib
matplotlib.use('Agg')

opp_mlt = '/raid1/cmg76/alphamelts/output/rocky-fo2/earth-tea23/'
# opp_mlt = '/home/claire/Works/min-fo2/alphamelts_output/earth-tea23/'

""" scatter plot of where Fe3+ is on ternary diagram """


def ternary_scatter(p_of_interest=None, T_of_interest=None, core_eff=88, Xf=3.0, component='Fe2O3', z_var='mgsi',
                    z_label=None, model='melts', cmap='rainbow', vmin=None, vmax=None, fontsize=14, offset=0.1, phases=['orthopyroxene', 'clinopyroxene', 'spinel'],
                    name=None, opp=None, save=False, ticksize=12, fig_path=fo2plt.figpath, **kwargs):
    """ plot distribution of component between 3 phases """

    fig, tax = ternary.figure(scale=100)
    fig.set_size_inches(10, 10)
    tax.bottom_axis_label(phases[0], fontsize=fontsize, offset=0)
    tax.right_axis_label(phases[1], fontsize=fontsize, offset=offset)
    tax.left_axis_label(phases[2], fontsize=fontsize, offset=offset)
    tax.set_title("Fe$^{3+}$ modality", fontsize=fontsize)

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=5)
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    tax.ticks(axis='lbr', linewidth=1, multiple=10, tick_formats="%.0f")

    def draw_point(name, z_var):
        if model == 'perplex':
            dat = pfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp,
                                         load_results_csv=True, **kwargs)
        elif model == 'melts':
            dat = mfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp,
                                         load_results_csv=True, **kwargs)
        if dat is None:
            return None
        if 'X_q' in dat.data.columns:  # don't print quartz saturation cases
            print('dropping case with quartz')
            return None

        if z_var is not None:
            z = [eval('dat.' + z_var)]

            cbar = colourbar(mappable=None, vector=[vmin, vmax], ax=plt.gca(), vmin=vmin, vmax=vmax, label=z_label, labelsize=fontsize,
                          ticksize=ticksize, labelpad=17, loc='right', cmap=cmap, c='k', pad=0.05)
        else:
            z = 'k'
        # print('z', z)

        d = dat.read_phase_comp(p_of_interest=p_of_interest, T_of_interest=T_of_interest, component=component,
                                verbose=False)

        if d is None:
            return None
        # get opx, cpx, sp or gt, normalise to 100
        d2 = {x:np.round(y, 3) for x,y in d.items() if (y!=0) and (~np.isnan(y))}  # round to remove fp imprecision
        if len(d2) > 3:
            raise NotImplementedError(name, 'more than 3 ferric hosts:', d2)
        if sum(d2.values()) == 0:
            print('name', d2, 'sum = 0')
            return None  # exit to outer function
        factor = 100 / sum(d2.values())
        for k in d2:
            d2[k] = d2[k] * factor
        print('d2', d2, 'sum', sum(d2.values()))

        xyz = []
        for k in phases:  # preserve order
            try:
                xyz.append(d2[k])
            except KeyError:
                xyz.append(0)  # e. g. no spinel in this composition?
        xyz = tuple(xyz)
        tax.scatter([xyz], marker='o', c=z, s=60, alpha=0.2, cmap=cmap, vmin=vmin, vmax=vmax, zorder=100)
        return 1

    if name is not None:
        draw_point(name, z_var)
    else:
        count = 0
        if model == 'melts':
            opp = opp_mlt + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'
        subfolders = rw.get_run_dirs(output_path=opp)
        for ii, sub in enumerate(subfolders):
            if len(os.listdir(sub)) > 1:  # 1 if contains nH_star e.g.
                name = os.path.basename(sub)
                add = draw_point(name, z_var)
                if add is not None:
                    count += add

    if save:
        plt.savefig(fig_path + 'ferric_ternary_' + str(p_of_interest) + 'GPa.png', bbox_inches='tight')
    print('num points', count)


ternary_scatter(p_of_interest=1, T_of_interest=1373.15, core_eff=88, Xf=3.0, component='Fe2O3', z_var='mgsi',
                    model='melts', cmap='viridis', vmin=0.5, vmax=1.6, phases=['orthopyroxene', 'clinopyroxene', 'spinel'],
                z_label='Mg/Si',
                    # name='Stolper', opp=mfug.output_parent_default
                save=True,fig_path='/raid1/cmg76/alphamelts/figs/')

# plt.show()
