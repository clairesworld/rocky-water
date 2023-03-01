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
import datetime

# for running remotely
import matplotlib

where = 'starlite'
ftype = '.pdf'

date = datetime.date.today()
print('date', str(date))

fig_path_starlite = '/home/claire/Works/min-fo2/figs_scratch/'

if where == 'apollo':
    opp_mlt = mfug.output_parent_default
    opp_perplex = pfug.output_parent_apollo
    perplex_path = '/raid1/cmg76/perple_x/'
    fig_path = '/raid1/cmg76/alphamelts/figs/'
    matplotlib.use('Agg')
elif where == 'starlite':
    opp_mlt = mfug.output_parent_mlt_earth
    opp_perplex = pfug.output_parent_default
    perplex_path = '/home/claire/Works/perple_x/'
    fig_path = fo2plt.figpath

""" scatter plot of where Fe3+ is on ternary diagram """


def ternary_scatter(p_of_interest=None, T_of_interest=None, core_eff=88, Xf=3.0, component='Fe2O3', z_var='mgsi',
                    z_label=None, model='melts', cmap='rainbow', vmin=None, vmax=None, fontsize=16, offset=0.1,
                    phases=['Opx', 'CPx', 'Sp'], marker='o', title=None, ftype=ftype,
                    absolute_abundance=True,
                    name=None, opp=None, save=False, mec=None, lw=1.5, ticksize=12, fig_path=fo2plt.figpath,
                    fig=None, tax=None, ax=None, verbose=False, **kwargs):
    """ plot distribution of component between 3 phases
    p_of_interest: GPa
    """
    missing_Sp = 0
    missing_Opx = 0
    missing_Cpx = 0
    missing_Gt = 0
    tracker = missing_Sp, missing_Opx, missing_Cpx, missing_Gt

    if model == 'perplex':
        phases_px = phases
    else:
        phases_px = [mfug.map_to_px_phase[k] for k in phases]  #

    if tax is None:
        fig, tax = ternary.figure(scale=100, ax=ax)
        # plot setup
        tax.bottom_axis_label(phases_px[0], fontsize=fontsize, offset=0)
        tax.right_axis_label(phases_px[1], fontsize=fontsize, offset=offset)
        tax.left_axis_label(phases_px[2], fontsize=fontsize, offset=offset)
        tax.set_title(title, fontsize=fontsize)

        # Draw Boundary and Gridlines
        tax.boundary(linewidth=2.0)
        tax.gridlines(color="blue", multiple=5)
        tax.get_axes().axis('off')
        tax.clear_matplotlib_ticks()
        tax.ticks(axis='lbr', linewidth=1, multiple=10, tick_formats="%.0f", fontsize=ticksize)

        # add colourbar
        # if z_var:
        # cbar = colourbar(mappable=None, vector=[vmin, vmax], ax=plt.gca(), vmin=vmin, vmax=vmax, label=z_label,
        #                  labelsize=fontsize,
        #                  ticksize=ticksize, labelpad=17, loc='right', cmap=cmap, c='k', pad=0.1, shrink=0.5,
        #                  size="3%")

    # fig.set_size_inches(10, 10)
    tax.set_background_color(color='w', zorder=-1000, alpha=0)

    def draw_point(name, z_var, p_of_interest, tracker):
        missing_Sp, missing_Opx, missing_Cpx, missing_Gt = tracker
        if model == 'perplex':
            dat = pfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp, perplex_path=perplex_path,
                                         load_results_csv=True, **kwargs)

        elif model == 'melts':
            dat = mfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp,
                                         load_results_csv=True, **kwargs)
            p_of_interest = p_of_interest * 1e4

        if dat is None:
            return None, tracker
        if 'X_q' in dat.data.columns:  # don't print quartz saturation cases
            print('dropping case with quartz:', dat.name)
            return None, tracker

        # print(dat.name, '\n', dat.data.head())

        if z_var:
            z = [eval('dat.' + z_var)]
        else:
            z = 'k'

        # d = dat.read_phase_main_components(p_of_interest=p_of_interest, T_of_interest=T_of_interest,
        #                                    component=component, absolute_abundance=absolute_abundance, verbose=False)
        d = dat.get_phase_composition_dict(p_of_interest=p_of_interest, T_of_interest=T_of_interest, component=component,
                                           phases=phases, to_absolute_abundance=absolute_abundance, verbose=True)

        if d is None:
            # general failure for this case and p of interest
            return None, tracker

        # get opx, cpx, sp or gt, normalise to 100
        d2 = {x: np.round(y, 3) for x, y in d.items()}  # if (y != 0) and (~np.isnan(y))}  # round to remove fp imprecision
        if len(d2) > 3:
            raise NotImplementedError(name, 'more than 3 ferric hosts:', d2)
        if sum(d2.values()) == 0:
            print(dat.name, d2, 'sum = 0')
            return None, tracker  # exit to outer function
        factor = 100 / sum(d2.values())
        for k in d2:
            d2[k] = d2[k] * factor
        # print('d2', d2, 'sum', sum(d2.values()))

        xyz = []
        for k in phases_px:  # preserve order
            try:
                if np.isnan(d2[k]):
                    if k == 'Sp':
                        missing_Sp += 1
                    elif k == 'Opx':
                        missing_Opx += 1
                    elif k == 'Cpx':
                        missing_Cpx += 1
                    elif k == 'Gt':
                        missing_Gt += 1

                xyz.append(d2[k])
            except KeyError:
                xyz.append(0)  # e. g. no spinel in this composition?

        # add point to axes
        xyz = tuple(xyz)
        tax.scatter([xyz], marker=marker, edgecolors=mec, linewidths=lw, c=z, s=90, alpha=0.4,
                    cmap=cmap, vmin=vmin, vmax=vmax, zorder=100)
        tracker = missing_Sp, missing_Opx, missing_Cpx, missing_Gt
        return 1, tracker

    if name is not None:
        draw_point(name, z_var, p_of_interest, tracker)
    else:
        count = 0
        if model == 'melts':
            opp = opp_mlt + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'
        elif model == 'perplex':
            opp = opp_perplex + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'
        subfolders = rw.get_run_dirs(output_path=opp)
        for ii, sub in enumerate(subfolders):
            if len(os.listdir(sub)) > 1:  # 1 if contains nH_star e.g.
                name = os.path.basename(sub)
                add, tracker = draw_point(name, z_var, p_of_interest, tracker)
                if add is not None:
                    count += add

    if save:
        plt.savefig(fig_path + 'ferric_ternary_' + str(p_of_interest) + 'GPa' + ftype, bbox_inches='tight')

    missing_Sp, missing_Opx, missing_Cpx, missing_Gt = tracker
    print('num points plotted', count)
    print('missing Sp:', missing_Sp)
    print('missing Opx:', missing_Opx)
    print('missing Cpx:', missing_Cpx)
    print('missing Gt:', missing_Gt)
    return fig, tax


phases_1GPa = ['orthopyroxene', 'clinopyroxene', 'spinel']
phases_4GPa = ['orthopyroxene', 'clinopyroxene', 'garnet']
fig, axes = plt.subplots(2, 2, figsize=(20, 20))
_, tax = ternary_scatter(p_of_interest=1, T_of_interest=1373.15, core_eff=88, Xf=3.0, component='Fe2O3', z_var='mgsi',
                         model='melts', cmap='viridis', vmin=0.69, vmax=1.6,
                         phases=phases_1GPa,
                         z_label='Mg/Si', mec='xkcd:midnight blue', lw=1.5, title='MELTS',
                         ax=axes[0][0], save=False)

# _, tax = ternary_scatter(p_of_interest=1, T_of_interest=1373, core_eff=88, Xf=3.0, component='Fe2O3', z_var='mgsi',
#                          model='perplex', cmap='viridis', vmin=0.69, vmax=1.6, phases=[mfug.map_to_px_phase[ph] for ph in phases_1GPa],
#                          mec='xkcd:brick red', marker='d', lw=1.5, title='Perple_x',
#                          ax=axes[0][1], save=False)

_, tax = ternary_scatter(p_of_interest=4, T_of_interest=1373.15, core_eff=88, Xf=3.0, component='Fe2O3', z_var='mgsi',
                         model='melts', cmap='viridis', vmin=0.69, vmax=1.6,
                         phases=phases_1GPa,
                         z_label='Mg/Si', mec='xkcd:midnight blue', lw=1.5,
                         ax=axes[1][0], save=False)

# _, tax = ternary_scatter(p_of_interest=3.9, T_of_interest=1373, core_eff=88, Xf=3.0, component='Fe2O3', z_var='mgsi',
#                          model='perplex', cmap='viridis', vmin=0.69, vmax=1.6, phases=[mfug.map_to_px_phase[ph] for ph in phases_4GPa],
#                          mec='xkcd:brick red', marker='d', lw=1.5,
#                          ax=axes[1][1], save=False)

fig.suptitle("Fe$^{3+}$ modality", fontsize=16)
fig.savefig(fig_path + 'ternary_subplots' + date + ftype)
plt.show()


""" test individual """
# _, tax = ternary_scatter(p_of_interest=4, T_of_interest=1373.15, core_eff=88, Xf=7.0, component='Fe2O3', z_var='mgsi',
#                     model='melts', cmap='viridis', vmin=0.69, vmax=1.6, phases=['orthopyroxene', 'clinopyroxene', 'garnet'],
#                    mec='r', marker='o', lw=1.5,
#                 save=True, fig_path=fig_path'#
#  )
