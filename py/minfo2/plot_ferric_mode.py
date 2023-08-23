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
from py.perplexdata import perplex_path_default
from py import main as rw
from py.useful_and_bespoke import colourbar
import datetime
import matplotlib  # for running remotely
from matplotlib import rc
import pickle
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmcrameri

""" scatter plot of where Fe3+ is on ternary diagram """

where = 'starlite'
ftype = '.pdf'

date = datetime.date.today()
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

phases_1GPa = ['orthopyroxene', 'clinopyroxene', 'spinel']
phases_4GPa = ['orthopyroxene', 'clinopyroxene', 'garnet']

if where == 'apollo':
    opp_mlt = mfug.output_parent_default
    opp_perplex = pfug.output_parent_apollo
    perplex_path = '/raid1/cmg76/perple_x/'
    fig_path = '/raid1/cmg76/alphamelts/figs/'
    matplotlib.use('Agg')
elif where == 'starlite':
    opp_mlt = mfug.output_parent_mlt_earth
    opp_perplex = pfug.output_parent_default + 'hypatia/'
    perplex_path = perplex_path_default
    fig_path = fo2plt.figpath


def ternary_scatter(p_of_interest=None, T_of_interest=None, core_eff=88, Xf=3.0, component='Fe2O3', v_var='mgsi',
                    model='melts', phases=['Opx', 'CPx', 'Sp'], cmap='rainbow', vmin=None, vmax=None,
                    marker='o', title=None, titlepad=50, fontsize=16, xoffset=0.1, offset=0.15,
                    markersize=110, markeralpha=0.6, show_cmap=False,
                    absolute_abundance=True,
                    name=None, opp=None, save=False, mec=None, lw=1.5, ticksize=12, fig_path=fo2plt.figpath, fname=None,
                    fig=None, tax=None, ax=None, verbose=False, pickle_to=None, pickle_from=None, ftype=ftype,
                    **kwargs):
    """ plot distribution of component between 3 phases
    p_of_interest: GPa
    """
    print('Starting ternary plot with phases', phases)

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
        _, tax = ternary.figure(scale=100, ax=ax)
        # plot setup
        tax.bottom_axis_label(phases_px[0], fontsize=fontsize, offset=xoffset)
        tax.right_axis_label(phases_px[1], fontsize=fontsize, offset=offset)
        tax.left_axis_label(phases_px[2], fontsize=fontsize, offset=offset)
        tax.set_title(title, fontsize=fontsize, pad=titlepad)
        # print('ax labels:', phases_px[0], phases_px[1], phases_px[2])

        # Draw Boundary and Gridlines
        tax.boundary(linewidth=2.0)
        tax.gridlines(color="k", multiple=5)
        tax.get_axes().axis('off')
        tax.clear_matplotlib_ticks()
        tax.ticks(axis='lbr', linewidth=1, multiple=20, tick_formats="%.0f", fontsize=ticksize, offset=0.02)

        # add colourbar
        # if z_var:
        # cbar = colourbar(mappable=None, vector=[vmin, vmax], ax=plt.gca(), vmin=vmin, vmax=vmax, label=z_label,
        #                  labelsize=fontsize,
        #                  ticksize=ticksize, labelpad=17, loc='right', cmap=cmap, c='k', pad=0.1, shrink=0.5,
        #                  size="3%")

    # fig.set_size_inches(10, 10)
    tax.set_background_color(color='w', zorder=-1000, alpha=0)

    def draw_point(name, v_var, p_of_interest, tracker, xyz0=None, v0=None):
        if xyz0 is not None:
            tax.scatter([xyz0], marker=marker, edgecolors=mec, linewidths=lw, c=v0, s=markersize, alpha=markeralpha,
                        cmap=cmap, vmin=vmin, vmax=vmax, zorder=100, colorbar=show_cmap)
            return 1, xyz0, v0, None
        else:
            print('xyz0', xyz0)
            missing_Sp, missing_Opx, missing_Cpx, missing_Gt = tracker
            if model == 'perplex':
                dat = pfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp, perplex_path=perplex_path,
                                             load_results_csv=True, verbose=verbose, **kwargs)
                phases0 = [pfug.map_to_JH_phase[ph] for ph in phases]

            elif model == 'melts':
                dat = mfug.init_from_results(name, X_ferric=Xf, output_parent_path=opp,
                                             load_results_csv=True, verbose=verbose, **kwargs)
                p_of_interest = p_of_interest * 1e4
                phases0 = phases

            if dat is None:
                return None, None, None, tracker
            if 'X_q' in dat.data.columns:  # don't print quartz saturation cases
                if verbose:
                    print('dropping case with quartz:', dat.name)
                return None, None, None, tracker

            # print(dat.name, '\n', dat.data.head())

            if v_var:
                v = [eval('dat.' + v_var)]
            else:
                v = 'k'  # hi

            d = dat.get_phase_composition_dict(p_of_interest=p_of_interest, T_of_interest=T_of_interest,
                                               component=component,
                                               phases=phases0, to_absolute_abundance=absolute_abundance,
                                               verbose=verbose)

            if d is None:
                # general failure for this case and p of interest
                return None, None, None, tracker

            # get opx, cpx, sp or gt, normalise to 100
            d2 = {x: np.round(y, 3) for x, y in
                  d.items()}  # if (y != 0) and (~np.isnan(y))}  # round to remove fp imprecision
            if len(d2) > 3:
                raise NotImplementedError(name, 'more than 3 ferric hosts:', d2)
            if sum(d2.values()) == 0:
                print(dat.name, d2, 'sum = 0')
                return None, None, None, tracker  # exit to outer function
            factor = 100 / sum(d2.values())
            for k in d2:
                d2[k] = d2[k] * factor

            xyz = []
            for k in phases_px:  # preserve order
                try:
                    xyz.append(d2[k])
                    if d2[k] == 0:
                        if k == 'Sp':
                            missing_Sp += 1
                        elif k == 'Opx':
                            missing_Opx += 1
                        elif k == 'Cpx':
                            missing_Cpx += 1
                        elif k == 'Gt':
                            missing_Gt += 1
                except KeyError:
                    xyz.append(0)  # e. g. no spinel in this composition?

            if np.nan in xyz:
                print('xyz', xyz, dat.name)

            # add point to axes
            xyz = tuple(xyz)
            tax.scatter([xyz], marker=marker, edgecolors=mec, linewidths=lw, c=v, s=markersize, alpha=markeralpha,
                        cmap=cmap, vmin=vmin, vmax=vmax, zorder=100, colorbar=show_cmap)
            tracker = missing_Sp, missing_Opx, missing_Cpx, missing_Gt
            # print('missing counts (Sp, Opx, Cpx, Gt):', tracker )
        return 1, xyz, v, tracker

    if name is not None:
        draw_point(name, v_var, p_of_interest, tracker)
    else:
        if pickle_from is None:
            xyz_list, v_list = [], []
            count = 0
            if model == 'melts':
                opp = opp_mlt + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'
            elif model == 'perplex':
                opp = opp_perplex + 'hypatia_' + str(core_eff) + 'coreeff_' + str(int(Xf)) + 'ferric_ext/'
            subfolders = rw.get_run_dirs(output_path=opp)
            for ii, sub in enumerate(subfolders):
                if len(os.listdir(sub)) > 1:  # 1 if contains nH_star e.g.
                    name = os.path.basename(sub)
                    add, xyzi, vi, tracker = draw_point(name, v_var, p_of_interest, tracker, xyz0=None, v0=None)
                    if add is not None:
                        count += add
                    xyz_list.append(xyzi)
                    v_list.append(vi)

            missing_Sp, missing_Opx, missing_Cpx, missing_Gt = tracker
            print('num points plotted', count)
            print('missing Sp:', missing_Sp)
            print('missing Opx:', missing_Opx)
            print('missing Cpx:', missing_Cpx)
            print('missing Gt:', missing_Gt)

            if pickle_to is not None:
                pickle.dump((xyz_list, v_list), open(pickle_to, "wb"))

        else:
            xyz_list, v_list = pickle.load(open(pickle_from, "rb"))
            for (xyz0, v0) in zip(xyz_list, v_list):
                print('xyz0', xyz0)
                if xyz0 is not None:
                    draw_point(name=None, v_var=None, p_of_interest=None, tracker=None, xyz0=xyz0, v0=v0)

    if fname is None:
        fname = 'ferric_ternary_' + str(p_of_interest) + 'GPa'
    if save:
        plt.savefig(fig_path + fname + ftype, bbox_inches='tight')
    return fig, tax


""" make 4x4 figure with perplex and melts columns, 1 and 4 GPa rows """

# fontsize = 36  # for matplotlib default
xoffset, offset = 0.1, 0.15  # for matplotlib default
fontsize = 42  # for latex sansserif
ticksize = 22  # for latex sansserif
titlepad = 70
# fig, axes = plt.subplots(2, 2, figsize=(20, 20))

# set up gridspec
fig = plt.figure(figsize=(20, 20))
# gs0 = gridspec.GridSpec(1, 2, width_ratios=(30, 1), wspace=0.2, figure=fig)
# gs00 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[0], height_ratios=(1, 1), width_ratios=(1, 1))
gs0 = gridspec.GridSpec(2, 2, width_ratios=(1, 1), height_ratios=(1, 1), wspace=0.06, hspace=0.15, figure=fig)
ax00 = fig.add_subplot(gs0[0, 0])
ax10 = fig.add_subplot(gs0[1, 0])
ax01 = fig.add_subplot(gs0[0, 1])
ax11 = fig.add_subplot(gs0[1, 1])

vmin, vmax = 0.8, 1.7
mec = 'xkcd:midnight blue'
cmap = cmcrameri.cm.bamako_r

# alphamelts column
ax00.text(0.05, 0.9, '1 GPa', transform=ax00.transAxes, size=fontsize, )
fig, tax0 = ternary_scatter(p_of_interest=1, T_of_interest=1373.15, core_eff=88, Xf=3.0, component='Fe2O3',
                            v_var='mgsi', model='melts', cmap=cmap, vmin=vmin, vmax=vmax,
                            fontsize=fontsize, xoffset=xoffset, offset=offset, ticksize=ticksize,
                            phases=phases_1GPa, title='Fe$_2$O$_3$ distribution in pMELTS', titlepad=titlepad,
                            show_cmap=False, save=True, mec=mec, lw=1.5, fname='ternary-tmp', fig=fig,
                            ax=ax00, pickle_from=fig_path + 'data/ternary0.p')

ax10.text(0.05, 0.9, '4 GPa', transform=ax10.transAxes, size=fontsize, )
fig, tax1 = ternary_scatter(p_of_interest=4, T_of_interest=1373.15, core_eff=88, Xf=3.0, component='Fe2O3',
                            v_var='mgsi', model='melts', cmap=cmap, vmin=vmin, vmax=vmax,
                            fontsize=fontsize, xoffset=xoffset, offset=offset, ticksize=ticksize,
                            phases=phases_4GPa, show_cmap=False, save=True, mec=mec, lw=1.5,
                            fname='ternary-tmp', fig=fig,
                            ax=ax10, pickle_from=fig_path + 'data/ternary1.p')

# perplex column
fig, tax2 = ternary_scatter(p_of_interest=1, T_of_interest=1373, core_eff=88, Xf=3.0, component='Fe2O3', v_var='mgsi',
                            model='perplex', cmap=cmap, vmin=vmin, vmax=vmax,
                            fontsize=fontsize, xoffset=xoffset, offset=offset, ticksize=ticksize,
                            phases=[mfug.map_to_px_phase[ph] for ph in phases_1GPa], marker='s',
                            title='Fe$_2$O$_3$ distribution in Perple_X',  titlepad=titlepad,
                            show_cmap=False, save=True, mec=mec, lw=1.5,
                            fname='ternary-tmp', fig=fig, ax=ax01,
                            pickle_from=fig_path + 'data/ternary2.p')

fig, tax3 = ternary_scatter(p_of_interest=3.9, T_of_interest=1373, core_eff=88, Xf=3.0, component='Fe2O3', v_var='mgsi',
                            model='perplex', cmap=cmap, vmin=vmin, vmax=vmax,
                            fontsize=fontsize, xoffset=xoffset, offset=offset, ticksize=ticksize,
                            phases=[mfug.map_to_px_phase[ph] for ph in phases_4GPa], marker='s', show_cmap=False,
                            save=True, mec=mec, lw=1.5, fname='ternary-tmp', fig=fig, ax=ax11,
                            pickle_from=fig_path + 'data/ternary3.p')

# fig.suptitle("Fe$^{3+}$ modality", fontsize=fontsize, y=0.9)

# # adjust spacing between subplots?? might be superfluous given call to gridspec
# plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.06, hspace=0.15)

# add colorbar
dum = np.linspace(vmin, vmax)
im = ax00.scatter(dum, dum, c=dum, cmap=cmap, s=0, vmin=vmin, vmax=vmax)

# cax = fig.add_subplot(gs0[0, 1])
# cbar = fig.colorbar(im, cax=cax, fraction=0.046,
#                     format="%.1f")
# cbar.ax.set_ylabel('Mg/Si', rotation=270, fontsize=fontsize, labelpad=50)
# cbar.ax.tick_params(axis="y", labelsize=ticksize)

cax = inset_axes(
    ax01,
    width="40%",  # width: 50% of parent_bbox width
    height="3%",  # height: 5%
    loc="upper right",
    borderpad=-4,
    # borderpad=0,
    # loc="lower right",
    # bbox_to_anchor=(1, 1.05, 1, 1),
)
cax.xaxis.set_ticks_position("bottom")
cbar = fig.colorbar(im, cax=cax, orientation="horizontal", format="%.1f", ticks=[vmin, 1.0, vmax])
cbar.ax.set_xlabel('Mg/Si', fontsize=ticksize+2, labelpad=-1)
cbar.ax.tick_params(axis="x", labelsize=ticksize-2)

fig.savefig(fig_path + 'ternary_subplots' + str(date) + ftype)
# for some reason this only saves axis labels if individual subplots save=True

# plt.show()

""" test individual """
# _, tax = ternary_scatter(p_of_interest=4, T_of_interest=1373.15, core_eff=88, Xf=7.0, component='Fe2O3', z_var='mgsi',
#                     model='melts', cmap='viridis', vmin=0.69, vmax=1.6, phases=['orthopyroxene', 'clinopyroxene', 'garnet'],
#                    mec='r', marker='o', lw=1.5,
#                 save=True, fig_path=fig_path'#
#  )
