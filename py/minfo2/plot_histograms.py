import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
import numpy as np
import os
import py.main as rw
import perplexfugacitydata as pfug
import meltsfugacitydata as mfug
from py.useful_and_bespoke import dark_background, colorize, colourbar, iterable_not_string
import cmcrameri
import datetime
import matplotlib  # for running remotely
from matplotlib import rc
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

""" compare fo2 hist for multiple runs """

today = datetime.date.today()
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
fformat = '.pdf'


def pop_hist_subplots(fig=None, axes=None, dirs=None, z_var=None, z_list=None, labelsize=16,
                      ticksize=12, labelpad=None, xlabel=None, p_of_interest=None, vmin=None, vmax=None,
                      xlim=None, ylim=None, labelpressure=True, **kwargs):
    print('vmin vmax', vmin, vmax)
    for zz, ax in enumerate(axes):
        # print('zz', zz, 'dirs[zz]', dirs[zz])
        z_val = z_list[zz]
        fig, ax = fo2plt.pop_hist([dirs[zz]], z_var=z_var, z_val=z_val, save=False, fig=fig, ax=ax, vmin=vmin,
                                  vmax=vmax,
                                  xlabel='', title=None, labelsize=labelsize, ticksize=ticksize,
                                  labelpad=labelpad, **kwargs)

        # set axis limits and ticks
        ax.set_yticks([])
        if zz < len(axes) - 1:
            ax.set_xticks([])
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ymid = np.ptp(np.array(ax.get_ylim())) / 2
        if z_var == 'core_eff':
            ax.set_yticks(ticks=[ymid], labels=[(100 - z_val) / 100])  # label y axis midpoints with z_var
        else:
            ax.set_yticks(ticks=[ymid], labels=[z_val / 100])  # label by z_var
        ax.tick_params(axis='y', width=0)  # hide tick lines on y axis
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # label pressure
    if labelpressure:
        axes[0].text(0.015, 0.8, str(int(p_of_interest)) + ' GPa', transform=axes[0].transAxes, size=labelsize,
                 va='top', ha='left')

    # make xlabel
    if xlabel:
        axes[-1].set_xlabel(xlabel,  # at ' + str(p_of_interest) + ' GPa',
                            fontsize=labelsize, labelpad=labelpad)
    return fig, axes


def pop_hist_allpressures_allmodels(z_var=None, ylabel=None, dirs=None, save=True,
                                    labelsize=16, ylabelpad=None, **kwargs):
    # set up gridspec for 2 pressures
    fig = plt.figure(figsize=(12, 5))  # too small?
    gs = gridspec.GridSpec(2, 2, height_ratios=(1, 1), width_ratios=(1, 1),
                           hspace=0.15, wspace=0.15, figure=fig)

    for col, model in enumerate(['melts', 'perplex']):

        # set values
        if model == 'melts':
            source = fo2plt.output_parent_mlt_earth
            model_str = 'pMELTS'
            if z_var == 'core_eff':
                xlim = -3, 1.5
                ylim = 0, 2
            elif z_var == 'X_ferric':
                xlim = -4.2, 2.4
                ylim = 0, 2.5
        elif model == 'perplex':
            source = fo2plt.output_parent_px
            model_str = r'Perple\_X'
            if z_var == 'core_eff':
                xlim = -5, 0.5
                ylim = 0, 1.5
            elif z_var == 'X_ferric':
                xlim = -6.4, 2
                ylim = 0, 1.5

        if z_var == 'core_eff':
            z_list = core_eff_list
            dirs = [source + 'hypatia_' + str(ce) + 'coreeff_' + str(Xf_constant) + 'ferric_ext/' for ce in z_list]
            ylabel = r'Fe$_{\rm mantle}$/(Fe$_{\rm core}$ + Fe$_{\rm mantle}$)'
            cmap = 'autumn'
            vmin, vmax = min(z_list), max(z_list) + 13
            z_scale = 100

        elif z_var == 'X_ferric':
            z_list = Xf_list
            dirs = [source + 'hypatia_' + str(core_eff_constant) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
                    z_list]
            ylabel = r'Fe$^{{\rm 3+}}/\Sigma$Fe'
            cmap = 'autumn_r'
            vmin, vmax = min(z_list) - 3.5, max(z_list)
            z_scale = 1

        nrows = len(z_list)
        gs0 = gridspec.GridSpecFromSubplotSpec(nrows, 1, subplot_spec=gs[0, col], height_ratios=[1] * nrows, hspace=0)
        gs1 = gridspec.GridSpecFromSubplotSpec(nrows, 1, subplot_spec=gs[1, col], height_ratios=[1] * nrows, hspace=0)
        axes0 = []
        [axes0.append(fig.add_subplot(gs0[row])) for row in range(nrows)]
        axes1 = []
        [axes1.append(fig.add_subplot(gs1[row])) for row in range(nrows)]
        # fig, axes = plt.subplots(len(dirs), 1, figsize=(6, 4))

        # 1 GPa
        fig, axes0 = pop_hist_subplots(dirs=dirs, x_var='delta_qfm_1GPa',
                                       p_of_interest=1, axes=axes0,
                                       model=model,
                                       z_var=z_var, z_list=z_list,
                                       fig=fig,
                                       labelsize=labelsize, verbose=False,
                                       z_scale=z_scale, ylabel=ylabel,
                                       xlim=xlim, ylim=ylim,
                                       cmap=cmap, vmin=vmin, vmax=vmax,
                                       **kwargs)
        # 4 GPa
        fig, axes1 = pop_hist_subplots(dirs=dirs, x_var='delta_qfm_4GPa',
                                       p_of_interest=4, axes=axes1, xlabel='$\Delta$FMQ',
                                       model=model,
                                       z_var=z_var, z_list=z_list,
                                       labelsize=labelsize, verbose=False,
                                       fig=fig,
                                       xlim=xlim, ylim=ylim,
                                       cmap=cmap, vmin=vmin, vmax=vmax,
                                       **kwargs)

        axes0[-1].tick_params(labelbottom=False)
        axes0[0].set_title(r'$f_{{\rm O}_2}$ distribution with ' + model_str, fontsize=labelsize)

        for axs in (axes0, axes1):
            for a in axs:
                a.xaxis.set_minor_locator(AutoMinorLocator())

    print('ylabelpad', ylabelpad)
    fig.supylabel(ylabel, fontsize=labelsize, x=ylabelpad)

    # plt.tight_layout()
    # plt.subplots_adjust(hspace=0)
    if save:
        fig.savefig(fo2plt.figpath + 'hist_' + z_var + str(today) + fformat, bbox_inches='tight',
                    facecolor=fig.get_facecolor(),
                    dpi=400)
    return fig


""" make plots below """




# set these
Xf_constant = 3  # only used in dirs
core_eff_constant = 88  # only used in dirs
Xf_list = [9, 7, 5, 3, 1]
core_eff_list = [70, 80, 88, 95, 99]
T_of_interest = 1473

labelsize = 18
ticksize = 12
legsize = 12
labelpad = 10
ylabelpad = 0.06  # not used in supylabel
nbins = 35
# z_var = 'core_eff'
z_var = 'X_ferric'
# publication version immediately below
# pop_hist_allpressures_allmodels(z_var=z_var,
#                                 labelsize=labelsize, ticksize=ticksize, legsize=legsize, labelpad=labelpad,
#                                 ylabelpad=ylabelpad,
#                                 nbins=nbins, lw=2, T_of_interest=T_of_interest,
#                                 save=True,
#                                 exclude_silica=True)

#
# def pop_hist_allpressures(z_var=None, z_list=None, model=None, model_str=None, ylabel=None, dirs=None, save=True,
#                           labelsize=16, ylabelpad=None, **kwargs):
#     # set up gridspec for 2 pressures
#     nrows = len(z_list)
#     fig = plt.figure(figsize=(6, 5))  # too small?
#     gs = gridspec.GridSpec(2, 1, height_ratios=(1, 1), hspace=0.15, figure=fig)
#     gs0 = gridspec.GridSpecFromSubplotSpec(nrows, 1, subplot_spec=gs[0], height_ratios=[1] * nrows, hspace=0)
#     gs1 = gridspec.GridSpecFromSubplotSpec(nrows, 1, subplot_spec=gs[1], height_ratios=[1] * nrows, hspace=0)
#     axes0 = []
#     [axes0.append(fig.add_subplot(gs0[row])) for row in range(nrows)]
#     axes1 = []
#     [axes1.append(fig.add_subplot(gs1[row])) for row in range(nrows)]
#     # fig, axes = plt.subplots(len(dirs), 1, figsize=(6, 4))
#
#     # 1 GPa
#     fig, axes0 = pop_hist_subplots(dirs=dirs, x_var='delta_qfm_1GPa',
#                                    p_of_interest=1, axes=axes0,
#                                    model=model,
#                                    z_var=z_var, z_list=z_list,
#                                    fig=fig,
#                                    labelsize=labelsize, verbose=False, **kwargs)
#     # 4 GPa
#     fig, axes1 = pop_hist_subplots(dirs=dirs, x_var='delta_qfm_4GPa',
#                                    p_of_interest=4, axes=axes1, xlabel='$\Delta$FMQ',
#                                    model=model,
#                                    z_var=z_var, z_list=z_list,
#                                    labelsize=labelsize, verbose=False,
#                                    fig=fig, **kwargs)
#
#     axes0[-1].tick_params(labelbottom=False)
#     fig.supylabel(ylabel, fontsize=labelsize, x=ylabelpad)
#     fig.suptitle(r'$f_{{\rm O}_2}$ distribution with ' + model_str, fontsize=labelsize, y=0.93)
#
#     # plt.tight_layout()
#     # plt.subplots_adjust(hspace=0)
#     if save:
#         fig.savefig(fo2plt.figpath + 'hist_' + z_var + '_' + model + str(today) + fformat, bbox_inches='tight',
#                     facecolor=fig.get_facecolor(),
#                     dpi=400)
#     return fig

# ylabelpad = 10  # not used in supylabel
# model = 'melts'
# # model = 'perplex'
# z_var = 'core_eff'
# # z_var = 'X_ferric'
#
# if model == 'melts':
#     source = fo2plt.output_parent_mlt_earth
#     model_str = 'pMELTS'
#     if z_var == 'core_eff':
#         xlim = -3, 1.5
#         ylim = 0, 2
#     elif z_var == 'X_ferric':
#         xlim = -4.2, 2.4
#         ylim = 0, 2.5
# elif model == 'perplex':
#     source = fo2plt.output_parent_px
#     model_str = r'Perple\_X'
#     if z_var == 'core_eff':
#         xlim = -5, 0.5
#         ylim = 0, 1.5
#     elif z_var == 'X_ferric':
#         xlim = -6.4, 2
#         ylim = 0, 1.5
#
# if z_var == 'core_eff':
#     z_list = core_eff_list
#     dirs = [source + 'hypatia_' + str(ce) + 'coreeff_' + str(Xf_constant) + 'ferric_ext/' for ce in z_list]
#     ylabel = r'Fe$_{\rm mantle}$/(Fe$_{\rm core}$ + Fe$_{\rm mantle}$)'
#     ylabelpad = 0.0
#     cmap = 'autumn'
#     vmin, vmax = min(z_list), max(z_list) + 5
#     z_scale = 100
# elif z_var == 'X_ferric':
#     z_list = Xf_list
#     dirs = [source + 'hypatia_' + str(core_eff_constant) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
#             z_list]
#     ylabel = r'Fe$^{{\rm 3+}}/\Sigma$Fe'
#     ylabelpad = 0  # 0.04
#     cmap = 'autumn_r'
#     vmin, vmax = min(z_list) - 1.5, max(z_list)
#     z_scale = 1

# xlim, ylim=(1e-7, 1e-6), None
# pop_hist_allpressures(z_var=z_var, z_list=z_list, z_scale=z_scale, model=model, ylabel=ylabel, dirs=dirs,
#                       xlim=xlim, ylim=ylim, model_str=model_str,
#                       labelsize=labelsize, ticksize=ticksize, legsize=legsize, labelpad=labelpad, ylabelpad=ylabelpad,
#                       nbins=nbins, cmap=cmap, vmin=vmin, vmax=vmax,  lw=2, T_of_interest=T_of_interest, save=True,
#                       exclude_silica=True)








import matplotlib.colors as col
import seaborn as sns
flat_huslmap = col.ListedColormap(sns.color_palette('husl',256))

cm = 1/2.54  # centimeters in inches
rc('font',**{'family':'serif','serif':['Times']})
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
z_list = Xf_list
fig = plt.figure(figsize=(14*cm, 12*cm))
gs = gridspec.GridSpec(1, 1,
                       hspace=0.15, wspace=0.15, figure=fig)
nrows = len(z_list)
gs0 = gridspec.GridSpecFromSubplotSpec(nrows, 1, subplot_spec=gs[0], height_ratios=[1] * nrows, hspace=0)
axes0 = []
[axes0.append(fig.add_subplot(gs0[row])) for row in range(nrows)]

model='melts'
xlabel = r'Upper mantle log$f_{\rm O_2}$ ($\Delta$FMQ)'
source = fo2plt.output_parent_mlt_earth
dirs = [source + 'hypatia_' + str(core_eff_constant) + 'coreeff_' + str(Xf) + 'ferric_ext/' for Xf in
        z_list]
xlim = -3.2, 2.4
ylim = 0, 1.5
ylabel = r'Fe$^{{\rm 3+}}/\Sigma$Fe'
cmap = 'autumn_r'
vmin, vmax = min(z_list) - 3.5, max(z_list)
z_scale = 1
labelsize = 16
tick_width = 2
major_length = 4
minor_length = 2
# cmap = flat_huslmap

fig, axes = pop_hist_subplots(dirs=dirs, x_var='delta_qfm_4GPa',
                               p_of_interest=4, axes=axes0, xlabel=xlabel,
                               model=model,  labelpressure=False,
                               z_var=z_var, z_list=z_list,
                               labelsize=labelsize, ticksize=16, legsize=16, verbose=False,
                               fig=fig, alpha=0.9,
                               xlim=xlim, ylim=ylim, show_sd='leftcorner',
                               cmap=cmap, vmin=vmin, vmax=vmax,)

axes[-1].xaxis.set_minor_locator(AutoMinorLocator())
axes[-1].tick_params(axis='x', which='major', length=major_length, width=tick_width)
axes[-1].tick_params(axis='x', which='minor', length=minor_length, width=tick_width)
axes[0].set_title(r'$f_{\rm O_2}$ modulation by bulk silicate composition:' + '\npMELTS calculations', fontsize=labelsize)
axes[2].set_ylabel(ylabel, fontsize=labelsize, labelpad=12)

for ax in axes:
    [spine.set_linewidth(tick_width) for spine in ax.spines.values()]


# fig, ax = dark_background(fig, ax)
fig.savefig('/home/claire/Works/min-fo2/shirt/images/histogram2.png', bbox_inches='tight', transparent=True,
            dpi=300)
# plt.show()
