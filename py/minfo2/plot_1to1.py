import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
import meltsfugacitydata as mfug
import perplexfugacitydata as pfug
import numpy as np
import os
import py.main as rw
from py.useful_and_bespoke import colourbar, cornertext
import pandas as pd
from matplotlib import rc
import datetime
from matplotlib.ticker import FormatStrFormatter

today = str(datetime.date.today())
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

core_eff = 88
Xf = 3

# to compare fo2 directly between melts and perplex, use Cr runs
dir1 = fo2plt.output_parent_mlt_earth + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext_Cr/'
dir2 = fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext_Cr/'


def fo2_1to1(dir1, dir2, x_var='logfo2_1GPa', T_iso=1373, z_var=None, cmap=None, c='k', vmin=None, vmax=None,
             xlabel=None, ylabel=None, alpha=0.5, mec=None,
             title=None, s=20, marker='o', model1=None, model2=None, verbose=False, zlabel=None, ticksize=10, dmm=True,
             c_dmm='r', xlims=None,
             labelsize=16, legsize=12, save=True, ffmt='.pdf', fname=None, exclude_names=[], exclude_silica=True,
             x_scale=1, fig_path=fo2plt.figpath, **kwargs):
    # get matching names

    # if cmap and (not vmin or not vmax):
    #     raise NotImplementedError('Cannot colour-code without explicit vmin and vmax')
    # if cmap:
    #     # get normalised cmap
    #     norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    #     colours = cm.ScalarMappable(norm=norm, cmap=cmap)
    if fname is None:
        fname = 'fo2_mdl'
    if xlabel is None:
        xlabel = x_var + ' (' + model1 + ')'
    if ylabel is None:
        ylabel = x_var + ' (' + model2 + ')'
    if zlabel is None:
        zlabel = z_var

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    # ax.set_title(title, fontsize=labelsize)
    ax = cornertext(ax, title, pos='top left', size=labelsize)

    # parse X_ferric
    spl = os.path.dirname(dir1).split('/')[-1].split('_')
    X_ferric = None
    for sp in spl:
        if 'ferric' in sp:
            X_ferric = int(''.join(filter(str.isdigit, sp)))

    # get directory names in folder
    subfolders = rw.get_run_dirs(output_path=dir1)
    df = pd.DataFrame(columns=['name', 'x1', 'x2', 'z'], index=range(len(subfolders)))  # tmp df for data to plot
    model1_nanlist = []
    idx = 0
    if subfolders:  # nonzero
        for ii, sub in enumerate(subfolders):
            name = os.path.basename(sub)
            if (len(os.listdir(sub)) > 0) and (name not in exclude_names):
                # load results data
                if model1 == 'melts':
                    dat = mfug.init_from_results(name, X_ferric=X_ferric, output_parent_path=dir1, T_final=T_iso,
                                                 load_results_csv=True, verbose=False)
                elif model1 == 'perplex':
                    dat = pfug.init_from_results(name, X_ferric=X_ferric, output_parent_path=dir1, T_iso=T_iso,
                                                 load_results_csv=True, verbose=False)

                if dat is not None:  # if results exist
                    if exclude_silica:
                        if 'X_q' in dat.data.columns:  # don't print quartz saturation cases
                            if verbose:
                                print('dropping case with quartz:', dat.name)
                            continue

                    if np.isnan(eval('dat.' + x_var)):
                        # print(name, '.....', x_var, '=', eval('dat.' + x_var), '@', T_iso, 'K', 'try processing again?')
                        model1_nanlist.append(dat.name)
                        df.x1.iloc[ii] = np.nan
                    else:
                        df.x1.iloc[ii] = eval('dat.' + x_var) * x_scale
                    df.name.iloc[ii] = name

                    if z_var is not None:
                        df.z.iloc[ii] = eval('dat.' + z_var)
                    if 'X_q' in dat.data.columns:
                        print('dropping silica failure!!!!!!!!!!!!')
                    # print('added row', ii, df.iloc[ii])
                    idx += 1

    # df.dropna(inplace=True, subset=['x1'])
    # print('df after dir1\n', df.head())

    # get matching runs
    model2_nanlist = []
    for ii in range(len(df)):
        name = df.name.iloc[ii]
        try:
            if np.isnan(name):
                continue  # skip this row
        except TypeError:
            pass  # name is a string

        if model2 == 'melts':
            dat = mfug.init_from_results(name, X_ferric=X_ferric, output_parent_path=dir2, verbose=False,
                                         load_results_csv=True, **kwargs)
        elif model2 == 'perplex':
            dat = pfug.init_from_results(name, X_ferric=X_ferric, output_parent_path=dir2, verbose=False,
                                         load_results_csv=True, **kwargs)

        # if the results exist
        if dat is not None:
            if exclude_silica:
                if 'X_q' in dat.data.columns:  # don't print quartz saturation cases
                    if verbose:
                        print('dropping case with quartz:', dat.name)
                    continue
            try:
                dat.read_fo2_results(verbose=False)
                if np.isnan(eval('dat.' + x_var)):
                    # print(name, '.....', x_var, '=', eval('dat.' + x_var), '@', T_iso, 'K, try processing again?')
                    model2_nanlist.append(dat.name)
                    df.x2.iloc[ii] = np.nan
                else:
                    df.x2.iloc[ii] = eval('dat.' + x_var) * x_scale
            except FileNotFoundError:
                if verbose:
                    print('results file not found:', dir2 + name)
    # df.dropna(inplace=True, subset=['x2'])

    print('df\n', df.head(), '\nn =', len(df.dropna(axis=0, inplace=False)))

    print('\n\n\n', model1, x_var, '\nnanlist (reprocess?)\n', model1_nanlist)
    print('\n\n\n', model2, x_var, '\nnanlist (reprocess?)\n', model2_nanlist)

    # set colours
    if z_var is not None:
        c = df.z
        if vmin is None:
            vmin = np.min(c)
        if vmax is None:
            vmax = np.max(c)
    else:
        c = c
        cmap = None

    sc = ax.scatter(df.x1, df.x2, c=c, edgecolors=mec, cmap=cmap, vmin=vmin, vmax=vmax, s=s, marker=marker, alpha=alpha)

    # make cbar
    if z_var is not None:
        # make dummy scatter with no alpha
        sc = ax.scatter(df.x1, df.x2, c=c, cmap=cmap, vmin=vmin, vmax=vmax, s=0, alpha=1)
        cbar = colourbar(sc, ax=ax, ticksize=ticksize, format="%.1f")
        cbar.ax.set_ylabel('Mg/Si', rotation=270, fontsize=labelsize, labelpad=20)

    # 1:1 line
    if xlims is None:
        xlims = ax.get_xlim()
    delta = [-2, -1, 0, 1, 2]
    lws = [0.3, 0.5, 0.9, 0.5, 0.3]
    for dl, lw in zip(delta, lws):
        ax.plot(xlims, np.array(xlims) + dl, c='k', ls=(1, (5, 10)), lw=lw, alpha=0.9)

    ax.set_xlim(xlims)
    ax.set_ylim(xlims)
    ax.tick_params(axis='both', labelsize=ticksize)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    if dmm:
        # add DMM reference point
        print('adding dmm')
        mdat = mfug.init_from_results('Stolper', output_parent_path=mfug.output_parent_default,
                                      load_results_csv=True, verbose=False)
        pdat = pfug.init_from_results('Stolper', output_parent_path=pfug.output_parent_default,
                                      load_results_csv=True, verbose=False)
        x, y = eval('mdat.' + x_var) * x_scale, eval('pdat.' + x_var) * x_scale

        ax.scatter(x, y, c=eval('mdat.' + z_var), cmap=cmap, vmin=vmin, vmax=vmax, edgecolors=c_dmm, marker='o', s=s)
        ax.text(x+0.07, y-0.07, 'DMM', ha='left', va='top', c=c_dmm, fontsize=labelsize)

    if save:
        fig.savefig(fig_path + fname + ffmt, bbox_inches='tight')
    return fig, ax


markersize = 40
labelsize = 16
legsize = 12
alpha = 0.3
vmin, vmax = 0.69, 1.6
mec = 'xkcd:midnight blue'
cmap = 'viridis'
c_dmm = 'xkcd:blood orange'
xlabel = 'log(fO$_2$), pMELTS'
ylabel = 'log(fO$_2$), Perple_x'
zlabel = 'Mg/Si'
ffmt = '.pdf'

# 1 GPa
fig, ax = fo2_1to1(dir1, dir2, x_var='logfo2_1GPa', z_var='mgsi', cmap=cmap, mec=mec, c_dmm=c_dmm,
                   xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,  x_scale=1, alpha=alpha,
                   title='1 GPa', s=markersize, marker='o', model1='melts', model2='perplex', vmin=vmin, vmax=vmax,
                   labelsize=labelsize, legsize=legsize, save=True, fname='logfo2_1GPa_mdls' + today, ffmt=ffmt,
                   exclude_names=[], exclude_silica=True, xlims=(-13.2, -9), verbose=True
                   )

# 4 GPa
fig, ax = fo2_1to1(dir1, dir2, x_var='logfo2_4GPa', z_var='mgsi', cmap=cmap, mec=mec, c_dmm=c_dmm,
                   xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, x_scale=1, alpha=alpha,
                   title='4 GPa', s=markersize, marker='o', model1='melts', model2='perplex', vmin=vmin, vmax=vmax,
                   labelsize=labelsize, legsize=legsize, save=True, fname='logfo2_4GPa_mdls' + today, ffmt=ffmt,
                   exclude_names=[], exclude_silica=True, xlims=(-11, -6))

# plt.show()


def fo2_1o1_subplot(**kwargs):
    """ vertically stacked pressures, 1 cbar """