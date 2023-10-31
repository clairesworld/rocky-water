import matplotlib.pyplot as plt
import numpy as np
import os
import pickle as pkl
import py.main as rw
from py.useful_and_bespoke import dark_background, colorize, colourbar, iterable_not_string
import datetime
import matplotlib  # for running remotely
from matplotlib import rc
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

perplex_output_path = '/home/g/guimond/Work/hypatia_compositions/'

def plot_phase_hist2d(x_name='mgsi', output_path=perplex_output_path, range_min=None, range_max=None, x_scale=1,
                      labelsize=16, legsize=10, pressure=5, c_hist='0.1', alpha_verts=0.1, max_mode=0, test=None,
                      save=True, fig_path=''):

    directory = 'hypatia1M_1600K_88Fe_hires'

    # set up scatter & histogram gridspec
    dats = rw.read_dir(output_path + directory + '/', verbose=True)
    if test is not None:
        dats = dats[:test]

    fig, ax = plt.subplots(figsize=(6,2))

    # find phases

    # all have gt and cpx?
    c_opx_qz = '#FEDB85'  #'xkcd:orangey yellow'
    c_ol_opx = '#E5FE85'  #'#CDFD15'
    c_ol = '#FE9E85'  ##FD4515'

    # plot histogram
    x = []
    for dat in dats:
        try:
            xi = eval('dat.' + x_name) * x_scale

            if (range_min is None) or (range_min is not None and xi >= range_min):
                if (range_max is None) or (range_max is not None and xi <= range_max):
                    x.append(xi)
        except KeyError:
            print(dat.name, 'does not have attribute', x_name)
            continue

        # get phase diagram
        idx = (np.abs(dat.df_comp['P(bar)'] - (pressure*1e4))).argmin()  # closest pressure
        row = dat.df_comp.iloc[idx]
        # print(dat.name, row)
        # print(('O' in row.index))

        has_ol = ('O' in row.index) and (row['O'] > max_mode)
        has_opx = ('Opx' in row.index) and (row['Opx'] > max_mode)
        has_sio2 = (('qtz' in row.index) and (row['qtz'] > max_mode)) or (('coes' in row.index) and (row['coes'] > max_mode))

        alpha_verts=1
        lw_verts=0.7
        if has_sio2 and has_opx and (not has_ol):
            ax.axvline(xi, c=c_opx_qz, alpha=alpha_verts, zorder=7, lw=lw_verts)
        elif has_opx and has_ol and (not has_sio2):
            ax.axvline(xi,c=c_ol_opx, alpha=alpha_verts, zorder=5, lw=lw_verts)
        elif has_ol and (not has_opx):
            ax.axvline(xi,  c=c_ol, alpha=alpha_verts, zorder=6, lw=lw_verts)
        else:
            print(dat.name, 'no composition found')
            print(row)

    # plot hist
    ax.hist(x, color=c_hist, edgecolor=c_hist, bins=50, density=True, stacked=True, zorder=10)

    # Hide the right and top spines
    ax.spines[['right', 'top', 'left']].set_visible(False)
    ax.get_yaxis().set_visible(False)

    # labels
    y_frac_labels = 0.9
    ax.text(0.25, y_frac_labels, 'opx+coes\n+cpx+gt', transform=ax.transAxes, va='top', ha='center', c='k', fontsize=legsize, zorder=11)
    ax.text(0.53, y_frac_labels, 'ol+opx\n+cpx+gt', transform=ax.transAxes, va='top', ha='center', c='k', fontsize=legsize, zorder=11)
    ax.text(0.85, y_frac_labels, 'ol\n+cpx+gt', transform=ax.transAxes, va='top', ha='center', c='k', fontsize=legsize, zorder=11)
    ax.text(0.05, 0.05, '{:d} GPa'.format(pressure), transform=ax.transAxes, va='bottom', ha='left', fontsize=labelsize, zorder=11)

    ax.set_xlim(0, 2)
    ax.set_ylim(0, 5.5)
    ax.set_xlabel('Mg/Si', fontsize=labelsize)

    print('ylim', ax.get_ylim())

    if not save:
        plt.show()
    else:
        fig.savefig(fig_path + 'mgsi_phases.pdf', bbox_inches='tight')
        return fig, ax


fig, ax = plot_phase_hist2d(x_name='mgsi', output_path=perplex_output_path, test=None, save=True,
                            fig_path='/home/g/guimond/Work/rocky-water/figs_scratch/')
