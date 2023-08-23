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

source = fo2plt.output_parent_mlt_earth
# source = fo2plt.output_parent_px
dirs = [source + 'hypatia_88coreeff_3ferric_ext' + ext + '/' for ext in ['', '_Cr']]

for ii, diri in enumerate(dirs):
    fo2plt.pop_hist(dirs=[diri], x_var='delta_qfm_1GPa', p_of_interest=1, figpath=fo2plt.figpath,
                    save=True, make_legend=False, T_of_interest=1373, fname=['pop_noCr', 'pop_Cr'][ii],
                    model='melts', show_sd='corner', c=['r', 'b'][ii],
                    xlabel='QFM', z_labels=['no Cr', 'Cr'][ii],
                    nbins=20, alpha=0.5, z_val=1,  # dummy
                    hist_kwargs={})
