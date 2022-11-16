import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
import meltsfugacitydata as mfug
import perplexfugacitydata as pfug

core_eff = 88
Xf = 3

dir1 = fo2plt.output_parent_mlt + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/'
dir2 = fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/'

fig, ax = fo2plt.fo2_1to1(dir1, dir2, x_var='logfo2_4GPa', z_var='mgsi', cmap=None, c='k', vmin=None, vmax=None,
                xlabel='log(fO$_2$), MELTS', ylabel='log(fO$_2$), Perple_x', zlabel='Mg/Si',
                title=None, s=20, marker='o', model1='melts', model2='perplex',
                labelsize=16, legsize=12, save=True, fname='1to1', exclude_names=[], exclude_silica=True, x_scale=1)

print('at dmm')

# add stolper
mdat = mfug.init_from_results('Stolper', X_ferric=Xf, output_parent_path=mfug.output_parent_default, load_results_csv=True, verbose=False)
pdat = pfug.init_from_results('Stolper', X_ferric=Xf, output_parent_path=pfug.output_parent_default, load_results_csv=True, verbose=False)
ax.scatter(mdat.logfo2_4GPa, pdat.logfo2_4GPa, c='r', marker='o', s=20)

plt.show()
