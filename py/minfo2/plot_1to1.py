import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt
import meltsfugacitydata as mfug
import perplexfugacitydata as pfug

core_eff = 88
Xf = 3
alpha=0.9

dir1 = fo2plt.output_parent_mlt + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext_Cr/'
dir2 = fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext_Cr/'

fig, ax = fo2plt.fo2_1to1(dir1, dir2, x_var='logfo2_4GPa', z_var='mgsi', cmap='viridis', c='k', vmin=None, vmax=None,
                          xlabel='log(fO$_2$), MELTS', ylabel='log(fO$_2$), Perple_x', zlabel='Mg/Si',
                          title='4 GPa', s=40, marker='o', model1='melts', model2='perplex',
                          labelsize=16, legsize=12, save=True, fname='logfo2_4GPa_mdls', ffmt='.png', exclude_names=[],
                          exclude_silica=True, x_scale=1, alpha=alpha)

fig, ax = fo2plt.fo2_1to1(dir1, dir2, x_var='logfo2_1GPa', z_var='mgsi', cmap='viridis', c='k', vmin=None, vmax=None,
                          xlabel='log(fO$_2$), MELTS', ylabel='log(fO$_2$), Perple_x', zlabel='Mg/Si',
                          title='1 GPa', s=40, marker='o', model1='melts', model2='perplex',
                          labelsize=16, legsize=12, save=True, fname='logfo2_1GPa_mdls', ffmt='.png', exclude_names=[],
                          exclude_silica=True, x_scale=1, xlims=(-13.5, -8.5), alpha=alpha)

plt.show()
