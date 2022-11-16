import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt

""" plot isothermal cross section of fo2 """

model = 'melts'
plot_kwargs = {}
core_eff = 88
Xf = 3
exclude_silica = True
exclude_names = [
    # '1M_88Ceff_HIP58237_999K'
    # '1M_88Ceff_2MASS16494226-1932340_999K'
]

if model == 'melts':
    output_parent_path = fo2plt.output_parent_mlt + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/'
    fname ='mlts_fo2_variation_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric'
elif model == 'perplex':
    output_parent_path = fo2plt.output_parent_px + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext/'
    fname ='px_fo2_variation_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric'

# # cmap_var, vmin, vmax = 'mg_number', 86, 94
cmap_var, vmin, vmax, cbar_label = 'mgsi', 0.7, 1.3, 'Mg/Si'
fo2plt.multicomp_xsection(output_parent_path=output_parent_path, cmap='spring', cmap_var=cmap_var, vmin=vmin, vmax=vmax,
                          has_cbar=True, cbar_label=cbar_label,
                          save=True, alpha=0.9, lw=0.5, p_min=1.1, ymin=-13, ymax=-6,
                          hist_y=True, exclude_names=exclude_names,
                          exclude_silica=exclude_silica,
                          model=model, fname=fname, verbose=True, **plot_kwargs)

plt.show()
