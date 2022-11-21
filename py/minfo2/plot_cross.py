import matplotlib.pyplot as plt
import oxygen_fugacity_plots as fo2plt

""" make cross plot of oxides / elements"""

core_eff = 88
Xf = 3
p_of_interest = 1
exclude_silica = True
model = 'melts'

if model == 'melts':
    source = fo2plt.output_parent_mlt
elif model == 'perplex':
    source = fo2plt.output_parent_px

if p_of_interest == 1:
    ylim = (-14, -9)
    phcomps = ['Ol', 'Opx', 'Cpx', 'Sp']

elif p_of_interest == 4:
    ylim = (-11.5, -7)
    phcomps = ['Ol', 'Opx', 'Cpx', 'Gt']

if not exclude_silica:
    phcomps.append('q')

output_parent_path = source + 'hypatia_' + str(core_eff) + 'coreeff_' + str(Xf) + 'ferric_ext_Cr/'

fo2plt.element_xplot(p_of_interest=p_of_interest, components=['MgO', 'SiO2', 'Al2O3', 'FeO', 'CaO'],
                     output_parent_path=output_parent_path,
                     ylim=ylim, linec='k', labelsize=16, save=True, fname='mlt_crossplot_oxides_' + str(p_of_interest),
                     model=model,
                     exclude_silica=exclude_silica)

fo2plt.element_xplot(p_of_interest=p_of_interest, components=phcomps, output_parent_path=output_parent_path,
                     ylim=ylim, linec='k', labelsize=16, save=True, fname='mlt_crossplot_phases_' + str(p_of_interest),
                     model=model,
                     exclude_silica=exclude_silica)

plt.show()
