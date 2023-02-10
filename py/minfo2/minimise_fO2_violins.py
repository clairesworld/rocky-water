import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator

""" Sean Jordan's implementation of FastChem to calculate fO2-dependent outgassing speciation """

fig, ax = plt.subplots(figsize=(12,10))
gs1 = gridspec.GridSpec(2,2)
ax1 = plt.subplot(gs1[0])
ax2 = plt.subplot(gs1[1])
ax3 = plt.subplot(gs1[2])
ax4 = plt.subplot(gs1[3])

for axis in [ax1,ax2,ax3,ax4]:

    if axis==ax1: 
        chem_in = './fO2_out/claire_project_violins/table_FMQ_solar_COSH_chem.dat'
        norm_in = './fO2_out/claire_project_violins/table_FMQ_solar_COSH_norm.dat'
    if axis==ax2: 
        chem_in = './fO2_out/claire_project_violins/table_FMQ_Venus_atm_COSH_chem.dat'
        norm_in = './fO2_out/claire_project_violins/table_FMQ_Venus_atm_COSH_norm.dat'
    if axis==ax3:
        chem_in = './fO2_out/claire_project_violins/table_IW-2_solar_COSH_chem.dat'
        norm_in = './fO2_out/claire_project_violins/table_IW-2_solar_COSH_norm.dat'
    if axis==ax4:
        chem_in = './fO2_out/claire_project_violins/table_IW-2_Venus_atm_COSH_chem.dat'
        norm_in = './fO2_out/claire_project_violins/table_IW-2_Venus_atm_COSH_norm.dat'

    df = pd.read_csv(chem_in, delimiter='\t')
    df_norm = pd.read_csv(norm_in, delimiter='\t')
    D=[]
    species_list = ['CO','CO2','CH4','H2O','H2S','SO2']
    colours      = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:pink','tab:brown','tab:olive','tab:cyan','tab:grey']
    p=0
    pos = [3,6,9,12,15,18]#[4,8,12,16,20,24]#[3,6,9,12,15,18]#[2,4,6,8,10]
    for species in species_list:
    
        norm = df_norm[species].values
        print(species, norm)
        observable = 1e-6/norm
        logobservable = np.log10(observable)
        print(observable)
        data = df[species].values/norm
        logdata = np.log10(df[species].values/norm)
        random_scatter=np.ones_like(logdata)
        for i in range(len(logdata)):
            random_scatter[i] = np.random.rand()
        #axis.scatter(random_scatter+pos[p]-2.5, logdata, s=1, marker='x', color='grey', alpha=0.5)
        xline = np.linspace(pos[p]-1,pos[p]+1,2)
        #axis.plot(xline,np.ones_like(xline)*logobservable, c='r')
        D.append(logdata)
        p+=1

    vp = axis.violinplot(D, positions = pos, widths=2, showmeans=False, showextrema=False, showmedians=False)#, points=1000)

    col=0
    for body in vp['bodies']:
        body.set_alpha(0.6)
        body.set_facecolor(colours[col])
        body.set_edgecolor(colours[col])
        col+=1
    
    pad_size = 12
    label_size = 24
    major_length = 4
    minor_length = 2
    tick_width = 2
    axes_label = 18
    leg_font = 12
    spec_font = 20
    
    [spine.set_linewidth(tick_width) for spine in axis.spines.values()]
    axis.tick_params(axis='both',which='both',pad=pad_size,labelsize=0,direction='in',top=True, right=True)
    axis.tick_params(axis='both',which='major',length=major_length,width=tick_width)
    axis.tick_params(axis='both',which='minor',length=minor_length,width=tick_width)
    axis.xaxis.set_minor_locator(AutoMinorLocator())
    axis.yaxis.set_minor_locator(AutoMinorLocator())
    axis.set_xticks(pos) 

    axis.set_ylim([-4,4])

    #axis.set_yscale('log')
    #axis.set_ylim([1e-5,1e4])

ax1.tick_params(axis='y',which='both',pad=pad_size,labelsize=label_size,direction='in',top=True, right=True)

ax3.tick_params(axis='both',which='both',pad=pad_size,labelsize=label_size,direction='in',top=True, right=True)
ax3.set_xticklabels([x.replace('2','$_2$').replace('4','$_4$') for x in species_list])

ax4.tick_params(axis='x',which='both',pad=pad_size,labelsize=label_size,direction='in',top=True, right=True)
ax4.set_xticklabels([x.replace('2','$_2$').replace('4','$_4$') for x in species_list])

ax1.set_ylabel('log($\dfrac{X_i}{X_{i, QFM}}$)',fontsize=label_size)
ax3.set_ylabel('log($\dfrac{X_i}{X_{i, IW-2}}$)',fontsize=label_size)

gs1.tight_layout(fig)
gs1.update(wspace=0.1, hspace=0.1)

plt.show()










