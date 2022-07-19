import pandas as pd
import numpy as np
import perplexdata as px
import main as rw
import parameters as p
import ask_hypatia as hyp

perplex_path = px.perplex_path_default
<<<<<<< HEAD
df = pd.read_csv('/home/claire/Works/rocky-water/py/earthsize_planets_FGKM.csv')
print(df.head(20))
host_names = hyp.random_star(n=-1)
core_eff = 0.8831461545794602
Tp = 1600

directory = 'earthsize_planets_1600K_88Fe/'
for index, row in df.iterrows():
    star = row['star_name']
    if (star == star) and star in host_names:  # is not NaN and in host_names
        print('running star', star)
        dat = rw.build_planet(Tp=Tp, M_p=row['M_p'] * p.M_E, star=star,
                              name=row['planet_name'].replace(' ', ''),
                              get_saturation=True, n=800, core_efficiency=core_eff,
                              output_parent_path=perplex_path + 'output/' + directory,
                              use_local_composition=False,
                              )

directory = 'earthsize_planets_1600K_88Fe_mgsimin/'
for index, row in df.iterrows():
    star = row['star_name']
    if (star == star) and star in host_names:  # is not NaN and in host_names
        print('running star', star)
        dat = rw.build_planet(Tp=Tp, M_p=row['M_p'] * p.M_E, star=star,
                              name=row['planet_name'].replace(' ', ''),
                              get_saturation=True, n=800, core_efficiency=core_eff,
                              output_parent_path=perplex_path + 'output/' + directory,
                              use_local_composition=False, get_hypatia_min=['mg'], get_hypatia_max=['si'],
                              )

directory = 'earthsize_planets_1600K_88Fe_mgsimax/'
=======
directory = 'earthsize_planets_1900K_88Fe/'
df = pd.read_csv('/home/claire/Works/rocky-water/py/earthsize_planets.csv')
print(df.head())

host_names = hyp.random_star(n=-1)

core_eff = 0.8831461545794602
Tp = 1900

>>>>>>> 780fc79a30ed92265a34ee2e0dc76d862124a3a0
for index, row in df.iterrows():
    star = row['star_name']
    if (star == star) and star in host_names:  # is not NaN and in host_names
        print('running star', star)
        dat = rw.build_planet(Tp=Tp, M_p=row['M_p'] * p.M_E, star=star,
                              name=row['planet_name'].replace(' ', ''),
<<<<<<< HEAD
                              get_saturation=True, n=800, core_efficiency=core_eff,
                              output_parent_path=perplex_path + 'output/' + directory,
                              use_local_composition=False, get_hypatia_min=['si'], get_hypatia_max=['mg'],
=======
                              get_saturation=True,
                              output_parent_path=perplex_path + 'output/' + directory,
                              use_local_composition=True
>>>>>>> 780fc79a30ed92265a34ee2e0dc76d862124a3a0
                              )