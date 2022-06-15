import pandas as pd
import numpy as np
import perplexdata as px
import main as rw
import parameters as p
import ask_hypatia as hyp

perplex_path = px.perplex_path_default
directory = 'earthsize_planets_1900K_88Fe/'
df = pd.read_csv('/home/claire/Works/rocky-water/py/earthsize_planets.csv')
print(df.head())

host_names = hyp.random_star(n=-1)

core_eff = 0.8831461545794602
Tp = 1900

for index, row in df.iterrows():
    star = row['star_name']
    if (star == star) and star in host_names:  # is not NaN and in host_names
        print('running star', star)
        dat = rw.build_planet(Tp=Tp, M_p=row['M_p'] * p.M_E, star=star,
                              name=row['planet_name'].replace(' ', ''),
                              get_saturation=True,
                              output_parent_path=perplex_path + 'output/' + directory,
                              use_local_composition=True
                              )