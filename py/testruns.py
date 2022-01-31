import numpy as np
import perple_x as px
import parameters as p
import eos
import os

# stars = ["HIP 98355", "sun"]
# px.build_planet(name='marstest',
#           M_p=0.1*p.M_E, stars='sun', core_efficiency=0.8,  # planet parameters
#           n=1000, tol=0.1,  # interior structure numerical parameters
#           overwrite=True, verbose=False  # run settings
#           )


# # test different pressure profiles based on mass - what is central value?
px.build_planet(name='test', M_p=p.M_E, test_CMF=0.33, test_oxides=px.wt_oxides_Earth,
                plot=True, store_all=True, get_saturation=True, overwrite=True,
                maxIter=50)


# # test Hakim EoS
# # pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, psurf, tsurf)
# _, density, alpha, cp = eos.EOS_all(234.399, 1700, 4)
# print('rho Bouchet', density)
# _, density, alpha, cp = eos.EOS_all(234.4, 1700, 4)
# print('rho Hakim', density)