import numpy as np
import perple_x as run
import parameters as p
import eos

# stars = ["HIP 98355", "sun"]
# run.start(name='marstest',
#           M_p=0.1*p.M_E, stars='sun', core_efficiency=0.8,  # planet parameters
#           n=1000, tol=0.1,  # interior structure numerical parameters
#           overwrite=True, verbose=False  # run settings
#           )


# # test different pressure profiles based on mass - what is central value?
M = [1, 5]
for m in M:
    dat = run.PerplexData(name='M' + str(m), M_p=m*p.M_E, star='sun', core_efficiency=None)
    dat.CMF = 0.33
    dat.wt_oxides = run.wt_oxides_Mars
    dat.iterate_structure(Psurf=1, Tsurf=1300, n=500, maxIter=100, tol=100, clean=True, overwrite=True)
    print('max pressure:', np.max(dat.pressure)*1e-12, 'TPa')
    dat.plot_structure(save=False)


# # test Hakim EoS
# # pressure, temperature = eos.pt_profile(n, radius, density, gravity, alpha, cp, psurf, tsurf)
# _, density, alpha, cp = eos.EOS_all(234.399, 1700, 4)
# print('rho Bouchet', density)
# _, density, alpha, cp = eos.EOS_all(234.4, 1700, 4)
# print('rho Hakim', density)