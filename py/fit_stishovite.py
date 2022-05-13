import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

""" use data from Panero+03 as covers T range 
DISCONTINUED - using Panero & Stixrude 2003 theoretical value instead
"""


def func(x, a, b):
    return a + b / x

T = [1597.7528089887642, 1958.4269662921347, 2012.3595505617977, 2153.9325842696626]  #, 2767.4157303370785]
c_h2o = [74.4444444444444, 80, 207.77777777777771, 127.77777777777771]  #, 84.44444444444434]

# other data
# Bolfan-Casanova+ 2000 - note Al-free
# T.append(1500+273)
# c_h2o.append(72)

# Liu+21
T.append(1900)
c_h2o.append(295)


popt, pcov = curve_fit(func, T, c_h2o)

plt.plot(T, c_h2o, 'k.')
plt.plot(T, func(T, *popt), 'r-',)
plt.show()


