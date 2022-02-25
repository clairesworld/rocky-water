import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, splrep, splev
import pandas as pd

df = pd.read_csv('/home/claire/Works/rocky-water/peridotite_solidus.csv', names=['P_GPa', 'T_K'])
yi = df['T_K'].to_numpy()
xi = df['P_GPa'].to_numpy()

# positions to inter/extrapolate
x = np.linspace(5, 200, 50)
# spline order: 1 linear, 2 quadratic, 3 cubic ...
order = 1
# do inter/extrapolation
s = InterpolatedUnivariateSpline(xi, yi, k=order)
y = s(x)

# example showing the interpolation for linear, quadratic and cubic interpolation
plt.figure()
plt.plot(xi, yi)
for order in range(1, 4):
    s = InterpolatedUnivariateSpline(xi, yi, k=order)
    y = s(x)
    plt.plot(x, y, label=str(order))
plt.scatter(xi, yi, c='k')

tck = splrep(xi, yi, k=1)
print(tck)
y2 = splev(x, tck)
plt.plot(x, y2, 'k')
plt.legend()





plt.show()

