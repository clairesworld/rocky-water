import perplexdata as px
import numpy as np
from scipy import interpolate
import pandas as pd


def read_qfm_os(T0, p0, path=px.perplex_path_default, fin='data_tables/fmqNl_fo2_oli.dat', verbose=False):
    """ read in oli's qfm calculations, T in K, p in bar """
    df = pd.read_csv(path + fin, delimiter=r"\s+", index_col=False, header=None,
                     names=['P(bar)', 'T(K)', 'logfo2', 'thing1', 'thing2', 'Rb'])

    # too many points to interpolate so find 4 nearest

    # get nearest pressures
    nT = 2023 - 1323  # number of temperature points
    tmp1 = p0 - df['P(bar)']  # want smallest nonnegative number so p1 < p0
    idx1 = tmp1.mask(tmp1 < 0).argmin()

    tmp2 = df['P(bar)'] - p0  # want smallest nonnegative number so p0 < p2
    idx2 = tmp2.mask(tmp2 < 0).argmin() + nT

    df = df.iloc[idx1:idx2 + 1]

    # get nearest temperatures
    tmp1 = T0 - df['T(K)']  # want smallest nonnegative number so p1 < p0
    idx1 = tmp1.mask(tmp1 < 0).argmin()
    T1 = df['T(K)'].iloc[idx1]

    tmp2 = df['T(K)'] - T0  # want smallest nonnegative number so p0 < p2
    idx2 = tmp2.mask(tmp2 < 0).argmin()
    T2 = df['T(K)'].iloc[idx2]

    df = df[df['T(K)'].isin([T1, T2])]

    # interpolate in this range
    if verbose:
        print('starting interpolation over df length', len(df), '| p_min', df['P(bar)'].iloc[0], 'p_max',
              df['P(bar)'].iloc[-1], 'p0', p0,
              '| T_min', df['T(K)'].iloc[0], 'T_max', df['T(K)'].iloc[-1], 'T0', T0)

    x_in = df['P(bar)'].to_numpy()
    y_in = df['T(K)'].to_numpy()
    z_in = df['logfo2'].to_numpy()
    # f = interpolate.griddata((x_in, y_in),
    f = interpolate.interp2d(x=x_in, y=y_in, z=z_in)
    logfo2 = f(p0, T0)
    return np.squeeze(logfo2)


T0 = 1630.5
p0 = 20399
fo2 = read_qfm_os(T0, p0, path=px.perplex_path_default, fin='data_tables/fmqNl_fo2_oli.dat')
print('log10(fo2)', fo2)
