from __future__ import unicode_literals
import numpy as np
import matplotlib.pyplot as plt
from useful_and_bespoke import dark_background, cornertext
import matplotlib.ticker as ticker
from datetime import date
from matplotlib import rc
from matplotlib.pyplot import rcParams
import parameters as p
from matplotlib import use as mpluse
import matplotlib.lines as mlines
import matplotlib

leg_loc = 'right'
rc('text', usetex=True)  # turn off for running over ssh
# rcParams['font.family'] = 'serif'
# rcParams['font.serif'] = 'CMU Serif'
legcolor='xkcd:off white'

today = date.today().strftime("%b-%d-%Y")
fig_path = '/home/claire/Desktop/rw-poster/'
rho_c = 2700
rho_w = 1000
M_E = 5.972e24  # earth mass in kg
R_E = 6371e3
g_E = 9.807  # earth grav
Y = 100e6
labelsize = 45  ##32
ticksize = 25
legsize = 25
lw = 5
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
rcParams['legend.title_fontsize'] = legsize


def plot_interior_capacity():
    c_int = 'xkcd:seafoam'
    Mp = np.array([0.1, 0.3, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5])
    med = np.array([0.36100584050789875, 1.1029287868637658, 1.5439122274024288, 2.4222631087220257, 3.213658563226869,
                    3.9726802221774142, 4.7222290938376155, 5.4505475983686225, 6.889111728730112,
                    8.305615141504067])  # medians in OM
    # std = np.array([0.07815052540158095, 0.3442004713636262, 0.902686817841038, 2.3024633913595154 , 3.3241903236080406,  4.2461395182413355, 5.093035310490696 ,  5.917910897159271 , 7.508831299776319 ,  9.06256520708619])  # standard deviations in OM
    med_wmf = med * p.TO / (Mp * p.M_E)
    # std_wmf = std * p.TO / (Mp * p.M_E)
    # up_wmf = med_wmf + std_wmf
    # lo_wmf = med_wmf - std_wmf
    lo = np.array([0.3154103080521038, 0.9918987071009212, 1.4130768816205463, 2.2483215268475187, 3.0177064265639078,
                   3.754590207348003, 4.476295226154339, 5.180187129302038, 6.547143954319714, 7.8912883520545565])
    up = np.array([0.416770658029093, 1.256317373359301, 1.7402325535580863, 2.655937233969786, 3.4734790571911316,
                   4.241381784970886, 5.019385870890552, 5.771796427861774, 7.289798571154144, 8.792363887608982])
    lo_wmf = lo * p.TO / (Mp * p.M_E)
    up_wmf = up * p.TO / (Mp * p.M_E)
    # ax.plot(Mp, med_wmf, c=c_int, lw=2, ls='--')
    for y in (up, lo):
        ax.plot(Mp, y, c=c_int, lw=0.5)
    h = ax.fill_between(Mp, lo_wmf, up_wmf, fc=c_int, ec=c_int, alpha=0.5, hatch='///', zorder=10,
                        label='If the entire mantle water capacity\nis brought to the surface')

    # leg_mantle = ax.legend(handles=[h, ],
    #                        ncol=1, frameon=False, fontsize=legsize,
    #                        # title=r'\textbf{Observed ocean masses}',
    #                        bbox_to_anchor=(1.04, 0.2),  # (0, 1.02, 1, 0.2),
    #                        loc='lower left',
    #                        borderaxespad=0)
    # ax.add_artist(leg_mantle)
    return h


def plot_trappist_wmf(cmf='all', err_c='0.8'):
    M_Tre = 0.692
    M_Trf = 1.039
    M_Trg = 1.321
    err_m = 's'
    err_kwargs = {'elinewidth': 1, 'capsize': 5, 'ms': 7, 'lw': 0, 'c': err_c, 'marker': err_m}
    lim_kwargs = err_kwargs.copy()
    lim_kwargs.update({'marker': '_'})
    annosize = legsize  # - 5
    anno_c = err_c  # ms
    if cmf == 'all':
        w_Tre_av = 0.3e-2
        w_Tre = 2.9e-2 + 1.7e-2
        e_Tre = [[w_Tre - 3e-5], [2.9e-2 + 1.7e-2 - w_Tre]]  # wmf error

        w_Trf_av = 1.9e-2
        w_Trf = 4.5e-2 + 1.8e-2
        e_Trf = [[w_Trf - 3e-5], [4.5e-2 + 1.8e-2 - w_Trf]]

        w_Trg_av = 0.72e-2
        w_Trg = 6.4e-2 + 2e-2
        e_Trg = [[w_Trg - 3e-5], [6.4e-2 + 2e-2 - w_Trg]]

        # TRAPPIST-1e
        ax.errorbar(M_Tre, w_Tre, yerr=e_Tre, xerr=None, uplims=True, zorder=3, **lim_kwargs)
        ax.errorbar(M_Tre, w_Tre_av, yerr=None, xerr=None, zorder=3, **err_kwargs)
        ax.annotate('e', (M_Tre + 0.05, w_Tre_av + 1e-3), fontsize=annosize,
                    color=anno_c, ha='center', va='bottom')

        # TRAPPIST-1f
        ax.errorbar(M_Trf, w_Trf, yerr=e_Trf, xerr=None, uplims=True, zorder=3, **lim_kwargs)
        ax.errorbar(M_Trf, w_Trf_av, yerr=None, xerr=None,zorder=3, **err_kwargs)
        ax.annotate('f', (M_Trf + 0.07, w_Trf_av + 1e-3), fontsize=annosize,
                    color=anno_c, ha='center', va='bottom')

        # TRAPPIST-1g
        ax.errorbar(M_Trg, w_Trg, yerr=e_Trg, xerr=None, uplims=True, zorder=3, **lim_kwargs)
        ax.errorbar(M_Trg, w_Trg_av, yerr=None, xerr=None, label='DUMMY', zorder=3, **err_kwargs)
        ax.annotate('g', (M_Trg + 0.15, w_Trg_av - 1e-3), fontsize=annosize,
                    color=anno_c, ha='center', va='top')

    elif cmf == 25:
        w_Tre_av = 0.3e-2
        e_Tre = [[0.3e-2], [1.8e-2]]  # wmf error

        w_Trf_av = 1.9e-2
        e_Trf = [[1.3e-2], [1.5e-2]]

        w_Trg_av = 3.5e-2
        e_Trg = [[1.3e-2], [1.6e-2]]

        # TRAPPIST-1e
        ax.errorbar(M_Tre, w_Tre_av + 1.8e-2, yerr=[[w_Tre_av + 1.8e-2 - 3e-5], [1.8e-2]], xerr=None, uplims=True, zorder=3, **lim_kwargs)
        ax.errorbar(M_Tre, w_Tre_av, yerr=None, xerr=None, zorder=3, **err_kwargs)
        ax.annotate('e', (M_Tre + 0.05, w_Tre_av + 1e-3), fontsize=annosize, zorder=3,
                    color=anno_c, ha='center', va='bottom')

        # TRAPPIST-1f
        ax.errorbar(M_Trf, w_Trf_av, yerr=e_Trf, xerr=None, zorder=3, **err_kwargs)
        ax.annotate('f', (M_Trf + 0.07, w_Trf_av + 1e-3), fontsize=annosize,zorder=3,
                    color=anno_c, ha='center', va='bottom')

        # TRAPPIST-1g
        ax.errorbar(M_Trg, w_Trg_av, yerr=e_Trg, xerr=None, label='DUMMY', zorder=3,**err_kwargs)
        ax.annotate('g', (M_Trg + 0.15, w_Trg_av - 1e-3), fontsize=annosize,zorder=3,
                    color=anno_c, ha='center', va='top')


def plot_ss():
    # Mars
    M_Mars = 6.39e23
    R_Mars = 3389.5e3
    h_peak_Mars = 21.9e3  # Olympus Mons

    M_E = p.M_E
    R_E = p.R_E
    h_peak_E = 8840  # Mt Everest

    # Venus
    M_V = 4.867e24
    R_V = 6051.8e3
    h_peak_V = 11000  # Maxwell Montes

    c_ss = 'xkcd:cobalt'
    p1 = plt.scatter(M_V / M_E, 4.727695769579759e+18 * 1000 / M_V, c=c_ss, marker='*', s=80, zorder=1000)
    p2 = plt.scatter(1, 4.641828875736904e+18 * 1000 / M_E, c=c_ss, marker='*', s=80, zorder=1000)
    p3 = plt.scatter(M_Mars / M_E, 4.013894774478276e+18 * 1000 / M_Mars, c=c_ss, marker='*', s=80, zorder=1000)
    plt.annotate(r'\textbf{Venus}', (M_V / M_E + 0.07, 4.727695769579759e+18 * 1000 / M_V + 0.0003), c='w',
                 fontsize=legsize,
                 bbox=dict(boxstyle="square,pad=0.1", fc=c_ss, ec=None, lw=0.5)
                 )  # Venus  u'\u2640', chr(0x263f+1)
    plt.annotate(r'\textbf{Earth}', (1 + 0.1, 4.641828875736904e+18 * 1000 / M_E - 0.00008), c='w', fontsize=legsize,
                 va='top',
                 bbox=dict(boxstyle="square,pad=0.1", fc=c_ss, ec=None, lw=0.5)
                 )  # Earth chr(0x263f+2)
    plt.annotate(r'\textbf{Mars}', (M_Mars / M_E + 0.01, 4.013894774478276e+18 * 1000 / M_Mars), c='w',
                 fontsize=legsize,
                 bbox=dict(boxstyle="square,pad=0.1", fc=c_ss, ec=None, lw=0.5)
                 )  # Mars chr(0x263f+3)

    leg_ss = ax.legend([p1, ],
                       ['Solar system planets', ], ncol=1, frameon=False, fontsize=legsize,
                       title=r'\textbf{Observed topographic capacity}',
                       bbox_to_anchor=(1.04, 0.6), borderaxespad=0,  labelcolor=[legcolor])
    leg_ss._legend_box.align = "left"


def radius_zeng(M_p, CMF=0.3):
    # input mass in kg, output radius in m
    # applicable to M_E <= 8 and CMF <= 0.4
    #     print('using Zeng radius model')
    return (1.07 - 0.21 * CMF) * (M_p / M_E) ** (1 / 3.7) * R_E


def grav(M, R):
    """Calculate acceleration due to gravity on a point mass in m s^-2"""
    return 6.674e-11 * M / R ** 2


def vol(R):
    return 4 / 3 * np.pi * R ** 3


def h_peak_rock(M=None, R=None, Y=100e6, rho_c=2700, C=1 / 2, **kwargs):
    # C is 1/3 to 1/2 - min stress difference supported elastically underneath load (Jeffreys' theorem)
    g = grav(M, R)
    return (C ** -1 * Y) / (rho_c * g)


def h_peak_dt(M=None, peak_scale=3.5, **kwargs):  # M in kg
    h_rms = (M / M_E) ** -0.44 * 10 ** 2.26
    return h_rms * peak_scale


def h_peak_dt_cold(M=None, peak_scale=3.5, **kwargs):  # M in kg
    h_rms = (M / M_E) ** -0.3434899 * (10 ** 3.1794847)
    return h_rms * peak_scale


def h_peak_dt_hot(M=None, peak_scale=3.5, **kwargs):  # M in kg
    h_rms = (M / M_E) ** -0.74629234 * (10 ** 2.45476118)
    return h_rms * peak_scale


def wmf_max(M_p, R_p, rho_m=3500, rho_w=1000, h_peak_fn=None, **kwargs):
    # M_p R_p in SI units
    V_p = vol(R_p)
    h = h_peak_fn(M=M_p, R=R_p, **kwargs)  # in m
    # m_cap = rho_w * (V_p - 4/3*np.pi*(R_p - h)**3)
    m_cap = rho_w * rho_m / (rho_m - rho_w) * 4 * np.pi * R_p ** 2 * h  # accounts for water loading
    return m_cap / M_p  # as mass fraction


# give an M, what is the max water content that would make it a waterworld?
masses = np.logspace(np.log10(0.1), np.log10(5), num=10)
gravs = np.ones_like(masses)
wmfs_200 = np.ones_like(masses)
wmfs_100 = np.ones_like(masses)
wmfs_50 = np.ones_like(masses)
wmfs_dt_cold = np.ones_like(masses)
wmfs_dt_hot = np.ones_like(masses)
for ii, m in enumerate(masses):
    M_p = m * M_E
    R_p = radius_zeng(M_p)  # in km
    gravs[ii] = grav(M_p, R_p) / g_E

    wmfs_100[ii] = wmf_max(M_p, R_p, h_peak_fn=h_peak_rock, Y=100e6, rho_c=2700)
    wmfs_50[ii] = wmf_max(M_p, R_p, h_peak_fn=h_peak_rock, Y=50e6, rho_c=2700)
    wmfs_200[ii] = wmf_max(M_p, R_p, h_peak_fn=h_peak_rock, Y=200e6, rho_c=2700)

    wmf_dt_cold = wmf_max(M_p, R_p, h_peak_fn=h_peak_dt_cold)
    wmfs_dt_cold[ii] = wmf_dt_cold

    wmf_dt_hot = wmf_max(M_p, R_p, h_peak_fn=h_peak_dt_hot)
    wmfs_dt_hot[ii] = wmf_dt_hot
wmfs_dt_avg = np.log10((10 ** wmfs_dt_hot + 10 ** wmfs_dt_cold) / 2)

""" DO PLOTTING """

# plot_ss()
plot_trappist_wmf(cmf=25, err_c='xkcd:off white')

handles, labels = ax.get_legend_handles_labels()
labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))

# add mantle water capacity
h_int = plot_interior_capacity()

p4 = plt.scatter(1, 1.4e21 / M_E, c='k', ec='k', marker='$\oplus$', s=500, zorder=1000)  # c='0.3' # Earth ocean
# ax.annotate('Earth, modern ocean', (1, 1.4e21 / M_E), fontsize=legsize,
#             color=anno_c, ha='center', va='bottom')

if leg_loc == 'top':
    bbox = (0.6, 1.06)
    loc = 'lower left'
elif leg_loc == 'right':
    bbox = (1.04, 0.7)
    loc = 'upper left'
leg_trapp = ax.legend([
    handles
    ,p4,
    h_int
],
    [
        'TRAPPIST-1 system, 25\% Fe',
        'Earth, modern ocean',
    'Entire mantle water capacity'
        ],
    ncol=1, frameon=False, fontsize=legsize,
    title=r'\textbf{Ocean mass fractions}',
    bbox_to_anchor=bbox,  # (0, 1.02, 1, 0.2),
    loc=loc,
    borderaxespad=0, labelcolor=[legcolor]*3)

c = '#b13c02ff'
c_cr = '#c5927fff'  # 'xkcd:sand'
c_dt = '#cc7832'  # 'xkcd:greey green'

# # l1, = plt.plot(masses, wmfs_200, c=c_cr, lw=6, label='Crust strength, 200 MPa', zorder=2)
l2, = plt.plot(masses, wmfs_100, c=c_cr, lw=6, label='Crust strength, 100 MPa', zorder=2)
# # l3, = plt.plot(masses, wmfs_50, c=c_cr, lw=1, label='Crust strength, 50 MPa', zorder=1)
l2 = mlines.Line2D([], [], color=c_cr, marker=None, lw=6, label='Crust strength, 100 MPa')  # dummy for leg

# l6 = plt.fill_between(masses, y1=wmfs_dt_hot, y2=wmfs_dt_cold, fc=c_dt, alpha=0.5, zorder=1)
# l6 = plt.plot(masses, wmfs_dt_avg, c=c_dt, lw=3, zorder=2, visible=False)
plt.errorbar(masses, wmfs_dt_hot, yerr=wmfs_dt_hot * 2e-1, c=c_dt, ls='-', lw=3, zorder=2,
             lolims=True)  # arbitrary constant error for plot
plt.errorbar(masses, wmfs_dt_cold, yerr=wmfs_dt_cold * 2e-1, c=c_dt, ls='-', lw=3, zorder=2, uplims=True)
l6 = mlines.Line2D([], [], color=c_dt, marker=None, lw=3, label='Dynamic topography limit')  # dummy for leg

plt.fill_between(masses, wmfs_dt_cold, [1] * len(masses), zorder=0, fc='#131a30', alpha=0.6)
plt.fill_between(masses, [1e-10] * len(masses), wmfs_dt_hot, zorder=0, fc='xkcd:stone', alpha=0.5)
ax = cornertext(ax, 'WATER\nWORLD', pos='top right', size=legsize+7, c='xkcd:off white', pad=0.04)  # pad=0.05
ax = cornertext(ax, 'BLUE\nMARBLE?', pos='bottom left', size=legsize+7, c='xkcd:off white', pad=0.04)

leg_trapp._legend_box.align = "left"

# ax.add_artist(leg_ss)
ax.add_artist(leg_trapp)

if leg_loc == 'top':
    bbox = (-0.15, 1.06)
    loc = 'lower left'
elif leg_loc == 'right':
    bbox = (1.04, 1)
    loc = 'upper left'
leg_scale = ax.legend(handles=[
    # # l1,
    l2,
    # # l3,
    l6
],
    ncol=1, frameon=False, fontsize=legsize, loc=loc,
    title=r'\textbf{Theoretical topographic capacity}',
    bbox_to_anchor=bbox, borderaxespad=0, labelcolor=[legcolor]*len(handles))
leg_scale._legend_box.align = "left"
ax.add_artist(leg_scale)



# dummy
# h = ax.fill_between([], [], [], fc=[], ec=[], alpha=0.5, hatch='///', zorder=10,
#                     label='If the entire mantle water capacity\nis brought to the surface')
#
# leg_mantle = ax.legend(handles=[h, ],
#                        ncol=1, frameon=False, fontsize=legsize,
#                        # title=r'\textbf{Observed ocean masses}',
#                        bbox_to_anchor=(1.04, 0.2),  # (0, 1.02, 1, 0.2),
#                        loc='lower left',
#                        borderaxespad=0)
# ax.add_artist(leg_mantle)

plt.xlabel(r'Planet mass $(M_\oplus$)', fontsize=labelsize, labelpad=20)
plt.ylabel('Waterworld limit\n(ocean mass fraction)', fontsize=labelsize, labelpad=20)
plt.loglog()

xticks = [0.1, 1, 5]
ax.set_ylim(2E-5, 1E-1)
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.tick_params(axis='x', labelsize=ticksize)
ax.tick_params(axis='y', labelsize=ticksize)

# plt.errorbar(masses[8], wmfs_dt_avg[8], yerr=[[wmfs_dt_avg[8] - wmfs_dt_cold[8]], [wmfs_dt_hot[8] - wmfs_dt_cold[8]]], c='k', ls='-', lw=1, zorder=2)

# change all spines
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(4)

# increase tick width
ax.tick_params(width=4)


plt.setp(leg_scale.get_title(), color='xkcd:off white')
plt.setp(leg_trapp.get_title(), color='xkcd:off white')
fig, ax = dark_background(fig, ax)

plt.savefig(fig_path + 'ww_lim.png', bbox_inches='tight', dpi=600,
            transparent=True,
            # facecolor=fig.get_facecolor()
            )
plt.show()
