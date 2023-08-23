import matplotlib.pyplot as plt
import numpy as np

with plt.xkcd():
    # Based on "Stove Ownership" from XKCD by Randall Munroe
    # https://xkcd.com/418/

    fig = plt.figure()
    ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
    ax.spines[['top', 'right']].set_visible(False)
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 1])

    prop = dict(arrowstyle="->,head_width=0.4,head_length=0.8", ec='k', fc='w', lw=2,
                shrinkA=4, shrinkB=4)

    plt.annotate("Relative\nsurface water\ncapacity", xy=(0.9, 0.1), xytext=(0.9, 0.8), arrowprops={**prop, 'relpos': (0.5, 0)}, ha='center')
    plt.annotate("Relative\ninterior water\ncapacity", xy=(0.9, 0.1), xytext=(0.1, 0.8), arrowprops=prop)
    plt.annotate("Water\noutgassing", xy=(0.9, 0.1), xytext=(0.1, 0.1), arrowprops={**prop, 'relpos': (0, 0.5)}, va='center')


    #
    # data = np.ones(100)
    # data[70:] -= np.arange(30)
    #
    # ax.annotate(
    #     'THE DAY I REALIZED\nI COULD COOK BACON\nWHENEVER I WANTED',
    #     xy=(70, 1), arrowprops=dict(arrowstyle='->'), xytext=(15, -10))
    #
    # ax.plot(data)

    ax.set_xlabel('Mg/Si', labelpad=15)
    ax.set_ylabel('Planet Mass', labelpad=10)
    # fig.text(
    #     0.5, 0.05,
    #     '"Stove Ownership" from xkcd by Randall Munroe',
    #     ha='center')

    fig.savefig('/home/claire/Works/thesis/v0/Conclusion/Figs/summary.pdf', bbox_inches='tight')

plt.show()