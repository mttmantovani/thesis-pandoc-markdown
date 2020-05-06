import matplotlib as mpl
mpl.use('pgf')

mpl_settings = {
    'text.usetex': False,
    'font.family': 'serif',
    'font.size': 11,
}

defaults = {
    'width': 6,
    'height': 3,
}

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams.update(mpl_settings)

from sys import argv
import os
import numpy as np

filename = os.path.splitext(argv[0])[0]

# Load and process data
x_FHeff, y_FHeff = np.loadtxt(argv[1], delimiter=',').T
x_FHfull, y_FHfull = np.loadtxt(argv[2], delimiter=',').T
x_n, y_n = np.loadtxt(argv[3], delimiter=',').T

fig, ax = plt.subplots(figsize=(defaults['width'], defaults['height']))
ins = inset_axes(ax,
                 width=2.,
                 height=1,
                 bbox_to_anchor=(0.97, 0.23),
                 bbox_transform=ax.transAxes,
                 loc=4)

ax.set_xscale('log')
ins.set_xscale('log')

ax.plot(x_FHfull, y_FHfull, ls='-', lw=1.5, c='red')
ax.plot(x_FHeff,
        y_FHeff,
        ls='none',
        marker='o',
        markeredgecolor='k',
        markerfacecolor='C1',
        markersize=5,
        clip_on=False,
        zorder=5)

ax.axvline(x=np.sqrt(5 / 8), c='gray', ls='--', lw=1)

ax.set_xlim(left=1e-1, right=12.44)
ax.set_xlabel(r'$\tilde{E}_J \Delta_0^2 Q\ [\omega_0]$')
ax.set_ylabel(r'$F_\mathrm{CP}$')

ins.plot(x_n, y_n, c='C0', ls='-', lw=1.5)
ins.set_ylabel(r'$F_n$')
ins.set_xlim(left=0.2)
ins.set_ylim(bottom=0.96)
ins.set_yticks([0.96, 0.98, 1])
ins.set_xlabel(r'$\tilde{E}_J \Delta_0^2 Q\ [\omega_0]$')

plt.tight_layout()

plt.savefig(f'{filename}.pgf')
