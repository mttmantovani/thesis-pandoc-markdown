import matplotlib as mpl
mpl.use('pgf')

mpl_settings = {
    'text.usetex': False,
    'font.family': 'serif',
    'font.size': 11,
}

defaults = {
    'width': 3,
    'height': 2.6,
}

import matplotlib.pyplot as plt

plt.rcParams.update(mpl_settings)

import os
from sys import argv

filename = os.path.splitext(argv[0])[0]

import numpy as np

# Load and process data
x_gamma, y_gamma = np.loadtxt(argv[1], delimiter=',').T
x_nav, y_nav = np.loadtxt(argv[2], delimiter=',').T

fig, ax = plt.subplots(figsize=(defaults['width'], defaults['height']))

ax.plot(x_nav, y_nav, ls='-', lw=1.5, label=r'$\langle n \rangle$')
ax.plot(x_gamma,
        y_gamma,
        ls='-',
        lw=1.5,
        label=r'$\Gamma_{\mathrm{CP}}\ [\kappa]$')
ax.axvline(x=0, c='gray', ls='--', lw=1)
ax.set_xlabel(r'$\delta\ [10^{-3}\omega_0]$')
ax.set_xlim(-1.5, 1.5)
ax.legend(loc=2, frameon=False)

plt.tight_layout()

plt.savefig(f'{filename}.pgf')
