import matplotlib as mpl
mpl.use('pgf')

mpl_settings = {
    'text.usetex': False,
    'font.family': 'serif',
    'font.size': 11,
}

defaults = {
    'width': 6,
    'height': 3.5,
}

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams.update(mpl_settings)

from sys import argv
import os
import numpy as np

filename = os.path.splitext(argv[0])[0]

# Load and data
x_gamma, y_gamma = np.loadtxt(argv[1], delimiter=',').T
x_Heff, y_Heff = np.loadtxt(argv[2], delimiter=',').T
x_Hfull, y_Hfull = np.loadtxt(argv[3], delimiter=',').T
x_oEJ2, y_oEJ2 = np.loadtxt(argv[4], delimiter=',').T
x_oEJ4, y_oEJ4 = np.loadtxt(argv[5], delimiter=',').T

# Colors

fig = plt.figure(figsize=(defaults['width'], defaults['height']))
gs = gridspec.GridSpec(2, 1, height_ratios=[1.5, 1])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')

ax1.plot(x_Hfull, y_Hfull, ls='-', lw=1.5, label=r'$H(t)$')
ax1.plot(x_Heff, y_Heff, ls='--', lw=1.5, c='r', label=r'$H_\mathrm{eff}$')
ax1.plot(x_oEJ2, y_oEJ2, ls=':', lw=1, c='k', label=r'$\mathcal{O}(E_J^2)$')
ax1.plot(x_oEJ4,
         y_oEJ4,
         ls='-.',
         lw=1,
         c='limegreen',
         label=r'$\mathcal{O}(E_J^4)$',
         zorder=0)
ax1.axvline(x=np.sqrt(5 / 8), c='gray', ls=':')
ax1.set_xticklabels([])
ax1.set_yticks([1e-5, 1])
ax1.set_xlim(0.2, 30)
ax1.set_ylabel(r'$\langle n \rangle$')
ax1.legend(loc=2, frameon=False)

ax2.plot(x_gamma, y_gamma, ls='-', c='C1', lw=1.5)
ax2.axvline(x=np.sqrt(5 / 8), c='gray', ls=':')
ax2.axhline(y=2, c='r', ls='--')
ax2.set_xlim(0.2, 30)
ax2.set_xlabel(r'$\tilde{E}_J \Delta_0^2 Q\ [\omega_0]$')
ax2.set_ylabel(r'$\frac{\Gamma_\mathrm{CP}}{\kappa \langle n \rangle}$')

plt.tight_layout()

plt.savefig(f'{filename}.pgf')
