import matplotlib as mpl
mpl.use('pgf')

mpl_settings = {'text.usetex': False,
                'font.family': 'serif',
                'font.size': 11,
               }

defaults = {'width': 5,
            'height': 3,
           }

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams.update(mpl_settings)

import os
from sys import argv

filename = os.path.splitext(argv[0])[0]

import numpy as np

# Load and process data
x_g2, y_g2 = np.loadtxt(argv[1], delimiter=',').T
x_M12, y_M12 = np.loadtxt(argv[2], delimiter=',').T


fig, ax = plt.subplots(figsize=(defaults['width'], defaults['height']))
ins = inset_axes(ax, width=2, height=0.7,
                 bbox_to_anchor=(0.85, 1),
                 bbox_transform=ax.transAxes,
                 )

ax.plot(x_g2, y_g2, ls='-', lw=1.5)
ax.axhline(y=0, c='gray', ls='--', lw=1)
ax.set_xlabel(r'$\Delta_0$')
ax.set_ylabel(r'$g^{(2)}(0)$')
ax.set_xlim(0.1, 1.5)
ax.set_ylim(top=1.4)
ax.set_xticks([0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4])
ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])

ins.plot(x_M12, y_M12, c='r', ls='-', lw=1.5)
ins.axhline(y=0, c='gray', ls='--', lw=1)
ins.set_xlabel(r'$\Delta_0$')
ins.set_ylabel(r'$M_{12}\ \left[\frac{\tilde{E}_J^2}{4\omega_J}\right]$')
ins.set_xlim(0., 1.5)
ins.set_ylim(bottom=-1)


plt.tight_layout()

plt.savefig(f'{filename}.pgf')


