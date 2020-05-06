import matplotlib as mpl
mpl.use('pgf')

preamble = [
    r'\usepackage{pgf}',
    r'\usepackage{amssymb}',
    r'\usepackage{amsmath}',
    r'\usepackage{newtxtext,newtxmath}'
]

mpl_settings = {
    'text.usetex': False,
    'font.family': 'serif',
    'pgf.texsystem': 'pdflatex',
    'font.size': 11,
    'pgf.preamble': preamble
}

defaults = {
    'width': 6,
    'height': 2.5,
}

import matplotlib.pyplot as plt

plt.rcParams.update(mpl_settings)

from sys import argv
import os
import numpy as np

filename = os.path.splitext(argv[0])[0]

# Define analytical function for Andreev current

def I(e, GS):
    return GS**2/(e**2 + 2*GS**2)

eps = np.linspace(-4e4, 4e4, 201)
current1 = I(eps, 1e4)
current2 = I(eps, 2e4)

# Plot
fig, ax = plt.subplots(figsize=(defaults['width'], defaults['height']))

ax.plot(eps, current1, ls='-', lw=1.5, c='C0', label=r'$\Gamma_S = 1\ [10^4\ \Gamma]$')
ax.plot(eps, current2, ls='--', lw=1.5, c='C1', label=r'$\Gamma_S =2\ [10^4\ \Gamma]$')

ax.legend(loc='best', frameon=False)

ax.set_xlim(left=eps.min(), right=eps.max())
ax.set_xlabel(r'$\bar{\epsilon}\ [10^4\ \Gamma]$')
ax.set_ylabel(r'$I_\alpha\ [e\Gamma]$')

ax.set_xticklabels([1e4*i for i in range(-4,5)])
ax.set_xticklabels(range(-4,5))
plt.tight_layout()

plt.savefig(f'{filename}.pdf')
