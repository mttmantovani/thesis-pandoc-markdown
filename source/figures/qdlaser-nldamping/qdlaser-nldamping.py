from sys import argv
import os
import matplotlib as mpl

mpl.use('pgf')

preamble = [
    r'\usepackage{pgf}', 
    r'\usepackage{amssymb}', 
    r'\usepackage{amsmath}',
    r'\usepackage{newtxtext,newtxmath}', 
    r'\usepackage{fontenc}', 
    r'\usepackage{inputenc}', 
    r'\DeclareUnicodeCharacter{2212}{-}'
]

mpl_settings = {
    'text.usetex': False,
    'font.family': 'serif',
    'pgf.texsystem': 'pdflatex',
    'font.size': 10,
    'pgf.preamble': preamble
}

TEXTWIDTH = 6.29  # Textwidth in inches

defaults = {
    'width': 5.,
    'height': 8/2.54,
}

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams.update(mpl_settings)
import numpy as np

filename = os.path.splitext(argv[0])[0]

def load_data():
    dataDict = {}
    datapaths = argv[1:]

    for idx, data in enumerate(datapaths):
        dataset = np.loadtxt(data, delimiter=',')
        dataDict[idx] = dataset

    return dataDict

def plot_data(**kwargs):

    # Load data
    dataDict = load_data()

    # Define figure layout
    fig = plt.figure(figsize=(defaults['width'], defaults['height']), linewidth=1)

    ax = fig.add_subplot()

    # Plot
    ax.plot(dataDict[0][:,0], -dataDict[0][:,2], c='royalblue', label=r'$\Gamma = 0.01\omega_0$')
    ax.plot(dataDict[0][:,0], -dataDict[0][:,4], c='turquoise', ls='--', label=r'$\Gamma=0.01\omega_0$ (RWA)')
    ax.plot(dataDict[0][:,0], -dataDict[0][:,1], c='darkgreen', label=r'$\Gamma = 0.1\omega_0$')
    ax.plot(dataDict[0][:,0], -dataDict[0][:,3], c='limegreen', ls='--', label=r'$\Gamma = 0.1\omega_0$ (RWA)')
    ax.axhline(y=1e-4, ls=':', color='k')
    ax.plot(12.7258, 1e-4, marker='o', c='darkgreen')
    ax.plot(19.4234, 1e-4, marker='o', c='darkgreen')
    ax.plot(16.629, 1e-4, marker='o', c='red', markerfacecolor='None')

    # Adjustments
    ax.set_ylim(0,3e-4)
    ax.set_xlim(0,40)
    ax.set_yticks(np.linspace(0,3e-4,7))
    ax.set_yticklabels([0, '', 1, '', 2, '', 3])
    ax.set_xlabel(r'$\bar{A}$')
    ax.set_ylabel(r'$-\gamma_\mathrm{eff}\ [\times 10^{-4}\omega_0]$')
 
    # Legend
    ax.legend(loc='upper right', frameon=False)

    # Letters and annotations
    ax.annotate('unstable', xy=(16.629,1e-4), xytext=(0,-40), xycoords='data', textcoords='offset points', arrowprops=dict(facecolor='k', edgecolor='k', shrinkA=1, shrinkB=5, arrowstyle='-|>'), ha='center', va='center')
    ax.annotate('stable', xy=(12.7258,1e-4), xytext=(22,1.5e-4), xycoords='data', textcoords='data', color='white', arrowprops=dict(facecolor='k', edgecolor='k', shrinkA=1, shrinkB=5, arrowstyle='-|>'), ha='center', va='center')
    ax.annotate('stable', xy=(19.4234,1e-4), xytext=(22,1.5e-4), xycoords='data', textcoords='data', arrowprops=dict(facecolor='k', edgecolor='k', shrinkA=1, shrinkB=5, arrowstyle='-|>'), ha='center', va='center')
    ax.annotate(r'$\kappa=10^{-4}\omega_0$', xy=(40,1e-4), xytext=(-2,2), xycoords='data', textcoords='offset points', ha='right', va='bottom')

    plt.subplots_adjust(left=0.08,bottom=0.13, right=0.98, top=0.92, wspace=0.4, hspace=0.08 )
    plt.savefig(f'{filename}.pdf', dpi=400, edgecolor='None')

if __name__ == '__main__':
    plot_data()
