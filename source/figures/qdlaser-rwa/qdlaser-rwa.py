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
    'font.size': 9,
    'pgf.preamble': preamble
}

TEXTWIDTH = 6.29  # Textwidth in inches

defaults = {
    'width': 7.7/2.54,
    'height': 6.5/2.54,
}

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams.update(mpl_settings)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np

filename = os.path.splitext(argv[0])[0]

def load_data():
    dataDict = {}
    parDict = {}
    varDict = {}
    datapaths = argv[1:]

    for idx, data in enumerate(datapaths):
        labelLst = list()
        dataset = np.loadtxt(data)
        for nidx, line in enumerate(open(data, 'r')):
            if nidx == 2:
                parDict[idx] = eval(line[1:].strip())
            if nidx == 4:
                varDict[idx] = eval(line[1:].strip())
            if nidx == 10:
                labelLst = eval(line[1:].strip())
        dataDict[idx] = dict((k, v) for k, v in zip(labelLst, dataset.T))
        dataDict[idx]['filename'] = data

    return parDict, varDict, dataDict

def plot_data(**kwargs):
    # Load data
    parDict, varDict, dataDict = load_data()

    # Define figure layout
    fig = plt.figure(figsize=(defaults['width'], defaults['height']))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.3, 0.7])
    ax1 = fig.add_subplot(gs[0, 0])  # Current
    ax2 = fig.add_subplot(gs[1, 0])  # Photon occupation
    axins = inset_axes(ax2, loc=1, width=1.3, height=0.65,
                       bbox_to_anchor=(0.8, 1), borderpad=0.5,
                       bbox_transform=ax2.transAxes
            )

    # Plot curves
    p = parDict # shorthand
    nB = lambda x: 1. / (np.exp(x) - 1)     # Bose function
    ABS = lambda i, sgn, n, m: 0.5 * sgn * np.sqrt((n * p[i]['wL'] + m * p[i]['wR'])**2 - 2 * p[i]['GS']**2) # Andreev addition energy

    ax1.plot(dataDict[0]['eps'], dataDict[0]['IL']/1e-3, lw=1., c='k', label=r'$k_BT=5\omega$')
    ax1.plot(dataDict[1]['eps'], dataDict[1]['IL']/1e-3, lw=1., c='orange', ls='--', label=r'$k_BT = 0$')
    ax2.plot(dataDict[0]['eps'], dataDict[0]['nOscL'], lw=1., c='k')
    ax2.axhline(y=nB(p[0]['wL'] / p[0]['TL']), color='0.6', ls=':', lw=1.)
    ax2.axvline(x=ABS(0,1,1,0), ymin=0.08, ymax=0.7, color='b', lw=0.5)

    axins.plot(dataDict[2]['GS'], dataDict[2]['nOscL']/nB(p[2]['wL']/p[2]['TL']), lw=1., c='k', label=r'$\Gamma=2\times 10^{-4}\omega$')
    axins.plot(dataDict[3]['GS'], dataDict[3]['nOscL']/nB(p[3]['wL']/p[3]['TL']), lw=1., c='ForestGreen', ls='--', label=r'$\Gamma=10^{-3}\omega$')


    # Adjustments
    ax1.set_xlim(-3, 3)
    ax1.set_ylim(bottom=0)
    ax1.set_ylabel(r'$I_\alpha\ [e\Gamma]$')
    ax1.set_yticks([0, 0.25, 0.5])
    ax1.set_yticklabels([0, '', 0.5])
    ax1.set_xticks(range(-3,4))
    ax1.set_xticklabels([])
    ax1.yaxis.set_label_coords(0., 0.25, transform=ax1.yaxis.get_ticklabels()[0].get_transform())

    ax2.set_xlim(-3, 3)
    ax2.set_ylim(bottom=0, top=14)
    ax2.set_xlabel(r'$\epsilon\ [\Gamma_S]')
    ax2.set_ylabel(r'$\bar{n}_\alpha$')
    ax2.set_yticks(range(0,16,2))
    ax2.set_yticklabels([0, '', 4, '', 8, '', 12, ''])
    ax2.yaxis.set_label_coords(0., 6, transform=ax2.yaxis.get_ticklabels()[0].get_transform())

    axins.set_xlim(0, 0.7)
    axins.set_ylim(0, 1)
    axins.set_xticks(np.linspace(0,0.7,3))
    axins.set_xticklabels([0,'',0.7])
    axins.set_yticks(np.linspace(0,1,3))
    axins.set_yticklabels([0,'',1])
    axins.set_xlabel(r'$\Gamma_S\ [\omega]$')
    axins.xaxis.set_label_coords(0.35,0, transform=axins.xaxis.get_ticklabels()[0].get_transform())
    axins.set_ylabel(r'$\frac{\bar{n}_{\alpha}\left(\epsilon_{c}\right)}{\bar{n}_{\alpha}^{\text{th}}}$')
    axins.yaxis.set_label_coords(0.,0.5, transform=axins.yaxis.get_ticklabels()[0].get_transform())


    # Legends
    ax1.legend(loc=3, bbox_to_anchor=(0.7, 0.5), labelspacing=0.1, 
               columnspacing=0.85, borderpad=0.35, handletextpad=0.1,
               ncol=1,handlelength=1.75, borderaxespad=0.,frameon=False)
    axins.legend(
        loc=3,
        bbox_to_anchor=(0.1, 0.45),
        labelspacing=0.1,
        columnspacing=0.85,
        borderpad=0.35,
        handletextpad=0.1,
        ncol=1,
        handlelength=1.4,
        borderaxespad=0.,
        frameon=False)


    # Letters and annotations
    ax1.annotate('(b)', xy=(0.01, 0.84), xycoords='axes fraction')
    ax2.annotate('(c)', xy=(0.01, 0.93), xycoords='axes fraction')
    ax2.annotate('cooling', xy=(2.398-0.2, 7.5), xycoords='data', rotation=90, color='blue', ha='center', va='center')
    ax2.annotate(r'$\epsilon_c$', xy=(2.398+0.2, 7.5), xycoords='data', color='blue', ha='center', va='center')


    plt.subplots_adjust(left=0.085,bottom=0.15, right=0.99, top=0.99, wspace=0.4, hspace=0.10 )
    plt.savefig(f'{filename}.pdf', facecolor='none', edgecolor='none', dpi=400)


if __name__ == '__main__':
    plot_data()
    