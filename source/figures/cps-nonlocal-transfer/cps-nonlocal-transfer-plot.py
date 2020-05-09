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
    'font.size': 11,
    'pgf.preamble': preamble
}

TEXTWIDTH = 6.29  # Textwidth in inches

defaults = {
    'width': 10/2.54,
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
            if nidx == 8:
                labelLst = eval(line[1:].strip())
        dataDict[idx] = dict((k, v) for k, v in zip(labelLst, dataset.T))
        dataDict[idx]['filename'] = data

    return parDict, varDict, dataDict

def plot_data(**kwargs):
    # Load data
    parDict, varDict, dataDict = load_data()

    # Define figure layout
    fig = plt.figure(figsize=(defaults['width'], defaults['height']))
    ax1 = fig.add_subplot()  # Current

    # Plot curves
    ax1.plot(dataDict[0]['x'], dataDict[0]['IL']/1e-4, lw=1., c='k', label=r'$\lambda=0.05\Gamma_S$')
    ax1.plot(dataDict[1]['x'], dataDict[1]['IL']/1e-4, lw=1., c='magenta', ls='--', label=r'$\lambda = 0.2\Gamma_S$')

    # Adjustments
    ax1.set_xlim(-4, 4)
    ax1.set_ylim(bottom=0)
    ax1.set_ylabel(r'$I_\alpha\ [e\Gamma]$')
    ax1.set_xlabel(r'$\epsilon\ [\Gamma_S]')
    ax1.set_yticks(np.linspace(0,0.5,5))
    ax1.set_yticklabels([0, '', '', '',0.5])
    ax1.set_xticks(range(-4,5))
    ax1.yaxis.set_label_coords(0., 0.25, transform=ax1.yaxis.get_ticklabels()[0].get_transform())

    # Legends
    ax1.legend(loc=3, bbox_to_anchor=(0.65, 0.8), labelspacing=0.1, 
               columnspacing=0.85, borderpad=0.35, handletextpad=0.1,
               ncol=1,handlelength=1.75, borderaxespad=0.,frameon=False)

    # Letters and annotations
    ax1.annotate('(b)', xy=(0.02, 0.94), xycoords='axes fraction')

    ax1.annotate(
    r'$\bar{\delta}\!=\!\omega_L\!-\!\omega_R$',
    xy=(0.63, 0.31),
    xytext=(1.4, 0.39),
    xycoords='data',
    textcoords='data',
    arrowprops=dict(arrowstyle='->', lw=0.4),
    bbox=dict(pad=-1, facecolor='none', edgecolor='none'),
    va='center',
    ha='center',
    fontsize=mpl_settings['font.size']-1,
    )

    ax1.annotate(
        r'$\bar{\delta}\!=\!\omega_R$',
        xy=(1.28, 0.305),
        xytext=(1.7, 0.305),
        xycoords='data',
        textcoords='data',
        arrowprops=dict(arrowstyle='->', lw=0.4),
        bbox=dict(pad=-1, facecolor='none', edgecolor='none'),
        verticalalignment='center',
        fontsize=mpl_settings['font.size']-1,
    )
    ax1.annotate(
        r'$\bar{\delta}\!=\!\omega_L$',
        xy=(2.36, 0.21),
        xytext=(3, 0.23),
        xycoords='data',
        textcoords='data',
        arrowprops=dict(arrowstyle='->', lw=0.4),
        bbox=dict(pad=-1, facecolor='none', edgecolor='none'),
        verticalalignment='center',
        fontsize=mpl_settings['font.size']-1,
    )
    ax1.annotate(
        r'$\bar{\delta}\!=\!2\omega_R$',
        xy=(2.88, 0.137),
        xytext=(3., 0.17),
        xycoords='data',
        textcoords='data',
        arrowprops=dict(arrowstyle='->', lw=0.4),
        bbox=dict(pad=-1, facecolor='none', edgecolor='none'),
        verticalalignment='center',
        fontsize=mpl_settings['font.size']-1,
    )



    plt.subplots_adjust(left=0.085,bottom=0.19, right=0.99, top=0.99, wspace=0.4, hspace=0.10 )
    plt.savefig(f'{filename}.pdf', facecolor='none', edgecolor='none', dpi=400)


if __name__ == '__main__':
    plot_data()
