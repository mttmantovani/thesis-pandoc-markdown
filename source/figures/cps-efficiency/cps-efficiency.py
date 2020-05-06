from sys import argv
import os
import matplotlib as mpl

mpl.use('pgf')

preamble = [
    r'\usepackage{pgf}', r'\usepackage{amssymb}', r'\usepackage{amsmath}',
    r'\usepackage{newtxtext,newtxmath}'
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
    'width': 6.4,
    'height': 2.5,
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
    gs = gridspec.GridSpec(1, 2, width_ratios=[0.3, 0.7])
    ax1 = fig.add_subplot(gs[0, 0])  # Local efficiency
    ax2 = fig.add_subplot(gs[0, 1])  # Nonlocal efficiency
    axins = inset_axes(
        ax2,
        loc=1,
        width=1,
        height=0.85,
        # bbox_to_anchor=(0.24, 0.19),
        borderpad=0.3,
        bbox_transform=ax2.transAxes)

    # Plot data
    ax1.plot(dataDict[0]['x'], abs(dataDict[0]['eta_L']), lw=1., c='orange')

    ax2.plot(dataDict[1]['x'], abs(dataDict[1]['eta_NL']), lw=1., c='g')

    nB = lambda x: 1. / (np.exp(x) - 1.)
    argnbL = float(parDict[1]['wL']) / float(parDict[1]['TOscL'])
    argnbR = float(parDict[1]['wR']) / float(parDict[1]['TOscR'])
    axins.plot(dataDict[1]['x'],
               dataDict[1]['nOscR'] / nB(argnbR),
               ls='-.',
               c='r',
               lw=1.,
               label=r'$n_R$')
    axins.plot(dataDict[1]['x'],
               dataDict[1]['nOscL'] / nB(argnbL),
               c='b',
               lw=1.,
               label='$n_L$')

    ax1.set_xlim(2.2, 2.6)
    ax1.set_ylim(0, 1.)
    ax1.set_xticks([2.2, 2.4, 2.6])
    ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax1.set_xlabel(r'$\epsilon\ [\Gamma_S]$')
    ax1.set_ylabel(r'$\eta_{\mathrm{loc}}^{(L)}$')

    ax2.set_xlim(0.706, 0.712)
    ax2.set_ylim(0, 1.)
    ax2.set_xticks([0.706, 0.709, 0.712])
    ax2.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax2.set_xlabel(r'$\epsilon\ [\Gamma_S]')
    ax2.set_ylabel(r'$\eta_{\mathrm{NL}}$')

    axins.set_xlim(0.706, 0.712)
    axins.set_ylim(0.7, 1.3)
    axins.set_xticks([0.706, 0.709, 0.712])
    axins.set_xticklabels(['', '', ''])
    #   axins.set_xticklabels(['0.706ins', '', '0.712ins'])
    axins.set_yticks([0.8, 1, 1.2])
    axins.set_xlabel(r'$\epsilon\ [\Gamma_S]$')
    axins.set_ylabel(r'$\frac{\bar{n}_\alpha}{n^{\mathrm{th}}_\alpha}$',
                     rotation=0,
                     labelpad=10)

    # This moves the xaxis label closer to the axis
    #axins.xaxis.set_label_coords(0.5, -0.1)
    # Legend

    axins.legend(frameon=False,
                 loc=1,
                 bbox_to_anchor=(0.7121, 1.28),
                 bbox_transform=axins.transData,
                 borderpad=0,
                 borderaxespad=0.2,
                 handletextpad=0.1,
                 ncol=1,
                 labelspacing=2.6,
                 columnspacing=0.7,
                 handlelength=1.38)

    # for ax in (ax1, ax2, ax3, axins):
    #     ax.tick_params(labelsize=plt_args['fontsize'],
    #                    pad=plt_args['labelpad'])

    # Lettering
    a_lett = ax1.annotate('(a)',
                          xy=(0., 1.),
                          xytext=(8, -8),
                          xycoords='axes fraction',
                          textcoords='offset points',
                          bbox=dict(fc='none', ec='none'),
                          ha='center',
                          va='center')
    b_lett = ax2.annotate('(b)',
                          xy=(0., 1.),
                          xytext=(8, -8),
                          xycoords='axes fraction',
                          textcoords='offset points',
                          bbox=dict(fc='none', ec='none'),
                          ha='center',
                          va='center')

    # axins.annotate('nav/nB',
    #                xy=(0.706, 1.16),
    #                xycoords='data',
    #                fontsize=plt_args['fontsize'])

    #ticklab = ax2.xaxis.get_ticklabels()[0]
    #trans=ticklab.get_transform()
    #ax2.xaxis.set_label_coords(0.0, 0.24, transform=trans)

    # plt.subplots_adjust(left=0.09,
    #                     right=0.99,
    #                     top=0.989,
    #                     bottom=0.1,
    #                     wspace=0.2,
    #                     hspace=0.35)
    plt.tight_layout()
    plt.savefig(f'{filename}.pdf', facecolor='none', edgecolor='none', dpi=400)


if __name__ == '__main__':
    plot_data()
