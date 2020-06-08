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
    'width': TEXTWIDTH,
    'height': 10/2.54,
}

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

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
    gs = fig.add_gridspec(1,2, wspace=0.15, hspace=0.3, width_ratios=[1,1.3])

    gs0 = gs[0].subgridspec(2,1)
    gs1 = gs[1].subgridspec(2,1, height_ratios=[0.05, 1], hspace=0.05)

    ax1 = fig.add_subplot(gs0[0])  # nav RWA-noRWA
    ax2 = fig.add_subplot(gs0[1])  # Fano RWA-noRWA
    ax3 = fig.add_subplot(gs1[1])  # pns vs g
    cax = fig.add_subplot(gs1[0])  # colorbar

    # Plot
    gs = np.unique(dataDict[0][:,0])
    ns = np.unique(dataDict[0][:,1])
    pns = dataDict[0][:,2].reshape(len(ns), len(gs))
    im = ax3.imshow(pns, extent=[gs.min(), gs.max(), ns.min(), ns.max()], 
               aspect='auto', origin='lower', vmax=0.03)
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[0, 0.01, 0.02, 0.03])
    cax.xaxis.tick_top()
    cax.set_xlabel('$p_n$')
    cax.xaxis.set_label_position('top')
    cax.xaxis.set_label_coords(0.015, 1.5, transform=cax.xaxis.get_ticklabels()[0].get_transform())


    # ax1.plot(dataDict[0][:,0], dataDict[0][:,3], ls='-', color='red', label='analytical')
    # ax1.plot(dataDict[0][:,0], dataDict[0][:,5], ls='--', color='limegreen', label='semiclassical')
    # ax1.plot(dataDict[0][:,0], dataDict[0][:,1], marker='.',  markersize=8, markevery=20, ls='None', color='dodgerblue', label='numerical')
    # ax1.axvline(x=0.005, color='gray', ls=':')

    # ax2.plot(dataDict[0][:,0], dataDict[0][:,4], ls='-', color='red')
    # ax2.plot(dataDict[0][:,0], dataDict[0][:,2],  marker='.', markersize=8, markevery=20, ls='-', color='dodgerblue')
    # ax2.axvline(x=0.005, color='gray', ls=':')

    # for i, ax in enumerate([ax3, ax4, ax5]):
    #     ax.bar(dataDict[i+1][:,0], dataDict[i+1][:,1], color='dodgerblue', alpha=0.8, width=0.8)
    #     ax.plot(dataDict[i+1][:,2], color='red', lw=1)
    #     ax.set_xticks([0, 50, 100])
    #     ax.set_xticklabels([0, '', 100])
    #     ax.set_xlim(left=-0.4,right=100)
    #     ax.set_xlabel('$n$')
    #     ax.set_ylabel('$p_n$')
    #     ax.xaxis.set_label_coords(50, 0., transform=ax.xaxis.get_ticklabels()[0].get_transform())

    # Adjustments
    ax3.set_xlabel(r'$g\ [\omega_0]$')
    ax3.set_ylabel('$n$')
    ax3.set_xticks([0, 0.1, 0.2, 0.3])
    ax3.yaxis.set_label_coords(0, 250, transform=ax3.yaxis.get_ticklabels()[0].get_transform())
    # ax1.set_xlim(0, 0.05)
    # ax1.set_ylim(bottom=0, top=35)
    # ax1.set_yticks([5*i for i in range(8)])
    # ax1.set_yticklabels([0] + ['' for i in range(6)] + [35])
    # ax1.set_xlabel(r'$g\ [\omega_0]$')
    # ax1.set_ylabel(r'$\bar{n}$')

    # ax2.set_xlim(0, 0.05)
    # ax2.set_xlabel(r'$g\ [\omega_0]$')
    # ax2.set_ylabel(r'$\mathcal{F}$')

    # ax3.set_yticks([0, 0.2, 0.4])
    # ax3.set_yticklabels([0, '', 0.4])
    # ax3.yaxis.set_label_coords(0., 0.2, transform=ax3.yaxis.get_ticklabels()[0].get_transform())

    # ax4.set_yticks([0, 0.04, 0.08])
    # ax4.set_yticklabels([0, '', 0.08])
    # ax4.yaxis.set_label_coords(0., 0.04, transform=ax4.yaxis.get_ticklabels()[0].get_transform())

    # ax5.set_yticks([0, 0.04, 0.08])
    # ax5.set_yticklabels([0, '', 0.08])
    # ax5.yaxis.set_label_coords(0., 0.04, transform=ax5.yaxis.get_ticklabels()[0].get_transform())


#    # Legend
#     ax1.legend(borderaxespad=0.2, loc='lower right', frameon=False)
# #
      # Letters and annotations
#      ax1.annotate(r'(a)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')
#     ax2.annotate(r'(b)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')
    ax3.annotate('(c)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top', color='white')
#     ax4.annotate(r'$g = 0.006\omega_0\$ (d)', xy=(1,1), xytext=(-2,-2), xycoords='axes fraction', textcoords='offset points', ha='right', va='top')
#     ax5.annotate(r'$g = 0.05\omega_0\$ (e)', xy=(1,1), xytext=(-2,-2), xycoords='axes fraction', textcoords='offset points', ha='right', va='top')

    plt.subplots_adjust(left=0.07,bottom=0.15, right=0.97, top=0.92, wspace=0.4, hspace=0.10 )
    plt.savefig(f'{filename}.pdf', dpi=400, edgecolor='k')


if __name__ == '__main__':
    plot_data()
