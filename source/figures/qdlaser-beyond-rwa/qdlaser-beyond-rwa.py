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
    gs = fig.add_gridspec(1,2, wspace=0.15, hspace=0.3, width_ratios=[1,1.3])

    gs0 = gs[0].subgridspec(2,1)
    gs1 = gs[1].subgridspec(2,1, height_ratios=[0.05, 1], hspace=0.05)

    ax1 = fig.add_subplot(gs0[0])  # nav RWA-noRWA
    ax2 = fig.add_subplot(gs0[1])  # Fano RWA-noRWA
    ax3 = fig.add_subplot(gs1[1])  # pns vs g
    cax = fig.add_subplot(gs1[0])  # colorbar
    axins = inset_axes(
    ax3,
    loc='lower center',
    width="45%",
    height="30%",
    # bbox_to_anchor=(0.6, 0.2),
    borderpad=2,
    bbox_transform=ax3.transAxes)

    # Plot
    ax1.plot(dataDict[0][:,0], dataDict[0][:,1], ls='--', c='gray', 
             label='RWA', zorder=3)
    ax1.plot(dataDict[1][:,0], dataDict[1][:,1], ls='-', c='dodgerblue', 
             label='without RWA')
    ax2.plot(dataDict[0][:,0], dataDict[0][:,2], ls='--', c='gray', 
             label='RWA', zorder=3)
    ax2.plot(dataDict[1][:,0], dataDict[1][:,2], ls='-', c='red', 
             label='without RWA')


    gs = np.unique(dataDict[3][:,0])
    ns = np.unique(dataDict[3][:,1])
    pns = dataDict[3][:,2].reshape(len(ns), len(gs))
    im = ax3.imshow(pns, extent=[gs.min(), gs.max(), ns.min(), ns.max()], 
               aspect='auto', cmap='Blues', origin='lower', vmax=0.03)
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[0, 0.01, 0.02, 0.03])
    cax.xaxis.tick_top()
    cax.set_xlabel('$p_n$')
    cax.xaxis.set_label_position('top')
    cax.xaxis.set_label_coords(0.015, 1.5, transform=cax.xaxis.get_ticklabels()[0].get_transform())

    ax3.axvline(x=0.0559, ls=':', color='dimgray')
    ax3.plot(0.0559, 250, marker='v', markeredgecolor='k', markerfacecolor='lightskyblue')

    axins.fill_between(dataDict[2][:,0], dataDict[2][:,1], cmap='Blues')
    axins.plot(30,0.0065, marker='v', markeredgecolor='k', markerfacecolor='lightskyblue')

    # Adjustments
    ax1.set_xticklabels([])
    ax1.set_ylim(bottom=0)
    ax1.set_xlim(0, 0.2)
    ax1.set_xticks([0, 0.05, 0.1, 0.15, 0.2])
    ax1.set_ylabel(r'$\bar{n}$')

    ax2.set_xlim(0, 0.2)
    ax2.set_ylim(bottom=0)
    ax2.set_xlabel(r'$g\ [\omega_0]$')
    ax2.set_xticks([0, 0.05, 0.1, 0.15, 0.2])
    ax2.set_xticklabels([0, '', 0.1, '', 0.2])
    ax2.set_ylabel(r'$\mathcal{F}$')


    ax3.set_xlabel(r'$g\ [\omega_0]$')
    ax3.set_ylabel('$n$')
    ax3.set_xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
    ax3.set_xticklabels([0, '', 0.1, '', 0.2, '', 0.3])
    ax3.yaxis.set_label_coords(0, 250, transform=ax3.yaxis.get_ticklabels()[0].get_transform())

    axins.set_xlim(0,500)
    axins.set_xticks([0, 250, 500])
    axins.set_xticklabels([0, '', 500])
    axins.set_xlabel('$n$', color='dimgray')
    axins.xaxis.set_label_coords(250, 0, transform=axins.xaxis.get_ticklabels()[0].get_transform())  
    axins.set_ylim(bottom=0)  
    axins.set_yticks([])
    axins.set_ylabel('$p_n$', color='dimgray')
    axins.tick_params(colors='dimgray')
    for spine in axins.spines.values():
        spine.set_edgecolor('dimgray')



    # Legend
    ax1.legend(loc='lower right', frameon=False)
    ax2.legend(loc='upper right', frameon=False)

    # Letters and annotations
    ax1.annotate('(a)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')
    ax2.annotate('(b)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')
    ax3.annotate('(c)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')


    plt.subplots_adjust(left=0.08,bottom=0.11, right=0.98, top=0.92, wspace=0.4, hspace=0.08 )
    plt.savefig(f'{filename}.pdf', dpi=400, edgecolor='None')


if __name__ == '__main__':
    plot_data()
