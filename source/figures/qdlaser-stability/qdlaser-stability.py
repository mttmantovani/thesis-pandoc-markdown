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
    'height': 11.5/2.54,
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
    gs = fig.add_gridspec(2,1, hspace=0.3, height_ratios=[1,1.4])
    gs0 = gs[1].subgridspec(3,2, wspace=0.3, hspace=0, width_ratios=[1.7,1])
    gs1 = gs[0].subgridspec(1,2, wspace=0.3)

    ax1 = fig.add_subplot(gs0[:,0])  # nav
    ax2 = fig.add_subplot(gs0[0,1])  # pn1
    ax3 = fig.add_subplot(gs0[1,1])  # pn2
    ax4 = fig.add_subplot(gs0[2,1])  # pn3
    ax5 = fig.add_subplot(gs1[1])  # stability semiclassical
    ax6 = fig.add_subplot(gs1[0]) # stability numerical


    # Nav
    des = np.unique(dataDict[0][:,0])
    gs = np.unique(dataDict[0][:,1])
    nav = dataDict[0][:,2].reshape(len(gs), len(des))
    im = ax1.imshow(nav.T, extent=[des.min(), des.max(), gs.min(), gs.max()], 
               aspect='auto', cmap='Spectral_r', origin='lower', vmin=0.)
    cb = fig.colorbar(im, ax=ax1)
    cb.ax.set_ylabel(r'$\bar{n}$', rotation=0)
    cb.ax.yaxis.set_label_coords(1.5, 50, transform=cb.ax.yaxis.get_ticklabels()[0].get_transform())

    ax1.plot([1], [0.13], linestyle='None', marker='^', markersize=7, markerfacecolor='red', markeredgewidth=0.9, markeredgecolor='k', zorder=3)
    ax1.plot([2], [0.13], linestyle='None', marker='*', markersize=9, markerfacecolor='chartreuse', markeredgewidth=0.9, markeredgecolor='k', zorder=3)
    ax1.plot([3], [0.29], linestyle='None', marker='o', markersize=7, markerfacecolor='yellow', markeredgewidth=0.9, markeredgecolor='k', zorder=3)

    # pn
    ax2.bar(dataDict[1][:,0], dataDict[1][:,1], color='dodgerblue')
    ax3.bar(dataDict[2][:,0], dataDict[2][:,1], color='dodgerblue')
    ax4.bar(dataDict[3][:,0], dataDict[3][:,1], color='dodgerblue')

    ax2.plot([90], [0.08], linestyle='None', marker='^', markersize=7, markerfacecolor='red', markeredgewidth=1, markeredgecolor='k')
    ax3.plot([90], [0.08], linestyle='None', marker='*', markersize=9, markerfacecolor='chartreuse', markeredgewidth=1, markeredgecolor='k')
    ax4.plot([90], [0.08], linestyle='None', marker='o', markersize=7, markerfacecolor='yellow', markeredgewidth=1, markeredgecolor='k')


    # Stability numerical
    Z1 = dataDict[0][:,4].reshape(len(gs), len(des))
    cmap = cm.get_cmap('Greens', np.max(Z1) + 1)

    im = ax5.imshow(Z1.T, cmap=cmap, extent=[des.min(), des.max(), gs.min(), gs.max()],
                       origin='lower', aspect='auto', zorder=0)

    contours = ax5.contour(Z1.T, [0.5, 1.5, 2.5, 3.5, 4.5, 5.5],
                       linewidths=1, colors='k',
                       extent=[des.min(), des.max(), gs.min(), gs.max()],
                       origin='lower', aspect='auto', zorder=1)

    ax5.plot([1], [0.13], linestyle='None', marker='^', markersize=7, markerfacecolor='red', markeredgewidth=0.9, markeredgecolor='k', zorder=3)
    ax5.plot([2], [0.13], linestyle='None', marker='*', markersize=9, markerfacecolor='chartreuse', markeredgewidth=0.9, markeredgecolor='k', zorder=3)
    ax5.plot([3], [0.29], linestyle='None', marker='o', markersize=7, markerfacecolor='yellow', markeredgewidth=0.9, markeredgecolor='k', zorder=3)

    # Stability semiclassical
    des2 = np.unique(dataDict[4][:,0])
    gs2 = np.unique(dataDict[4][:,1])
    Z2 = dataDict[4][:,2].reshape(len(gs2), len(des2))

    im2 = ax6.imshow(Z2.T, cmap=cmap, extent=[des2.min(), des2.max(), gs2.min(), gs2.max()],
                       origin='lower', aspect='auto', zorder=0)

    contours2 = ax6.contour(Z2.T, [0.5, 1.5, 2.5, 3.5, 4.5, 5.5],
                       linewidths=1, colors='k',
                       extent=[des2.min(), des2.max(), gs2.min(), gs2.max()],
                       origin='lower', aspect='auto', zorder=1)
    
    ax6.plot([1], [0.13], linestyle='None', marker='^', markersize=7, markerfacecolor='red', markeredgewidth=0.9, markeredgecolor='k', zorder=3)
    ax6.plot([2], [0.13], linestyle='None', marker='*', markersize=9, markerfacecolor='chartreuse', markeredgewidth=0.9, markeredgecolor='k', zorder=3)
    ax6.plot([3], [0.29], linestyle='None', marker='o', markersize=7, markerfacecolor='yellow', markeredgewidth=0.9, markeredgecolor='k', zorder=3)


#   # Adjustments
    for ax in (ax1, ax5, ax6):
        ax.set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
        ax.set_yticklabels([0, '', 0.1, '', 0.2, '', 0.3])
        ax.set_xlabel(r'$\Delta\epsilon\ [\omega_0]$')
        ax.set_ylabel(r'$g\ [\omega_0]$')
      #  ax.yaxis.set_label_coords(0., 0.15, transform=ax.yaxis.get_ticklabels()[0].get_transform())

    for ax in (ax2, ax3, ax4):
        ax.set_ylim(bottom=0, top=0.1)
        ax.set_yticks([])
        ax.set_xlim(0,100)

    for ax in (ax2,ax3):
        ax.set_xticks([])
    
    ax3.set_ylabel('Fock distribution ($p_n$)')
    ax4.set_xticks([0, 50, 100])
    ax4.set_xticklabels([0, '', 100])
    ax4.set_xlabel('$n$', va='top')
    ax4.xaxis.set_label_coords(50, 0, transform=ax4.xaxis.get_ticklabels()[0].get_transform())

    # Letters and annotations
    ax1.annotate('(c)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')
    ax2.annotate('(d)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')
    ax3.annotate('(e)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')
    ax4.annotate('(f)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')
    # ax5.annotate('(e)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')
    ax5.annotate('(b) numerical', xy=(1,0), xytext=(-5,5), xycoords='axes fraction', textcoords='offset points', ha='right', va='bottom')
    # ax6.annotate('(f)', xy=(0,1), xytext=(2,-2), xycoords='axes fraction', textcoords='offset points', ha='left', va='top')
    ax6.annotate('(a) semiclassical', xy=(1,0), xytext=(-5,5), xycoords='axes fraction', textcoords='offset points', ha='right', va='bottom')


    # Numbers
    stable_pos = {'1': [(1, 0.08), (1.3, 0.2), (2.5,0.1),
                        (3.5, 0.15)],
                  '2': [(1.65, 0.136), (2.25, 0.225), (3.5, 0.225)],
                  '3': [(2.75, 0.27), (3.78, 0.287 )],
    }

    for ax in (ax5, ax6):
        for k in stable_pos:
            for v in stable_pos[k]:
                ax.annotate(k, xy=v, xycoords='data', va='center', ha='center', fontstyle='italic')
        ax.annotate('3', xy=(2.5,0.17), xytext=(2.5,0.225 ), xycoords='data', va='center', ha='center', fontstyle='italic', arrowprops=dict(arrowstyle='->', lw=1, color='w'),
        bbox=dict(pad=-1, facecolor='none', edgecolor='none'))


    ax5.annotate('4', xy=(3.28,0.295), xytext=(3.3,0.23 ), xycoords='data', va='center', ha='center', fontstyle='italic', arrowprops=dict(arrowstyle='->', lw=1, color='w'),
        bbox=dict(pad=-1, facecolor='none', edgecolor='none'))
    ax6.annotate('4', xy=(3.4,0.295), xytext=(3.3,0.25 ), xycoords='data', va='center', ha='center', fontstyle='italic', arrowprops=dict(arrowstyle='->', lw=1, color='w'),
        bbox=dict(pad=-1, facecolor='none', edgecolor='none'))

    plt.subplots_adjust(left=0.08,bottom=0.1, right=0.98, top=0.98)
    plt.savefig(f'{filename}.pdf', dpi=400)


if __name__ == '__main__':
    plot_data()
