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
    'height': 5/2.54,
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
    fig = plt.figure(figsize=(defaults['width'], defaults['height']))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.15,1,1])
    ax1 = fig.add_subplot(gs[0,0])  # n vs g
    ax2 = fig.add_subplot(gs[0,1])  # n vs (g, deps)
    ax3 = fig.add_subplot(gs[0,2])  # #stable sols


    # Ax1
    ax1.plot(dataDict[0][:,0], dataDict[0][:,1], ls='-', label=r'$\gamma_\mathrm{sr}/\omega_0 = 0$')
    ax1.plot(dataDict[2][:,0], dataDict[2][:,1], ls='--', color='red', label=r'$\gamma_\mathrm{sr}/\omega_0= 10^{-3}$')
    ax1.plot(dataDict[1][:,0], dataDict[1][:,1], ls=':', color='green', label=r'$\gamma_\mathrm{sr}/\omega_0 = 10^{-2}$')

    #Ax 2 and 3
    des = np.unique(dataDict[3][:,0])
    gs = np.unique(dataDict[3][:,1])

    ph = dataDict[3][:,2].reshape(len(des), len(gs))
    peaks = dataDict[3][:, 4].reshape(len(des), len(gs))

    im1 = ax2.imshow(ph.T, extent=[des.min(), des.max(), gs.min(), gs.max()], 
            origin='lower', aspect='auto', interpolation='lanczos', 
            cmap='Spectral_r', vmax=35)
    cbar = fig.colorbar(im1, ax=ax2)

    im2 = ax3.imshow(peaks.T, cmap=cm.get_cmap('Greens', np.max(peaks) + 1),   extent=[des.min(), des.max(), gs.min(), gs.max()], origin='lower', aspect='auto', zorder=0)

    #cont2 = ax2.contour(peaks.T, [0.5, 1.5, 2.5], linewidths=1, colors='k',
    #                extent=[des.min(), des.max(), gs.min(), gs.max()],
    #                origin='lower', aspect='auto', zorder=1)



    # Adjustments
    ax1.set_xlim(0, 0.15)
    ax1.set_ylim(bottom=9)
    ax1.set_xlabel(r'$g\ [\omega_0]$')
    ax1.set_ylabel(r'$\bar{n}$')

    ax2.set_xlabel(r'$\Delta\epsilon\ [\omega_0]$')
    ax2.set_ylabel(r'$g\ [\omega_0]$')
    ax2.set_xticks(list(range(0,5)))

    ax3.set_xlabel(r'$\Delta\epsilon\ [\omega_0]$')
    ax3.set_ylabel(r'$g\ [\omega_0]$')
    ax3.set_xticks(list(range(0,5)))

    # Legend
    ax1.legend(borderaxespad=0.2, loc='lower right', frameon=False)

    # Letters and annotations
    ax1.annotate('(a)', xy=(0.02, 0.92), xycoords='axes fraction')
    ax2.annotate('(b)', xy=(0.02, 0.92), xycoords='axes fraction')
    ax2.annotate(r'$\bar{n}$', xy=(1.24, 0.92), xycoords='axes fraction')
    ax3.annotate('(c)', xy=(0.02, 0.92), xycoords='axes fraction')

    ax3.annotate('1', xy=(1,0.1), xycoords='data', va='center', ha='center', fontstyle='italic', )
    ax3.annotate('1', xy=(2.5,0.15), xycoords='data', va='center', ha='center', fontstyle='italic', )
    ax3.annotate('1', xy=(3.5,0.25), xycoords='data', va='center', ha='center', fontstyle='italic', )
    ax3.annotate('1', xy=(1.9,0.25), xycoords='data', va='center', ha='center', fontstyle='italic', )
    ax3.annotate('2', xy=(1.9,0.19), xycoords='data', va='center', ha='center', fontstyle='italic', )
    ax3.annotate('2', xy=(2.6,0.28), xycoords='data', va='center', ha='center', fontstyle='italic', )
    ax3.annotate('3', xy=(2.5,0.23), xytext=(3.1,0.18), xycoords='data', va='center', ha='center', fontstyle='italic', arrowprops=dict(arrowstyle='->', lw=0.4),
        bbox=dict(pad=-1, facecolor='none', edgecolor='none'))


    plt.subplots_adjust(left=0.07,bottom=0.23, right=0.99, top=0.97, wspace=0.4, hspace=0.10 )
    plt.savefig(f'{filename}.pdf', dpi=400)





if __name__ == '__main__':
    plot_data()
