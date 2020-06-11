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
    'height': 11/2.54,
}

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams.update(mpl_settings)
import numpy as np
from scipy.interpolate import interp1d

filename = os.path.splitext(argv[0])[0]

def arrowed_spines(fig, ax):

    xmin, xmax = ax.get_xlim() 
    ymin, ymax = ax.get_ylim()

    # removing the default axis on all sides:
    for side in ['bottom','right','top','left']:
       ax.spines[side].set_visible(False)

    # removing the axis ticks
    ax.set_xticks([]) # labels 
    ax.set_yticks([])
    ax.xaxis.set_ticks_position('none') # tick markers
    ax.yaxis.set_ticks_position('none')

    # get width and height of axes object to compute 
    # matching arrowhead length and width
    dps = fig.dpi_scale_trans.inverted()
    bbox = ax.get_window_extent().transformed(dps)
    width, height = bbox.width, bbox.height

    # manual arrowhead width and length
    hw = 1./25.*(ymax-ymin) 
    hl = 1./25.*(xmax-xmin)
    lw = 1. # axis line width
    ohg = 0.2 # arrow overhang

    # compute matching arrowhead length and width
    yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width 
    yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height

    # draw x and y axis
    ax.arrow(xmin, 0, xmax-xmin, 0., fc='k', ec='k', lw = lw, 
             head_width=hw, head_length=hl, overhang = ohg, 
             length_includes_head= True, clip_on = False) 

    ax.arrow(0, ymin, 0., ymax-ymin, fc='k', ec='k', lw = lw, 
             head_width=yhw, head_length=yhl, overhang = ohg, 
             length_includes_head= True, clip_on = False)

def telegraph_current(x):
    """ Dummy function for bistable current. """
    if x <= 0.2: 
        return 0.75
    if x > 0.2 and x <=0.4:
        return 0.2
    if x > 0.4 and x <=0.6:
        return 0.75
    if x > 0.6:
        return 0.2

def telegraph(tmin=0, tmax=0.8, c1=0.5, c2=telegraph_current, noise=0.03, npoints=250):
    """ Generate dummy data for telegraph plots. """
    x = np.linspace(tmin, tmax, npoints)
    y1 = c1 + np.random.normal(0, noise, x.shape)
    y2 = np.vectorize(c2)(x) + np.random.normal(0, noise, x.shape)

    return x, y1, y2

def tsm_rates(data_noise=None, data_tsm=None,  data_prob=None, npoints=150, kappa=1e-3):
    """ Two-state-model rates of switching. """
    noise = data_noise
    tsm = data_tsm
    probs = data_prob

    x   = np.linspace(tsm[:,0][0], tsm[:,0][-1], npoints)

    # Interpolating functions
    S0 = interp1d(noise[:,0], noise[:,2], kind='cubic')
    I1 = interp1d(tsm[:,0], tsm[:,1], kind='cubic')
    I2 = interp1d(tsm[:,0], tsm[:,2], kind='cubic')
    p1 = interp1d(probs[:,0], probs[:,1], kind='cubic')
    p2 = interp1d(probs[:,0], probs[:,2], kind='cubic')

    variance = p1(x)*p2(x)*(I1(x) - I2(x))**2

    sum_rates_tsm = 10.*variance/S0(x)
    sum_rates_tsm[-4:] = np.nan

    return x, kappa/sum_rates_tsm


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
    gs = fig.add_gridspec(2,1, hspace=0.15, height_ratios=[1.,1.7])

    gs0 = gs[1].subgridspec(2,2, wspace=0.2, height_ratios=[0.35,1])
    gs1 = gs[0].subgridspec(1,2)

    ax1 = fig.add_subplot(gs0[0,0]) # p_i
    ax2 = fig.add_subplot(gs0[1,0]) # current
    ax3 = fig.add_subplot(gs0[:,1]) # Rates
    ax4 = fig.add_subplot(gs1[0])
    ax5 = fig.add_subplot(gs1[1])


    # Plot
    ax1.plot(dataDict[0][:,0], dataDict[0][:,1], c='k', label=r'$\mathcal{P}_\mathrm{I}$')
    ax1.plot(dataDict[0][:,0], dataDict[0][:,2], c='r', label=r'$\mathcal{P}_\mathrm{II}$')


    ax2.plot(dataDict[1][:,0], 10*dataDict[1][:,1], label=r'$I\ [e\Gamma_R]$')
    ax2.plot(dataDict[2][:,0], 10*dataDict[2][:,1], c='k', ls='--', label=r'$I_\mathrm{I}$')
    ax2.plot(dataDict[2][:,0], 10*dataDict[2][:,2], c='r', ls='--', label=r'$I_\mathrm{II}$')
    ax2.plot(dataDict[2][:,0], 10*dataDict[2][:,3], c='orange', ls='--', label=r'$I_\mathrm{tsm}$')
    ax2.fill_between(dataDict[2][:,0], 10.*dataDict[2][:,1] + 0.03,  10.*dataDict[2][:,1] - 0.03, alpha=0.5, facecolor='silver', edgecolor=None, zorder=0)
    ax2.fill_between(dataDict[2][:,0], 10.*dataDict[2][:,2] + 0.03,  10.*dataDict[2][:,2] - 0.03, alpha=0.3, facecolor='red', edgecolor=None, zorder=0)
    ax2.plot([1], .285, linestyle='None', marker='^', 
        markersize=9, markerfacecolor='red', markeredgewidth=0.9, 
        markeredgecolor='k', zorder=5)
    ax2.plot([2], .21, linestyle='None', marker='*', 
        markersize=11, markerfacecolor='chartreuse', markeredgewidth=0.9, 
        markeredgecolor='k', zorder=5)

    de, sum_rates = tsm_rates(data_noise=dataDict[4], data_tsm=dataDict[2], data_prob=dataDict[0])

    
    ax3.plot(dataDict[3][:,0], -1e-3*dataDict[3][:,1], marker='o', ls='None', markerfacecolor='None', label=r'$\frac{\kappa}{|\lambda_1|}$')
    ax3.plot(de, sum_rates, color='orange', ls='--', label=r'$\frac{\kappa}{W_{\mathrm{I}\rightarrow\mathrm{II}} + W_{\mathrm{II}\rightarrow\mathrm{I}}}$')

    t, c1, c2 = telegraph()
    ax4.plot(t, c1, c='k'), ax4.axhline(y=0.5, ls='--', c='k')
    ax5.plot(t, c2, c='k'), ax5.axhline(y=0.2, ls='--', c='k')
    ax5.axhline(y=0.75, ls='--', c='k')
    ax4.plot([0.05], 0.95, linestyle='None', marker='^', 
        markersize=9, markerfacecolor='red', markeredgewidth=0.9, 
        markeredgecolor='k', zorder=5)
    ax5.plot([0.05], 0.95, linestyle='None', marker='*', 
        markersize=11, markerfacecolor='chartreuse', markeredgewidth=0.9, 
        markeredgecolor='k', zorder=5)



    # Adjustments
    ax1.set_xlim(left=0, right=2.4)
    ax1.set_xticks(np.linspace(0,2,5))
    ax1.set_xticklabels([])
    ax1.set_yticks(np.linspace(0,1,5))
    ax1.set_yticklabels([0,'','','',1])
    ax1.set_ylim(top=1.1)

    ax2.set_xlim(left=0, right=2.4)
    ax2.set_xticks(np.linspace(0,2,5))
    ax2.set_xticklabels([0, '', 1, '', 2])
    ax2.set_xlabel(r'$\Delta\epsilon\ [\omega_0]$')
    ax2.set_yticks([0, 0.1, 0.2, 0.3])

    ax3.set_xlim(1, 3)
    ax3.set_xticks(np.linspace(1,3,5))
    ax3.set_xlabel(r'$\Delta\epsilon\ [\omega_0]$')
    
    ax4.set_ylim(0,1), ax4.set_xlim(0,0.9)
    ax5.set_ylim(0,1), ax5.set_xlim(0,0.9)
    ax4.set_ylabel('$I(t)$')
    ax5.set_ylabel('$I(t)$')
    ax4.set_xlabel('$t$')
    ax5.set_xlabel('$t$')
    ax4.xaxis.set_label_coords(0.95,0.13)
    ax5.xaxis.set_label_coords(0.95,0.13)
    arrowed_spines(fig, ax4)
    arrowed_spines(fig, ax5)



    # Legends
    ax1.legend(loc='center left', frameon=False)
    ax2.legend(loc='best', borderpad=2, frameon=False)
    ax3.legend(loc='upper left', frameon=False)

    # Letters and annotations
    ax1.annotate('(c)', xy=(0,1), xytext=(-20,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='top')
    ax2.annotate('(d)', xy=(0,1), xytext=(-20,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='top')
    ax3.annotate('(e)', xy=(0,1), xytext=(-10,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='top')
    ax4.annotate('(a)', xy=(0,1), xytext=(-20,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='top')
    ax5.annotate('(b)', xy=(0,1), xytext=(-10,0), xycoords='axes fraction', textcoords='offset points', ha='center', va='top')




    ax4.arrow(0.45,0.5,0,-0.1,color='red', head_width=0.02, zorder=5)
    ax4.arrow(0.45,0.5,0,0.1,color='red', head_width=0.02, zorder=5)
    ax4.annotate(r'$\Delta I$', xy=(0.45,0.6), xytext=(0,5), xycoords='data', textcoords='offset points', ha='center', va='bottom', color='red')

    ax5.arrow(0.3,0.2,0,-0.06,color='red', head_width=0.02, zorder=5)
    ax5.arrow(0.3,0.2,0,0.06,color='red', head_width=0.02, zorder=5)
    ax5.annotate(r'$\Delta I_\mathrm{I}$', xy=(0.3,0.26), xytext=(0,5), xycoords='data', textcoords='offset points', ha='center', va='bottom', color='red')
    ax5.arrow(0.5,0.75,0,-0.06,color='red', head_width=0.02, zorder=5)
    ax5.arrow(0.5,0.75,0,0.06,color='red', head_width=0.02, zorder=5)
    ax5.annotate(r'$\Delta I_\mathrm{II}$', xy=(0.5,0.69), xytext=(0,-5), xycoords='data', textcoords='offset points', ha='center', va='top', color='red') 
    ax5.arrow(0.85,0.475,0,-0.475/2,color='red', head_width=0.02, zorder=5)
    ax5.arrow(0.85,0.475,0,0.475/2,color='red', head_width=0.02, zorder=5)
    ax5.annotate(r'$I_\mathrm{II} - I_\mathrm{I}$', xy=(0.85,0.475), xycoords='data', ha='center', va='center', color='red', backgroundcolor='w', zorder=6)

    plt.subplots_adjust(left=0.06,bottom=0.10, right=0.98, top=0.99, wspace=0.4, hspace=0.08 )
    plt.savefig(f'{filename}.pdf', dpi=400)


if __name__ == '__main__':
    plot_data()
