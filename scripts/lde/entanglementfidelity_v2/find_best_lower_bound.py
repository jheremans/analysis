import os, sys, time
import pickle
import pprint
import numpy as np
from scipy.special import erf

from matplotlib import pyplot as plt
from matplotlib import rcParams

from analysis.lib.fitting import fit
from analysis import config

import MLE

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'
config.datadir = r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde'
savedir = os.path.join(config.outputdir,
                time.strftime('%Y%m%d')+'-ldefidelity')

rcParams['mathtext.default'] = 'regular'
rcParams['legend.numpoints'] = 1

### set stuff

fidelity_srcdir = '20121128-ldefidelity'

state = 'psi2'
ZZidx = 0
XXidx = 1
XmXidx = 2

### load the fidelity data (best guess, i.e., sqrt-term=0)
fid = {}
f = np.load(os.path.join(config.outputdir,fidelity_srcdir,
    'fidelities_dtmin0.npz'))
for k in f.keys():
    fid[k] = f[k]
f.close()

corr = {}
f = np.load(os.path.join(config.outputdir,fidelity_srcdir,
    'correlations_dtmin0.npz'))
for k in f.keys():
    corr[k] = f[k]
f.close()

# get the desired raw correlations
_c = corr['raw'+state+'correlations']

for i,win1 in enumerate(corr['win1vals']):
    win1idx = i
    win1 = corr['win1vals'][win1idx]
    win2idx = argmin(abs(corr['win2vals']-win1))
    win2 = corr['win2vals'][win2idx]

    c = _c[:,win1idx,win2idx,...]
    dtvals = corr['dtvals']
    edges = np.append(0,dtvals)
    meandtvals = edges[:-1] + (edges[1:]-edges[:-1])/2.

    # do an MLE for each dt bin
    F = np.array([])
    u_F = np.array([])

    FLB = np.array([])
    u_FLB = np.array([])

    # max probability for not having entanglement
    PNE = np.array([])
    sigmas = np.array([])
    N = np.array([])

    for i,dt in enumerate(dtvals):
        emily = MLE.Emily()
        emily.prob_pts = 101
        emily.state = state
        emily.N_ZZ = c[i,...,ZZidx,:]
        emily.N_XX = c[i,...,XXidx,:]
        emily.N_XmX = c[i,...,XmXidx,:]
        emily.fac_exact = 1
        emily.likelihood()
        emily.analysis(savefn='MLE_win1%d_win2%d_dt%d' \
                % (win1,win2,dt), save=False)
        
        F = np.append(F, emily.F)
        u_F = np.append(u_F, emily.u_F)
        
        FLB = np.append(FLB, emily.FLB)
        u_FLB = np.append(u_FLB, emily.u_FLB)
        
        PNE = np.append(PNE, 0.5 * abs(1. - \
                erf((emily.FLB-0.5)/np.sqrt(2.)/emily.u_FLB)))
        sigmas = np.append(sigmas, (emily.FLB-0.5)/emily.u_FLB)
        N = np.append(N, emily.N_ZZ.sum()+emily.N_XX.sum()+emily.N_XmX.sum())
        
        plt.close('all')

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)

    ax.plot(np.append(0,dtvals), np.append(F[0],F), 'k--',
            drawstyle='steps', label='best guess')

    #ax.errorbar(meandtvals, F, fmt='o', mec='k', mfc='w',
    #        yerr=u_F, ecolor='k')

    ax.plot(np.append(0,dtvals), np.append(FLB[0],FLB), 'r',
            drawstyle='steps', label='Lower bd')
    ax.errorbar(meandtvals, FLB, fmt='o', mec='r', mfc='w',
            yerr=u_FLB, ecolor='r')

    for i,p in enumerate(PNE):
        ax.text(meandtvals[i], 0.99, '$P(NE) \leq$ %.2f%% ($%.1f\sigma$)' % (p*100.,sigmas[i]), ha='center',
                va='top', color='r', rotation='vertical', size='small')
        ax.text(meandtvals[i], 0.21, '$N = %d$' % (N[i]), ha='center', va='bottom',
                color='k', rotation='vertical', size='small')
        ax.text(meandtvals[i], 0.31, '$(%.3f \pm %.3f)$%%' %(FLB[i], u_FLB[i]),
                ha='center', va='bottom',
                color='r', rotation='vertical', size='small')
        ax.text(meandtvals[i], 0.46, '$(%.3f \pm %.3f)$%%' %(F[i], u_F[i]),
                ha='center', va='bottom',
                color='k', rotation='vertical', size='small')


    ax.set_ylim(0.2,1.0)
    ax.set_xlabel('dt')
    ax.set_ylabel('Fidelity')
    ax.legend(loc=4)
    ax.set_title('Psi2, win1 = %d, win2 = %d' % (win1, win2))

    fig.savefig(os.path.join(savedir, 'Numbers_psi2_F_vs_dtslice_win1%d_win2%d.png' \
            % (win1,win2)))


    
# now the psi1 fidelity for a first window of 150bins
state = 'psi1'
_c = corr['raw'+state+'correlations']

for i,win2 in enumerate(corr['win2vals']):
    win1 = 150
    win1idx = argmin(abs(corr['win1vals']-win1))
    win2idx = i
    win2 = corr['win2vals'][win2idx]

    c = _c[:,win1idx,win2idx,...]
    dtvals = corr['dtvals']
    edges = np.append(0,dtvals)
    meandtvals = edges[:-1] + (edges[1:]-edges[:-1])/2.

    # do an MLE for each dt bin
    F = np.array([])
    u_F = np.array([])

    FLB = np.array([])
    u_FLB = np.array([])

    # max probability for not having entanglement
    PNE = np.array([])
    sigmas = np.array([])
    N = np.array([])

    for i,dt in enumerate(dtvals):
        emily = MLE.Emily()
        emily.prob_pts = 101
        emily.state = state
        emily.N_ZZ = c[i,...,ZZidx,:]
        emily.N_XX = c[i,...,XXidx,:]
        emily.N_XmX = c[i,...,XmXidx,:]
        emily.fac_exact = 1
        emily.likelihood()
        emily.analysis(savefn='MLE_win1%d_win2%d_dt%d' \
                % (win1,win2,dt), save=False)
        
        F = np.append(F, emily.F)
        u_F = np.append(u_F, emily.u_F)
        
        FLB = np.append(FLB, emily.FLB)
        u_FLB = np.append(u_FLB, emily.u_FLB)
        
        PNE = np.append(PNE, 0.5 * abs(1. - \
                erf((emily.FLB-0.5)/np.sqrt(2.)/emily.u_FLB)))
        sigmas = np.append(sigmas, (emily.FLB-0.5)/emily.u_FLB)
        N = np.append(N, emily.N_ZZ.sum()+emily.N_XX.sum()+emily.N_XmX.sum())
        
        plt.close('all')

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)

    ax.plot(np.append(0,dtvals), np.append(F[0],F), 'k--',
            drawstyle='steps', label='best guess')

    #ax.errorbar(meandtvals, F, fmt='o', mec='k', mfc='w',
    #        yerr=u_F, ecolor='k')

    ax.plot(np.append(0,dtvals), np.append(FLB[0],FLB), 'r',
            drawstyle='steps', label='Lower bd')
    ax.errorbar(meandtvals, FLB, fmt='o', mec='r', mfc='w',
            yerr=u_FLB, ecolor='r')

    for i,p in enumerate(PNE):
        ax.text(meandtvals[i], 0.99, '$P(NE) \leq$ %.2f%% ($%.1f\sigma$)' % (p*100.,sigmas[i]), ha='center',
                va='top', color='r', rotation='vertical', size='small')
        ax.text(meandtvals[i], 0.21, '$N = %d$' % (N[i]), ha='center', va='bottom',
                color='k', rotation='vertical', size='small')
        ax.text(meandtvals[i], 0.31, '$(%.3f \pm %.3f)$%%' %(FLB[i], u_FLB[i]),
                ha='center', va='bottom',
                color='r', rotation='vertical', size='small')
        ax.text(meandtvals[i], 0.46, '$(%.3f \pm %.3f)$%%' %(F[i], u_F[i]),
                ha='center', va='bottom',
                color='k', rotation='vertical', size='small')


    ax.set_ylim(0.2,1.0)
    ax.set_xlabel('dt')
    ax.set_ylabel('Fidelity')
    ax.legend(loc=4)
    ax.set_title('Psi1, win1 = %d, win2 = %d' % (win1, win2))

    fig.savefig(os.path.join(savedir, 'Numbers_psi1_F_vs_dtslice_win1%d_win2%d.png' \
            % (win1,win2)))










