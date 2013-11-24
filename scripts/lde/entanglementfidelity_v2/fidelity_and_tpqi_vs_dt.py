import os, sys, time
import pickle
import pprint
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rcParams

from analysis.lib.fitting import fit
from analysis import config

import MLE

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'
config.datadir = r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde'
savedir = os.path.join(config.outputdir,
                time.strftime('%Y%m%d')+'-ldefidelity_manydts')

rcParams['mathtext.default'] = 'regular'
rcParams['legend.numpoints'] = 1

### set stuff

# what to do
do_analysis = True
do_howmuchdata = True

# sources
tpqi_srcdir = '20121214-tpqi_manydts'
fidelity_srcdir = '20121205-ldefidelity_manydts'
ap_srcdir = '20121205-ldefidelity_manydts'

### tools
def idx(arr, val):
    return argmin(abs(arr-val))

### load the tpqi visibility data
tpqi = {}
f = np.load(os.path.join(config.outputdir,tpqi_srcdir,
    'cumulative_tpqi_ch0640_win1000.npz'))
for k in f.keys():
    tpqi[k] = f[k]
f.close()

### load the correlations
corr = {}
f = np.load(os.path.join(config.outputdir,fidelity_srcdir,
    'correlations_dtmin0.npz'))
for k in f.keys():
    corr[k] = f[k]
f.close()

p_afterpulsing = np.loadtxt(os.path.join(config.outputdir, ap_srcdir, 'apratios.dat'))
u_p_afterpulsing = np.zeros(len(p_afterpulsing))

ZZidx = 0
XXidx = 1
XmXidx = 2

win1 = 150
win2s = [75, 150]
states = ['psi1', 'psi2']

for state,win2 in zip(states, win2s):
    # get the desired raw correlations
    _c = corr['raw'+state+'correlations']
    win1idx = argmin(abs(corr['win1vals']-win1))
    win2idx = argmin(abs(corr['win2vals']-win2))
    c = _c[:,win1idx,win2idx,...]
    dtvals = corr['dtvals']
    edges = np.append(0,dtvals)
    meandtvals = edges[:-1] + (edges[1:]-edges[:-1])/2.

    lde_ZZ_events = np.array([])
    lde_XX_events = np.array([])
    lde_XmX_events = np.array([])
    tpqi_events = np.array([])

    if do_howmuchdata:
        dts = len(dtvals)

        for i,dt in enumerate(dtvals):
            lde_ZZ_events = np.append(lde_ZZ_events, c[i,...,ZZidx,:].sum())
            lde_XX_events = np.append(lde_XX_events, c[i,...,XXidx,:].sum())           
            lde_XmX_events = np.append(lde_XmX_events, c[i,...,XmXidx,:].sum())
            tpqi_events = np.append(tpqi_events, 
                    tpqi['zprebinned'][dts-i-1:dts+i+1].sum())

        
        if not os.path.exists(savedir):
            os.makedirs(savedir)

        np.savez(os.path.join(savedir, state+'_dataamount_vs_dtslice_win1%d_win2%d' \
                % (win1, win2)),
                lde_ZZ_events=lde_ZZ_events,
                lde_XX_events=lde_XX_events,
                lde_XmX_events=lde_XmX_events,
                tpqi_events=tpqi_events)
            
    
    if do_analysis:
        # do an MLE for each dt bin
        F = np.array([])
        u_F = np.array([])
        FLB = np.array([])
        u_FLB = np.array([])
        for i,dt in enumerate(dtvals):
            emily = MLE.Emily()
            emily.prob_pts = 101
            emily.state = state
            emily.N_ZZ = c[i,...,ZZidx,:]
            emily.N_XX = c[i,...,XXidx,:]
            emily.N_XmX = c[i,...,XmXidx,:]
            emily.fac_exact = 1
            emily.likelihood()
            emily.analysis(save=False)
            F = np.append(F, emily.F)
            u_F = np.append(u_F, emily.u_F)
            FLB = np.append(FLB, emily.FLB)
            u_FLB = np.append(u_FLB, emily.u_FLB)
            plt.close('all')

        ### the psi2 figure
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # expected value from the TPQI visibility
        ax.plot(tpqi['binedges'], np.append(tpqi['fidelity'][0], tpqi['fidelity']),
                'b', drawstyle='steps', label='TPQI')
        ax.errorbar(tpqi['meandts'], tpqi['fidelity'], fmt='o', mec='b',
                mfc='w', yerr=tpqi['u_fidelity'], ecolor='b')

        # include afterpulsing
        if state == 'psi1':
            F_tpqi_ap = p_afterpulsing * 0.25 + \
                    (1-p_afterpulsing) * tpqi['fidelity']
            u_F_tpqi_ap = ( ((0.25-tpqi['fidelity'])*u_p_afterpulsing)**2 + \
                    ((1-p_afterpulsing)*tpqi['u_fidelity'])**2 )**0.5
            ax.plot(tpqi['binedges'], np.append(F_tpqi_ap[0], F_tpqi_ap),
                    'r', drawstyle='steps', label='TPQI+AP')
            ax.errorbar(tpqi['meandts'], F_tpqi_ap, fmt='o', mec='r',
                    mfc='w', yerr=u_F_tpqi_ap, ecolor='r')
        else:
            F_tpqi_ap = np.zeros(len(tpqi['fidelity']))
            u_F_tpqi_ap = np.zeros(len(tpqi['fidelity']))

        # fidelity best guess from spin-spin correlations
        ax.plot(np.append(0,dtvals), np.append(F[0],F), 'k',
                drawstyle='steps', label='LDE best guess')
        ax.errorbar(meandtvals, F, fmt='o', mec='k', mfc='w',
                yerr=u_F, ecolor='k')

        # the lower bd
        ax.plot(np.append(0,dtvals), np.append(FLB[0],FLB), 'k--',
                drawstyle='steps', label='LDE lower bd')

        ax.legend()
        ax.set_xlabel('$dt_{max}$')
        ax.set_ylabel('Fidelity')
        ax.set_title(state)

        if not os.path.exists(savedir):
            os.makedirs(savedir)

        fig.savefig(os.path.join(savedir, state+'_F_vs_dtslice_win1%d_win2%d.png' \
                % (win1,win2)))

        np.savez(os.path.join(savedir, state+'_F_vs_dtslice_win1%d_win2%d' \
                % (win1, win2)), 
                tpqi_binedges = tpqi['binedges'],
                tpqi_fidelity = tpqi['fidelity'],
                tpqi_dt = tpqi['meandts'],
                u_tpqi_fidelity = tpqi['u_fidelity'],
                tpqi_fidelity_with_ap = F_tpqi_ap,
                u_tpqi_fidelity_with_ap = u_F_tpqi_ap,
                dt = meandtvals,
                F = F,
                u_F = u_F,
                FLB = FLB,
                u_FLB = u_FLB)



    










