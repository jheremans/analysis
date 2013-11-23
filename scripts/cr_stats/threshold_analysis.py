import os, sys
import string
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import ssro
from analysis.lib.m2.ssro import sequence
reload(sequence)


folder = None
timestamp =  None#'111857'# None#'095551'#'180057'# None


if folder == None:
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('AdwinSSRO')

a = ssro.SSROAnalysis(folder)
b = sequence.SequenceAnalysis(folder)

runs = 1
sweep_probes = 1 # 2 if also analyze_probe, 1 if only preselect.

debug = False
if debug:
    runs =2 

ths = [ [ [ ] for i in range (runs)] for s in range (sweep_probes) ] 
taus = [ [ [ ] for i in range (runs)] for s in range (sweep_probes) ] 
u_taus =[ [ [ ] for i in range (runs)] for s in range (sweep_probes) ] 
means = [ [ [ ] for i in range (runs)] for s in range (sweep_probes) ] 
percentage_passes = [ [ [ ] for i in range (runs)] for s in range (sweep_probes) ] 


favo_th = 15
favo_mean = []
favo_cr_pass = []
favo_k = []
favo_u_k = []



# analyze the data

for i,g in enumerate(a.g.items()):    
    gn = g[0]
    if 'instrument_settings' in gn:
        continue
    if i>7:
        continue

    a.get_run(gn)
    pre = int(string.split(gn, '_')[2])
    pro = int(string.split(gn, '_')[4])
    run = int(string.split(gn, '_')[6])
    analyze_probe = str(string.split(gn, '_')[8])

    if debug:   
        if run > 2:
            continue

    if analyze_probe == 'True':
        th = pro
        probe = 1
    else:
        th = pre
        probe = 0

    _t, _c = a.readout_relaxation(a.ro_time, a.ro_counts, a.reps, a.binsize, name=gn,
            plot=False, ret=True)

    idx0 = np.argmax(_c)
    idx1 = -1 
    t,c = _t[idx0:idx1]/1e3, _c[idx0:idx1]

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot(t,c, 'o')

    res = fit.fit1d(t,c, common.fit_exp_decay_with_offset, 0, c[0], 10, 
            ret=True, do_print=False, fixed=[])
    
    if res != False:
        plot.plot_fit1d(res, t, ax=ax, plot_data=False)
        ths[probe][run-1].append(th)
        taus[probe][run-1].append(res['params_dict']['tau'])
        u_taus[probe][run-1].append(res['error_dict']['tau'])

    ax.set_xlabel('t ($\mu$s)')
    ax.set_ylabel('cts (Hz)')
    ax.set_title(a.default_plot_title + ', ' + gn)
    fig.savefig(os.path.join(folder, 'sp_relaxation_Th_pre={}_Th_pro={}_run_{}_sweep_{}'.format(pre,pro, run, th)))

    b.get_cr_results(gn, plot=True)

    hist,mean,var = b.get_mean_cr_cts()
    means[probe][run-1].append(mean[0])

    stats = b.adwingrp(gn)['statistics'].value
    fail = stats[2]
    percentage_pass = 5000./(5000. + fail) * 100 #5000 is the number of succesfull measurements.
    percentage_passes[probe][run-1].append(percentage_pass)

    plt.close('all')


# make plots

for r in np.arange(runs):
    #print 80*'='
    #print 'run number {}'.format(r+1)
    #print 80*'='

    for s in np.arange(sweep_probes):
        if s == 1:
            p = 'probe'
        else:
            p = 'preselect'

        #print 80*'='
        #print p+' threshold'
        #print 80*'='

        sortidxs = np.argsort(ths[s][r])
        ths[s][r] = np.array(ths[s][r])[sortidxs][:]
        taus[s][r] = np.array(taus[s][r])[sortidxs][:]
        u_taus[s][r] = np.array(u_taus[s][r])[sortidxs][:]
        means[s][r] = np.array(means[s][r])[sortidxs][:]

        percentage_passes[s][r] = np.array(percentage_passes[s][r])[sortidxs][:]

        ks = 1./taus[s][r] * 1e3
        u_ks = 1./taus[s][r]**2 * u_taus[s][r] * 1e3

        fig, ([ax1,ax2,ax3]) = plt.subplots(3,1, figsize = (8,10))
        #ax1 = plt.subplot2grid((2,2),(0,0))
        ax1.errorbar(ths[s][r], ks, yerr=u_ks, fmt='o', capsize=0, elinewidth=2)

        ax1.set_xlabel('{} threshold'.format(p))
        ax1.set_ylabel('relaxation rate (kHz)')
        ax1.set_title(a.timestamp + '\n{} threshold for run {}, relaxation rate'.format(p,r+1))

        fig.savefig(os.path.join(folder, 'Relaxation_rate_th_{}_run_{}'.format(p,r+1)))


        #fig = a.default_fig(figsize=(6,4))
        #ax2 = plt.subplot2grid((2,2),(0,1))#a.default_ax(fig)
        ax2.errorbar(ths[s][r], means[s][r], yerr=means[s][r]/np.sqrt(5000), fmt='o', capsize=0, elinewidth=2)
        #ax2.plot(ths[s][r], means[s][r],'o')
        ax2.set_xlabel('{} threshold'.format(p))
        ax2.set_ylabel('mean CR counts after sequence')
        ax2.set_title('\n' +'\n{} threshold for run {}, mean CR cts'.format(p,r+1))


        fig.savefig(
            os.path.join(folder, 'post-CR_sum_vs_sweepparam_th_{}_run_{}.png'.format(p,r+1)),
            format='png')

        #fig = a.default_fig(figsize=(6,4))
        #ax3 = plt.subplot2grid((2,2),(1,1),colspan =2)#a.default_ax(fig)
        ax3.plot(ths[s][r], percentage_passes[s][r],'o')
        ax3.set_xlabel('{} threshold'.format(p))
        ax3.set_ylabel('percentage CR passes')
        ax3.set_title('\n'+ '\n{} th for run {}, CR pass prob'.format(p,r+1))

        fig.savefig(
            os.path.join(folder, 'percentage_CR_pass_vs_sweepparam_th_{}_run_{}.png'.format(p,r+1)),
            format='png')

        #gather info about my favourite threshold combination(which occurs during the probe sweep):
        if runs > 2 and sweep_probes > 1:
            if s == 1:
                th_index = list(ths[s][r]).index(favo_th)
                favo_mean.append(means[s][r][th_index])
                favo_cr_pass.append(percentage_passes[s][r][th_index])
                favo_k.append(1/ taus[s][r][th_index] *1e3)
                favo_u_k.append(1./taus[s][r][th_index]**2 * u_taus[s][r][th_index] * 1e3)

# plot info for my favourite threshold

if runs > 2 and sweep_probes > 1:
    fig, ([ax1,ax2,ax3]) = plt.subplots(3,1, figsize = (8,10))
    ax1.errorbar(np.arange(runs)+1, favo_k, yerr=favo_u_k, fmt='o', capsize=0, elinewidth=2)

    ax1.set_xlabel('run number')
    ax1.set_ylabel('relaxation rate (kHz)')
    ax1.set_title(a.timestamp + '\n Relaxation rate vs. run for th:{}/{}'.format(30,favo_th))

    fig.savefig(os.path.join(folder, 'Relaxation_rate_vs_run_for_ths_{}_{}'.format(30,favo_th)))

    ax2.plot(np.arange(runs)+1, favo_mean,'o')
    ax2.set_xlabel('run number')
    ax2.set_ylabel('mean CR counts after sequence')
    ax2.set_title('\n' + 'Mean post CR cts vs. run for th:{}/{}'.format(30,favo_th))

    fig.savefig(
        os.path.join(folder, 'Mean_post_CR_cts_vs_run_for_ths_{}_{}'.format(30,favo_th)),
        format='png')

    ax3.plot(np.arange(runs)+1, favo_cr_pass,'o')
    ax3.set_xlabel('run number')
    ax3.set_ylabel('percentage CR passes')
    ax3.set_title('\n'+ 'CR passes vs. run for th:{}/{}'.format(30,favo_th))


    fig.savefig(
        os.path.join(folder, 'CR_passes_vs_run_for_ths_{}_{}'.format(30,favo_th)),
        format='png')





a.finish()
b.finish() 
        
               
