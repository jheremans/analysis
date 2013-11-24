import os, sys
import string
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib.fitting import fit, common
from measurement.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.m2.ssro import ssro


folder = None
timestamp =None#'095551'#'180057'# None

if folder == None:
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('AdwinSSRO')

a = ssro.SSROAnalysis(folder)
b = sequence.SequenceAnalysis(folder)

runs = 8
sweep_probes = 2 # 2 if also analyze_probe, 1 if only preselect.

ths = [ [ [ ] for i in range (runs)] for s in range (sweep_probes) ] 
taus = [ [ [ ] for i in range (runs)] for s in range (sweep_probes) ] 
u_taus =[ [ [ ] for i in range (runs)] for s in range (sweep_probes) ] 
means = [ [ [ ] for i in range (runs)] for s in range (sweep_probes) ] 
precentage_passes = [ [ [ ] for i in range (runs)] for s in range (sweep_probes) ] 

for g in a.g.items():    
    gn = g[0]
    if 'instrument_settings' in gn:
        continue
    
    a.get_run(gn)
    pre = float(string.split(gn, '_')[2])
    pro = float(string.split(gn, '_')[4])
    run = float(string.split(gn, '_')[6])
    analyze_probe = float(string.split(gn, '_')[8])
    
    if analyze_probe:
        th = pro
        probe = 1
    else:
        th = pre
        probe = 0

    _t, _c = a.readout_relaxation(a.ro_time, a.ro_counts, a.reps, a.binsize, name=gn,
            plot=False, ret=True)

    idx0 = argmax(_c)
    idx1 = -1   qq
    t,c = _t[idx0:idx1]/1e3, _c[idx0:idx1]

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot(t,c, 'o')

    res = fit.fit1d(t,c, common.fit_exp_decay_with_offset, 0, c[0], 10, 
            ret=True, do_print=False, fixed=[])
    
    if res != False:
        plot.plot_fit1d(res, t, ax=ax, plot_data=False)
        ths[probe][run].append(th)
        taus[probe][run].append(res['params_dict']['tau'])
        u_taus[probe][run].append(res['error_dict']['tau'])

    ax.set_xlabel('t ($\mu$s)')
    ax.set_ylabel('cts (Hz)')
    ax.set_title(a.default_plot_title + ', ' + gn)
    fig.savefig(os.path.join(folder, 'sp_relaxation_Th_pre={}_Th_pro={}_run_{}_sweep_{}'.format(int(pre),int(pro), run, th)))
    plt.close('all')


    b.get_cr_results(gn, plot=True)
    plt.close('all')

    mean = a.get_mean_cr_cts()
    means[probe][run].append(mean)
    
    stats = a.adwingrp(name)['statistics'].value
    fail = stats[2]
    percentage_pass = 5000./(5000. + fail) * 100 #5000 is the number of succesfull measurements.
    percentage_passes[probe][run].append(percentage_pass)


plt.close('all')
a.finish()
b.finish()

print ths
print taus
print means
print percentage_passes


for r in np.arange(runs):
    for a in np.arange(sweep_probes):
        if a == 1:
            p = probe
        else:
            p = preselect

        sortidxs = argsort(ths[a][r])
        ths[a][r] = np.array(ths[a][r])[sortidxs][:]
        taus[a][r] = np.array(taus[a][r])[sortidxs][:]
        u_taus[a][r] = np.array(u_taus[a][r])[sortidxs][:]

        ks[a][r] = 1./taus[a][r] * 1e3
        u_ks[a][r] = 1./taus[a][r]**2 * u_taus[a][r] * 1e3

        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
        ax.errorbar(ths[a][r], ks[a][r], yerr=u_ks[a][r], fmt='o', capsize=0, elinewidth=2)

        ax.set_xlabel('{} threshold'.format(p))
        ax.set_ylabel('relaxation rate (kHz)')
        ax.set_title(a.timestamp + '\n{} threshold for run {}'.format(p,r))

        fig.savefig(os.path.join(folder, '{}_threshold_calib_run_{}'.format(p,r)))


        fig = a.default_fig(figsize=(6,4))
        ax = a.default_ax(fig)
        ax.plot(ths, means[a][r],'o')
        ax.set_xlabel('sweep threshold')
        ax.set_ylabel('mean CR counts after sequence')

        fig.savefig(
            os.path.join(folder, 'post-CR_sum_vs_sweepparam.png'),
            format='png')

        fig = a.default_fig(figsize=(6,4))
        ax = a.default_ax(fig)
        ax.plot(ths[a][r], percentage_passes[a][r],'o')
        ax.set_xlabel('sweep threshold')
        ax.set_ylabel('percentage CR passes')

        fig.savefig(
            os.path.join(folder, 'percentage_CR_pass_vs_sweepparam.png'),
            format='png')