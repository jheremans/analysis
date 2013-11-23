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
from analysis.lib.m2.ssro import sequence
reload(sequence)

folder = None
timestamp = None#'095551'#'180057'# None

if folder == None:
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('AdwinSSRO')

a = ssro.SSROAnalysis(folder)
b = sequence.SequenceAnalysis(folder)

runs = 1
debug = False
if debug:
    runs =2 

max_ts = [  ] 
taus = [  ] 
u_taus =[  ] 
means = [ ] 
u_means = []
percentage_passes = [ ] 

# analyze the data

for i,g in enumerate(a.g.items()):    
    gn = g[0]
    if 'instrument_settings' in gn:
        continue

    a.get_run(gn)
    max_t = int(string.split(gn, '_')[2])

    _t, _c = a.readout_relaxation(a.ro_time, a.ro_counts, a.reps, a.binsize, name=gn,
            plot=False, ret=True)

    idx0 = argmax(_c)
    idx1 = -1 
    t,c = _t[idx0:idx1]/1e3, _c[idx0:idx1]

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot(t,c, 'o')

    res = fit.fit1d(t,c, common.fit_exp_decay_with_offset, 0, c[0], 10, 
            ret=True, do_print=False, fixed=[])
    
    if res != False:
        plot.plot_fit1d(res, t, ax=ax, plot_data=False)
        max_ts.append(max_t)
        taus.append(res['params_dict']['tau'])
        u_taus.append(res['error_dict']['tau'])

    ax.set_xlabel('t ($\mu$s)')
    ax.set_ylabel('cts (Hz)')
    ax.set_title(a.default_plot_title + ', ' + gn)
    fig.savefig(os.path.join(folder, 'sp_relaxation_t_max={}'.format(max_t)))

    b.get_cr_results(gn, plot=False)

    hist, mean, u_mean = b.get_mean_cr_cts()
    means.append(mean[0])
    u_means.append(u_mean[0])
    
    stats = b.adwingrp(gn)['statistics'].value
    fail = stats[2]
    pts = np.linspace(5000,10000,21)[i-1]
    percentage_pass = (pts/(pts + fail)) * 100 #5000 is the number of succesfull measurements.
    percentage_passes.append(percentage_pass)

    plt.close('all')


# make plots

sortidxs = argsort(max_ts)
max_ts = np.array(max_ts)[sortidxs][:]
taus = np.array(taus)[sortidxs][:]
u_taus = np.array(u_taus)[sortidxs][:]
means= np.array(means)[sortidxs][:]
percentage_passes = np.array(percentage_passes)[sortidxs][:]

print percentage_passes
print max_ts

ks = 1./taus * 1e3
u_ks = 1./taus**2 * u_taus * 1e3

fig, ([ax1,ax2,ax3]) = plt.subplots(3,1, figsize = (8,10))
#ax1 = plt.subplot2grid((2,2),(0,0))
ax1.errorbar(max_ts, ks, yerr=u_ks, fmt='o', capsize=0, elinewidth=2)

ax1.set_xlabel('max t (us)')
ax1.set_ylabel('relaxation rate (kHz)')
ax1.set_title(a.timestamp + '\n relaxation rate')

fig.savefig(os.path.join(folder, 'Relaxation_rate'))


#fig = a.default_fig(figsize=(6,4))
#ax2 = plt.subplot2grid((2,2),(0,1))#a.default_ax(fig)
ax2.plot(max_ts, means,'o')
ax2.set_xlabel('max t (us)' )
ax2.set_ylabel('mean CR counts after sequence')
ax2.set_title('\n'+a.timestamp + ' mean CR cts')

fig.savefig(
    os.path.join(folder, 'post-CR_sum_vs_sweepparam.png'),
    format='png')

#fig = a.default_fig(figsize=(6,4))
#ax3 = plt.subplot2grid((2,2),(1,1),colspan =2)#a.default_ax(fig)
ax3.plot(max_ts, percentage_passes,'o')
ax3.set_xlabel('max t (us)')
ax3.set_ylabel('percentage CR passes')
ax3.set_title('\n'+a.timestamp + ' CR pass prob')

fig.savefig(
    os.path.join(folder, 'percentage_CR_pass_vs_sweepparam.png'),
    format='png')



a.finish()
b.finish() 
        
               
