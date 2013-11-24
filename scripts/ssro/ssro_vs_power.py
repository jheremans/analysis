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
timestamp =  None

if folder == None:
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('AdwinSSRO')

a = ssro.SSROAnalysis(folder)

pwrs = []
taus = []
u_taus = []

for g in a.g.items():    
    gn = g[0]
    if 'instrument_settings' in gn:
        continue
    
    a.get_run(gn)
    #a.get_run('ms0')
   
    pwr = float(string.split(gn, '_')[-2])
    _t, _c = a.readout_relaxation(a.ro_time, a.ro_counts, a.reps, a.binsize, name=gn,
            plot=False, ret=True)

    idx0 = argmax(_c)+1
    idx1 = -1   
    t,c = _t[idx0:idx1]/1e3, _c[idx0:idx1]

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot(t,c, 'o')

    res = fit.fit1d(t,c, common.fit_exp_decay_with_offset, 0, c[0], 10, 
            ret=True, do_print=True, fixed=[])
    
    if res != False:
        plot.plot_fit1d(res, t, ax=ax, plot_data=False)
        pwrs.append(pwr)
        taus.append(res['params_dict']['tau'])
        u_taus.append(res['error_dict']['tau'])

    ax.set_xlabel('t ($\mu$s)')
    ax.set_ylabel('cts (Hz)')
    ax.set_title(a.default_plot_title + ', ' + gn)
    fig.savefig(os.path.join(folder, 'readout_relaxation_P=%dnW' % pwr))

plt.close('all')
a.finish()

sortidxs = argsort(pwrs)
pwrs = np.array(pwrs)[sortidxs][:]
taus = np.array(taus)[sortidxs][:]
u_taus = np.array(u_taus)[sortidxs][:]

ks = 1./taus * 1e3
u_ks = 1./taus**2 * u_taus * 1e3


res = fit.fit1d(pwrs, ks, common.fit_saturation, 120, 50,
        ret=True, do_print=True)

fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111)
ax.errorbar(pwrs, ks, yerr=u_ks, fmt='o', capsize=0, elinewidth=2)

if res != False:
    plot.plot_fit1d(res, np.linspace(0,pwrs[-1],101), ax=ax,
            plot_data=False)

ax.set_xlabel('RO power (nW)')
ax.set_ylabel('relaxation rate (kHz)')
ax.set_title(a.timestamp + '\nSaturation power')

fig.savefig(os.path.join(folder, 'saturation_power'))


