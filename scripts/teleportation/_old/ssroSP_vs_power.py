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
timestamp = None

if folder == None:
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data('AdwinSSRO_FTSP_vs_P')

a = ssro.SSROAnalysis(folder)
a.plot_format = 'png'

pwrs = []
fs = []
u_fs = []

tindex = 16

for g in a.g.items():
    gn = g[0]
    a.get_run(gn)
   
    pwr = int(string.split(gn, '_')[-1][:-2])
    fid_dat = a.fidelity(a.ro_counts, a.reps, a.binsize, 0, name=gn,
            plot=False, save=False, ret=True)
    
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    ax.errorbar(fid_dat[:,0], fid_dat[:,1], yerr=fid_dat[:,2],
            fmt='o', capsize=0, elinewidth=2)

    pwrs.append(pwr)
    fs.append(fid_dat[tindex,1])
    u_fs.append(fid_dat[tindex,2])

    ax.set_xlabel('t ($\mu$s)')
    ax.set_ylabel('F ($m_s$=0)')
    ax.set_title(a.default_plot_title + ', ' + gn)
    fig.savefig(os.path.join(folder, 'Fid_ms0_SPP=%dnW' % pwr))

    a.spinpumping(a.sp_time, a.sp_counts, a.reps, a.binsize, name=gn)

plt.close('all')
a.finish()

sortidxs = argsort(pwrs)
pwrs = np.array(pwrs)[sortidxs]

fig = plt.figure(figsize=(4,4))

ax = fig.add_subplot(111)
ax.errorbar(pwrs, fs, yerr=u_fs, fmt='o', capsize=0, elinewidth=2)

ax.set_xlabel('SP power (nW)')
ax.set_ylabel('(uncorrected) SP fidelity')
ax.set_title(a.timestamp + '\nSP fidelity')

fig.savefig(os.path.join(folder, 'SP_fidelity'))


