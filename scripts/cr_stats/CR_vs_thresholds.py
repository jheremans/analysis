import os, sys
import numpy as np
from matplotlib import pyplot as plt

import h5py
import logging

from analysis.lib.m2 import m2
from measurement.lib.tools import toolbox
from analysis.lib.m2.ssro import sequence

### params
timestamp = None#'095551'# '181958'#'20130908164621' # None


### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()


pts = 10
pres = np.linspace(3,30,pts)#np.ones(pts)*30#np.ones(pts)*30#np.linspace(6,46,pts)np.linspace(10,30,pts)
pros = pres#[1,2,3,4,5,6,8,10,15,20,25,30]# np.linspace(5,50,pts)#
analyze_probe = False
if analyze_probe:
    ths = pros
else:
    ths = pres
means = []
percentage_passes = []

for i,th in enumerate(ths):
    
    if analyze_probe:
        name = 'th_pres_{}_probe_{}'.format(pres[0],th)
    else:
        name = 'th_pres_{}_probe_{}'.format(th,th)
    
    a = sequence.SequenceAnalysis(folder)
    a.get_cr_results(name, plot=True)
    plt.close('all')

    mean = a.get_mean_cr_cts()
    means.append(mean)
    
    stats = a.adwingrp(name)['statistics'].value
    fail = stats[2]
    percentage_pass = 5000./(5000. + fail) * 100 #5000 is the number of succesfull measurements.
    percentage_passes.append(percentage_pass)

    a.finish()

fig = a.default_fig(figsize=(6,4))
ax = a.default_ax(fig)
ax.plot(ths, means,'o')
ax.set_xlabel('sweep threshold')
ax.set_ylabel('mean CR counts after sequence')
if save:
    fig.savefig(
        os.path.join(folder, 'post-CR_sum_vs_sweepparam.png'),
        format='png')

fig = a.default_fig(figsize=(6,4))
ax = a.default_ax(fig)
ax.plot(ths[:], percentage_passes[:],'o')
ax.set_xlabel('sweep threshold')
ax.set_ylabel('percentage CR passes')
if save:
    fig.savefig(
        os.path.join(folder, 'percentage_CR_pass_vs_sweepparam.png'),
        format='png')