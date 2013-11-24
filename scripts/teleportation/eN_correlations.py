import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error, tomography

timestamp = None

states = ['$\Phi^+$', '$\Phi^-$', '$\Psi^+$', '$\Psi^-$']
correlation_names = ['0,0', '0,-1', '-1,0', '-1,-1']
fc = 'RoyalBlue'

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_correlations(name='adwindata')

corr = a.normalized_correlations
u_corr = a.u_normalized_correlations

for i in range(corr.shape[0]):
    
    fig = a.default_fig()
    ax = a.default_ax()

    pos = np.arange(4)
    ax.bar(pos, corr[i,:], 0.8, color=fc, yerr=u_corr[i,:], ecolor='k', 
        capsize=5)

    ax.set_xticks(pos + 0.4)
    ax.set_xlim(-0.1, 3.9)
    ax.set_xticklabels(correlation_names)
    ax.set_ylim(0,1)

    ax.set_ylabel('Probability')
    ax.set_xlabel('Correlation ($m_I$ (N), $m_s$ (e))')

    if len(states) == corr.shape[0]:
        ax.text(0.2, 0.8, states[i], fontsize=20)

a.finish()
