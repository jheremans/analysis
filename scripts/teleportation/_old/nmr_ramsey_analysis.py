import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, ramsey
from analysis.lib.tools import plot
from analysis.lib.math import error


timestamp = None # '20130107231602'
guess_f1 = 1./7000.
guess_A1 = 0.5
guess_phi1 = 0.

guess_tau = 30000
guess_a = 0.5


### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('NMRRamsey')

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results()
a.get_N_ROC(0.99, 0.02, 0.94, 0.005)
ax = a.plot_results_vs_sweepparam(ret='ax', )

fit_result = fit.fit1d(a.sweep_pts, a.p0.reshape(-1), ramsey.fit_ramsey_gaussian_decay,
        guess_tau, guess_a, (guess_f1, guess_A1, guess_phi1), fixed=[],
        do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(0,a.sweep_pts[-1],201), ax=ax,
        plot_data=False)

plt.savefig(os.path.join(folder, 'electronrabi_analysis.pdf'),
        format='pdf')

### FFT
# p0_fft = np.fft.fft(a.p0.reshape(-1), n=32)
# frq = np.fft.fftfreq(32, d=(a.sweep_pts[1]-a.sweep_pts[0])/1e3)
#  
# fig = a.default_fig()
# ax = a.default_ax(fig)
# ax.plot(frq, p0_fft, 'o-')
