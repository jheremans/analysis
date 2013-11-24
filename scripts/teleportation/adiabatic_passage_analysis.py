import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error

timestamp = None # '20130107231602'

o = fit.Parameter(0, 'o')
A = fit.Parameter(1, 'A')
F0 = fit.Parameter(1e6, 'F0')

def fitfunc(x) : 
    return o() + (A() - o()) * exp(-2*pi*(F0() ** 2) * x**-1)
fitfunc_str = 'o + (A-o) * exp(-2pi F0**2 / (df/dt)'
p0 = [o, F0]

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('AdiabaticPassage')

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax', )

x = a.sweep_pts.reshape(-1)[:][::-1]
y = a.p0.reshape(-1)[:][::-1]

f0 = a.g.attrs['passage_start_mod_frq']
f1 = a.g.attrs['passage_stop_mod_frq']
frate = (f1 - f0)/ (x*1e-6) # this is now in SI units -- Hz/s
frate2 = frate * 1e-12 # this is now MHz/us

fig = a.default_fig(figsize=(6,4))
ax = a.default_ax(fig)

ax.plot(frate, y, 'o')
ax.set_ylabel('P(ms=0)')
ax.set_xlabel('df/dt (MHz/us)')

fit_result = fit.fit1d(frate, y, None, p0=p0, fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(frate[0],frate[-1],201), ax=ax,
        plot_data=False)
xticks = ax.get_xticks()
ax.set_xticklabels(xticks*1e-12)