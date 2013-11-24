import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi
from analysis.lib.fitting import fit
from measurement.lib.tools import toolbox
from analysis.lib.m2 import m2
from analysis.lib.m2.ssro import ssro
from analysis.lib.tools import plot

timestamp = None

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('165731')

a = sequence.SequenceAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='ssro')
a.get_electron_ROC()
ax = a.plot_result_vs_sweepparam(ret='ax', name='ssro')

x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

a = fit.Parameter(0, 'a')
o = fit.Parameter(0.5, 'o')
fitfunc_str = ''

def fitfunc_linear(x):
    return o() + a() * x
fit_result = fit.fit1d(x,y, None, p0=[o,a], fixed= [1], fitfunc=fitfunc_linear,
    fitfunc_str=fitfunc_str, do_print = True, ret=True)

plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
    plot_data=False)

