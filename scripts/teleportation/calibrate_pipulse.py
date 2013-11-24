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

timestamp = None # '20130605232959'

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
#a.get_N_ROC(1.00, 0.02, 0.94, 0.01, 0.96, 0.01)
ax = a.plot_results_vs_sweepparam(ret='ax', name='adwindata')

x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

x0 = fit.Parameter(0.02, 'x0')
of = fit.Parameter(1., 'of')
a = fit.Parameter(0., 'a')
fitfunc_str = '(1-of) + a (x-x0)**2'

def fitfunc(x):
    return (1.-of()) + a() * (x-x0())**2

fit_result = fit.fit1d(x,y, None, p0=[of, a, x0], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False)


