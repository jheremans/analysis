import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi

timestamp = None

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('MBIElectronRabi')

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
#a.get_N_ROC(0.99, 0.02, 0.94, 0.01, 0.96, 0.01)
ax = a.plot_results_vs_sweepparam(ret='ax', )

x = a.sweep_pts.reshape(-1)[1:-1]
y = a.p0.reshape(-1)[1:-1]

x0 = fit.Parameter(170, 'x0') #Guess of where this is the right value
of = fit.Parameter(0.68136, 'of') #This is the value we need for 60deg, 300 deg and 420 deg CORPSE
a = fit.Parameter(0.3/60, 'a')
fitfunc_str = '-a (x-x0)'

def fitfunc(x):
    return of()-a()*(x-x0()) #Only plus or minus depends on which part of corpse

fit_result = fit.fit1d(x,y, None, p0=[x0,a], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False)


