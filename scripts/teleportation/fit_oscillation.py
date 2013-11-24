import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.fitting import fit
from analysis.lib.tools import plot
from analysis.lib.m2.ssro import mbi
reload(mbi)

from analysis.lib.tools import toolbox

timestamp = None#'20130710134201'

g_f = 1./180#2.19290
g_A = 0.5
g_o = 0.5
g_x0 = 0# 191

f = fit.Parameter(g_f, 'f')
A = fit.Parameter(g_A, 'A')
o = fit.Parameter(g_o, 'o')
x0 = fit.Parameter(g_x0, 'x0')

def fitfunc(x):
    return o() + A() * np.cos(2*pi*(f()*(x - x0())))

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax', )

x = a.sweep_pts
y = a.p0[:,0] #a.p0[:,1]   #.reshape(-1)

fit_result = fit.fit1d(x, y, None,
        fitfunc=fitfunc, p0=[f,x0,A,o], fixed=[], do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(a.sweep_pts[0],a.sweep_pts[-1],201), ax=ax,
        plot_data=False)

