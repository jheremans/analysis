import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.fitting import fit
from analysis.lib.tools import plot, toolbox
from analysis.lib.m2.ssro import mbi
reload(mbi)

timestamp = None#'20130710135342'

g_f = 1./360
g_A = 0.5
g_o = 0.5
g_x0= 40/180


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
a.get_correlations(name = 'adwindata')
ax = a.plot_results_vs_sweepparam(ret='ax', mode = 'correlations')

x = a.sweep_pts

y = a.normalized_correlations[:,2]

fit_result= fit.fit1d(x, y, None,
        fitfunc=fitfunc, p0=[x0,A,o,f], fixed=[], do_print=True, ret=True)




plot.plot_fit1d(fit_result, np.linspace(a.sweep_pts[0],a.sweep_pts[-1],201), ax=ax,
        plot_data=False)

