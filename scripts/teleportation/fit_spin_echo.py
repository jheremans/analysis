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

timestamp = '20130529132217'

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('CORPSEEcho')

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
#a.get_N_ROC(1.00, 0.02, 0.94, 0.01, 0.96, 0.01)
ax = a.plot_results_vs_sweepparam(ret='ax', name='adwindata')

x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

x0 = fit.Parameter(5000, 'x0')
a = fit.Parameter(1., 'a')
o = fit.Parameter(0.5, 'o')
c = fit.Parameter(1000, 'c')
f = fit.Parameter(0.0004, 'f')
phi = fit.Parameter(0., 'phi')
fitfunc_str = 'o + a * exp((x-x0)**2/c**2) * cos(2*pi*(f*x + phi))'

def fitfunc(x):
    return o() + a() * np.exp( -(x-x0())**2 / (c()**2)) * np.cos(2*pi*(f()*x + phi()))

fit_result = fit.fit1d(x,y, None, p0=[x0,a,c,f,phi], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False)


