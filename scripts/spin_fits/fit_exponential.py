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
    folder = toolbox.latest_data()

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax')

x = a.sweep_pts.reshape(-1)[:]
y = a.normalized_ssro.reshape(-1)[:]

a = fit.Parameter(1, 'a')
o = fit.Parameter(0., 'o')
c = fit.Parameter(1000, 'c')
fitfunc_str = 'o - a exp(-x/c)'

def fitfunc(x):
    return o() + a() * np.exp( -(x) / (c()))

fit_result = fit.fit1d(x,y, None, p0=[o,c,a], fixed=[], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False)


