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

#if timestamp != None:
#   folder = toolbox.data_from_time(timestamp)
#else:
#    folder = toolbox.latest_data()
    
#a = mbi.MBIAnalysis(folder)
#a.get_sweep_pts()
#a.get_readout_results()
#a.get_electron_ROC()
#a.get_N_ROC(0.99, 0.02, 0.94, 0.01, 0.96, 0.01)
#ax = a.plot_results_vs_sweepparam(ret='ax', )

x,y = loadtxt(r'D:\measuring\data\20130408\184036_Saturation_curve_Sil2-PSB_SM\184036_Saturation_curve_Sil2-PSB_SM.dat',unpack=True)

A = fit.Parameter(900000, 'x0')
xsat = fit.Parameter(100, 'of')
fitfunc_str = 'A * x / (x + xsat)'

def fitfunc(x):
    return A() * x / (x + xsat())
    
fig = plt.figure()
ax = fig.add_subplot(111)
    
fit_result = fit.fit1d(x,y, None, p0=[A,xsat], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=True)

ax.set_xlabel('Power ($\mu$W)')

