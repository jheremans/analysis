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

name='140933_LaserFrequencyScan_red_scan_sil2_power_1nW'

V,x_dat,y_dat = loadtxt(r'D:/measuring/data/20130416/' + name + '/' + name + '.dat',unpack=True)

x=x_dat[240:360]
y=y_dat[240:360]

a = fit.Parameter(40000, 'a')
of = fit.Parameter(800,'of')
f = fit.Parameter(57.2, 'f')
s = fit.Parameter(0.2, 's')
a2 = fit.Parameter(-20000, 'a2')
f2 = fit.Parameter(57.1, 'f2')
s2 = fit.Parameter(0.2, 's2')
fitfunc_str = 'of + a*exp(-(x-f)^2/2s^2)+  a*exp(-(x-f)^2/2s^2)'

def fitfunc(x):
    return of() + a() * np.exp( -(x - f())**2./(2.*s()**2.)) \
            + a2() * np.exp( -(x - f2())**2./(2.*s2()**2.))

#def fitfunc(x):
#    return of() + 1./(pi*s()*(1+((x-f())/s())**2))


fig = plt.figure()
ax = fig.add_subplot(111)
    
fit_result = fit.fit1d(x,y, None, p0=[a,of,f,s,a2,f2,s2], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],401), ax=ax,
        plot_data=True)

ax.set_xlabel('frequency (GHz)')

