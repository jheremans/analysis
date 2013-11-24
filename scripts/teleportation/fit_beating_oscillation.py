import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.fitting import fit
from analysis.lib.tools import plot
from analysis.lib.m2.ssro import sequence
reload(mbi)

from analysis.lib.tools import toolbox

timestamp = None#'20130729125612'

g_f = 1/110#2.19290
g_f2 = 1/90
g_f3 = 1/100
g_A = 0.3
g_A2=0.3
g_A3 = 0.3
g_o = 0.5
g_x0 = 0
g_x02 = 0
g_x03 = 0

f = fit.Parameter(g_f, 'f')
A = fit.Parameter(g_A, 'A')
o = fit.Parameter(g_o, 'o')
x0 = fit.Parameter(g_x0, 'x0')

f2 = fit.Parameter(g_f2,'f2')
A2 = fit.Parameter(g_A2, 'A2')
x02 = fit.Parameter(g_x02, 'x02')

f3 = fit.Parameter(g_f3,'f3')
A3 = fit.Parameter(g_A3, 'A3')
x03 = fit.Parameter(g_x03, 'x03')




def fitfunc(x):
    return o() + A() * np.cos(2*pi*(f()*(x - x0()))) + A2() * np.cos(2*pi*(f2()*(x-x02())))+ A3() * np.cos(2*pi*(f3()*(x-x03())))



if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('114816')

a = sequence.SequenceAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='ssro')
a.get_electron_ROC()
ax = a.plot_result_vs_sweepparam(ret='ax',name='ssro')

x = a.sweep_pts
y = a.p0 #a.p0[:,1]   #.reshape(-1)

fit_result = fit.fit1d(x, y, None,
        fitfunc=fitfunc, p0=[f,x0,A,o,f2,A2,x02,f3,A3,x03], fixed=[1,6], do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(a.sweep_pts[0],a.sweep_pts[-1],201), ax=ax,
        plot_data=False)

