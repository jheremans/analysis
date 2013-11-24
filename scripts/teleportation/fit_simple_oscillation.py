import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.fitting import fit
from analysis.lib.tools import plot
from analysis.lib.m2.ssro import sequence
reload(sequence)

from analysis.lib.tools import toolbox

timestamp = '171832'#None #'20130801132824'

g_f = 1./360.#8.593e-3#2.19290
g_A = 0.5
g_o = 0.5
g_x0 = 0


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

a = sequence.SequenceAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='ssro')
a.get_electron_ROC()
ax = a.plot_result_vs_sweepparam(ret='ax',name='ssro')

x = a.sweep_pts
y = a.p0 #a.p0[:,1]   #.reshape(-1)

fit_result = fit.fit1d(x, y, None,
        fitfunc=fitfunc, p0=[f,x0,A,o], fixed=[1], do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(a.sweep_pts[0],a.sweep_pts[-1],201), ax=ax,
        plot_data=False)

