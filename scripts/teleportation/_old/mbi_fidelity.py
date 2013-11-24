import os, sys
import numpy as np
import h5py
import logging
import sympy

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi
from measurement.lib.tools import toolbox
from analysis.lib.tools import plot
from analysis.lib.math import error

### settings
timestamp = None # 

        
### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('MBI_fidelity')

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax',name='adwindata')

x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]
u_y = a.u_p0.reshape(-1)[:]

ax.set_ylim(0.4,1.05)

for (_x,_y,_u_y) in zip(x,y,u_y):
    ax.annotate("{:.3f} +/- {:.3f}".format(_y,_u_y), (_x,_y), 
        xytext=(0,-10), textcoords='offset points',
        ha='center', va='top')

pol = error.Formula()
a0, am1, ap1 = sympy.symbols('a0, am1, ap1')
pol.formula = am1 / (a0 + ap1 + am1)
pol.values[a0] = y[3]-y[1]
pol.values[am1] = y[3]-y[0]
pol.values[ap1] = y[3]-y[2]
pol.uncertainties[a0] = u_y[1]
pol.uncertainties[am1] = u_y[0]
pol.uncertainties[ap1] = u_y[2]

print 'Spin polarization = %.3f +/- %.3f' \
        % (float(pol.value()), float(pol.uncertainty()))
