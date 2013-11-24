import os, sys
import numpy as np
import h5py
import logging
import sympy

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, esr
from analysis.lib.tools import plot
from analysis.lib.math import error

### settings
timestamp = None # 

guess_offset = 1.0
guess_A_min1 = 0.5
guess_A_plus1 = 0.
guess_A_0 = 0.
guess_x0 = 2827.093
guess_sigma = 0.435
guess_Nsplit = 2.187

splitting = 2.187



### fitfunction
A_min1 = fit.Parameter(guess_A_min1, 'A_min1')
A_plus1 = fit.Parameter(guess_A_plus1, 'A_plus1')
A_0 = fit.Parameter(guess_A_0, 'A_0')
o = fit.Parameter(guess_offset, 'o')
x0 = fit.Parameter(guess_x0, 'x0')
sigma = fit.Parameter(guess_sigma, 'sigma')
Nsplit = fit.Parameter(guess_Nsplit, 'Nsplit')

def fitfunc(x):
    # return o() - A_min1()*np.exp(-((x-(x0()-splitting-Nsplit()))/sigma())**2) \
    #         - A_min1()*np.exp(-((x-(x0()+splitting-Nsplit()))/sigma())**2) \
    #         - A_plus1()*np.exp(-((x-(x0()-splitting+Nsplit()))/sigma())**2) \
    #         - A_plus1()*np.exp(-((x-(x0()+splitting+Nsplit()))/sigma())**2) \
    #         - A_0()*np.exp(-((x-(x0()+Nsplit()))/sigma())**2) \
    #         - A_0()*np.exp(-((x-(x0()-Nsplit()))/sigma())**2) 
    return o() - A_min1()*np.exp(-((x-(x0()-Nsplit()))/sigma())**2) \
            - A_plus1()*np.exp(-((x-(x0()+Nsplit()))/sigma())**2) \
            - A_0()*np.exp(-((x-x0())/sigma())**2) \
          
### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('PostInitDarkESR')

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
ax = a.plot_results_vs_sweepparam(ret='ax',name='adwindata')

x = a.sweep_pts
y = a.p0.reshape(-1)

# try fitting
fit_result = fit.fit1d(x, y, None, p0 = [A_min1, A_plus1, A_0, sigma, o],
        fitfunc = fitfunc, do_print=True, ret=True, fixed=[])
plot.plot_fit1d(fit_result, x, plot_data=False, ax=ax)

ax.set_ylim(0.4,1.05)

plt.savefig(os.path.join(folder, 'mbi_darkesr_analysis.pdf'),
        format='pdf')

pol = error.Formula()
a0, am1, ap1 = sympy.symbols('a0, am1, ap1')
pol.formula = am1 / (a0 + ap1 + am1)
pol.values[a0] = A_0()
pol.values[am1] = A_min1()
pol.values[ap1] = A_plus1()
pol.uncertainties[a0] = fit_result['error_dict']['A_0']
pol.uncertainties[am1] = fit_result['error_dict']['A_min1']
pol.uncertainties[ap1] = fit_result['error_dict']['A_plus1']

print 'Spin polarization = %.3f +/- %.3f' \
        % (float(pol.value()), float(pol.uncertainty()))
