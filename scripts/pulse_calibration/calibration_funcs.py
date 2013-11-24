import os, sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pprint

from analysis.lib.fitting import fit
from analysis.lib.tools import plot
from measurement.lib.tools import toolbox

def calibrate_pulse_amplitude(x, y, ax, *args):
    x0 = fit.Parameter(args[0], 'x0')
    of = fit.Parameter(args[1], 'of')
    a = fit.Parameter(args[2], 'a')
    fitfunc_str = '(1-of) + a (x-x0)**2'
    
    def fitfunc_parabolic(x):
        return (1.-of()) + a() * (x-x0())**2
    
    fit_result = fit.fit1d(x,y, None, p0=[of, a, x0], fitfunc=fitfunc_parabolic,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result