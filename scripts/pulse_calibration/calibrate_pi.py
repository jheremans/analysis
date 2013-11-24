import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from measurement.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error

import calibration_funcs as funcs

### parameters
timestamp = None
msmt_type = 'sequence'
guess_x0 = 0.08
guess_of = 1
guess_a = 0

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()

if msmt_type == 'sequence':
    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results('ssro')
    a.get_electron_ROC()
    ax = a.plot_result_vs_sweepparam(ret='ax')

    x = a.sweep_pts
    y = a.p0

elif msmt_type == 'mbi':
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    ax = a.plot_results_vs_sweepparam(name='adwindata', ret='ax')

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

else:
    raise Exception('Unknown msmt type')

res = funcs.calibrate_pulse_amplitude(x, y, ax, guess_x0, guess_of, guess_a)

