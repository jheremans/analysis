import os, sys
import numpy as np
import h5py
import logging

from analysis.lib.m2 import m2
from measurement.lib.tools import toolbox
from analysis.lib.m2.ssro import sequence

### params
timestamp = None # '20130908164621' # None


### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()


a = sequence.SequenceAnalysis(folder)
a.get_sweep_pts()
# a.get_readout_results(name)
a.get_cr_results('ssro')
a.plot_cr_vs_sweep(ionization_crit=3)
# a.get_electron_ROC()
# a.plot_result_vs_sweepparam()
a.finish()
