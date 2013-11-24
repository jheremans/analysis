import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro, mbi
from analysis.lib.tools import toolbox
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot
from analysis.lib.math import error, tomography

# idealstate = np.array([1,0])#1/np.sqrt(2)*np.array([1,1])#

timestamp = None # '20130616155042'


### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
# a.get_N_ROC(0.95,0.03,0.04,0.03,0.98,0.01,0.99,0.01)#(0.99, 0.02, 0.94, 0.01, 0.96, 0.01)
a.get_electron_ROC()

a.save()    
# a.plot_results_vs_sweepparam()

mz, u_mz = a.p0.reshape(-1)[0], a.u_p0.reshape(-1)[0]
mx, u_mx = a.p0.reshape(-1)[1], a.u_p0.reshape(-1)[1]
my, u_my = a.p0.reshape(-1)[2], a.u_p0.reshape(-1)[2]

# dm, u_dm = tomography.measured_single_qubit_dm(mx,my,mz,u_mx,u_my,u_mz)
# fid, u_fid = tomography.fidelity(idealstate, dm, u_dm)

# plotting
# plot.dm(dm)
fig = a.default_fig()
ax = a.default_ax()

plot.msmnt_outcomes(np.array([mx,my,mz]), np.array([u_mx, u_my, u_mz]),
        xticks=['X', 'Y', 'Z'], ax=ax)

a.finish()
