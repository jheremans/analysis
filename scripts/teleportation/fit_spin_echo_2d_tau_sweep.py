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

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


timestamp = None
pts = 3

fidelities = []
tau_pi_pi = []

for i in np.arange(pts):

    name = 'twod_tau_sweep_{}_'.format(i)
    
    folder = toolbox.latest_data(name)

    a = sequence.SequenceAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='ssro')
    a.get_electron_ROC()
    ax = a.plot_result_vs_sweepparam(ret='ax', name='ssro')
    
    fid = a.p0.reshape(-1)[:]

    if i == 0:
        fidelities = fid
    else:
        fidelities = np.vstack((fidelities,fid))

    tau_pi_pi = tau_pi_pi.append(a.g.attrs['tau_pi_to_pi']/1e-6)

tau_pi2_pi = a.sweep_pts.reshape(-1)[:]
tau_pi2_pi, tau_pi_pi = np.meshgrid(tau_pi2_pi, tau_pi_pi)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.xlabel('pi/2 - pi tau (us)')
ax.ylabel('pi - pi tau (us)')

surface = ax.plot_surface(deltat, fet , fidelities, rstride=1, cstride=1, cmap=cm.coolwarm)

fig.colorbar(surface, aspect=5)


# testing

#deltat = np.linspace(-400,100,pts)
#fet = np.linspace(-1,1,pts)

#deltat,fet = np.meshgrid(fet, deltat)

#for i in np.arange(pts):
#    if i == 0:
#        fidelities = np.ones(pts)*i
#   else:
#        fidelities = np.vstack((fidelities,np.ones(pts)*i))

