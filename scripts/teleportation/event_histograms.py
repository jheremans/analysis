import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from measurement.lib.tools import toolbox

FOLDER = toolbox.latest_data('out-003')
FN = toolbox.measurement_filename(FOLDER)
BINEDGES = np.arange(3001)-0.5

def get_photon_times_ns(FN):
    f = h5py.File(FN, 'r')
    channel = f['/HH_channel-1'].value
    special = f['/HH_special-1'].value
    sync_time = f['/HH_sync_time-1'].value
    f.close()

    # TODO this could maybe be done more efficient with cython
    is_not_special = special==0
    is_channel_0 = channel==0
    is_channel_1 = channel==1

    is_photon_0 = np.logical_and(is_not_special, is_channel_0)
    is_photon_1 = np.logical_and(is_not_special, is_channel_1)

    return sync_time[is_photon_0]*1e-3, sync_time[is_photon_1]*1e-3

def make_hist(photons):
    h, b = np.histogram(photons, bins=BINEDGES)
    return h, b

photons_0, photons_1 = get_photon_times_ns(FN)
h0, b0 = make_hist(photons_0)

fig0, ax0 = plt.subplots(1,1)
ax0.semilogy(b0[:-1], h0)



