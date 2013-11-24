import numpy as np
import h5py

import T2_filter

### Analysis settings
FN = r'D:\measuring\data\20130823\002947_Teleportation_testing_lt2-9_optical_rabi\002947_Teleportation_testing_lt2-9_optical_rabi.hdf5'
sweep_pts = linspace(0.5,1.5,11)
sweep_name = 'EOM pulse amplitude (V)'
t0 = 10.145e6
t1 = 10.19e6
bins = 200

def get_hist(fn=FN):
    f = h5py.File(FN, 'r')
    channel = f['/HH_channel-1'].value
    special = f['/HH_special-1'].value
    sync_time = f['/HH_sync_time-1'].value
    f.close()
    h2d, binedges = T2_filter.hist2d(channel, special, sync_time,
        100, len(sweep_pts), t0, t1, bins, 1)
    
    return h2d, binedges


### script
h2d, binedges = get_hist()

fig, ax = subplots(1,1, figsize=((4,4)))
ax.imshow(h2d, origin='lower', interpolation='nearest', cmap='PiYG',
    extent=(t0/1e3, t1/1e3, 0, len(sweep_pts)), aspect='auto')

ax.set_yticks(np.arange(len(sweep_pts))+0.5)
ax.set_yticklabels(sweep_pts)

ax.set_xlabel('time (ns)')
ax.set_ylabel(sweep_name)
