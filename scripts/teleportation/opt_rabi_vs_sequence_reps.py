import numpy as np
import h5py
from matplotlib import plot as plt



### Analysis settings
FN = r'D:\measuring\data\20130830\002947_Teleportation_testing_lt1-4_optical_rabi\002947_Teleportation_testing_lt2-9_optical_rabi.hdf5'

r = 1 # after which rep to look



f = h5py.File(FN, 'r')
#channel = f['/HH_channel-1'].value
sync_time = f['/HH_sync_time-1'].value
sync_nr = f['/HH_sync_nr-1'].value
f.close()




hist = np.zeros(max(sync_time))
for i in sync_time[where(sync_nr%1000==r)]:
    hist[i] += 1

fig, ax = subplots(1,1, figsize=((4,4)))

plt.plot(hist)



ax.set_xlabel('sync time')
ax.set_ylabel('counts')
