import os
import numpy as np
import h5py
from matplotlib import pyplot as plt
from measurement.lib.tools import toolbox

# settings
timestamp = None


### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('scan2d')

# get data
fn = os.path.join(folder, os.path.split(folder)[1]+'.hdf5')

d = h5py.File(fn, 'r')
x = d['x'].value
y = d['y'].value
z = d['countrate'].value
d.close()

fig = plt.figure()
ax = fig.add_subplot(111)

ax.imshow(z, origin='lower', interpolation='None')
ax.set_yticklabels([])
ax.set_xticklabels([])

xrng = abs(x[-1]-x[0])
yrng = abs(y[-1]-y[0])
xbar = 1.
xlen = xbar/xrng*len(x)
xstart = 0.1/xrng*len(x)
ystart = 0.1/yrng*len(x)

ax.hlines([ystart], xstart, xstart+xlen, colors='w', linestyles='solid', lw=3)
ax.text((xstart+xstart+xlen)/2., ystart*1.05, '%.1f $\mu$m' % xbar, color='w',
        ha='center', va='bottom')

fig.savefig(os.path.splitext(fn)[0]+'.png')

