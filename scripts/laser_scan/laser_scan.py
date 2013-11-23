import os
import numpy as np
from matplotlib import pyplot as plt
from measurement.lib.tools import toolbox
from analysis.lib.nv import nvlevels

# settings
timestamp = None

Ex = 43.8
Ey = 53.17

lines = nvlevels.get_ES_ExEy(Ex,Ey)

### script
if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('LaserFrequencyScan')

# get data
fn = os.path.join(folder, os.path.split(folder)[1]+'.dat')
d = np.loadtxt(fn)

_folder, _file = os.path.split(fn)
_folder, _ = os.path.split(_folder)
_folder, _date = os.path.split(_folder)
title = _date+'/'+_file[:-4]

frq = d[:,1]
cts = d[:,2]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(frq, cts, 'k-', lw=2)

ax.set_xlabel('frequency (GHz) - 470.4 THz')
ax.set_ylabel('counts (Hz)')
ax.set_title(title)

y0,y1 = ax.get_ylim()
ax.vlines(lines, 0, y1, colors='r', lw=1)

ax.set_ylim(0,y1)
ax.text(lines[2], y1*0.95, '$E_y$', ha='center', va='top',
        color='r', backgroundcolor='w')
ax.text(lines[3], y1*0.95, '$E_x$', ha='center', va='top',
        color='r', backgroundcolor='w')

ax.set_xlim(min(frq), max(frq))


