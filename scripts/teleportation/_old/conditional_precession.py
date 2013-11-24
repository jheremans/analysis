import numpy as np
from matplotlib import pyplot as plt
from analysis.lib.m2.ssro import mbi
from measurement.lib.tools import toolbox

# settings
msmnt_folders = [
        '20130207175652',
        '20130207180726',
        ]

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(111)

for f in msmnt_folders:
    a = mbi.ConditionalPrecessionAnalysis(toolbox.data_from_time(f))
    a.get_sweep_pts()
    a.get_readout_results()

    of = min(a.normalized_ssro[:,0])
    ax.errorbar(a.sweep_pts, a.normalized_ssro[:,0], fmt='o-',
            yerr=a.u_normalized_ssro[:,0], label=a.name + ' (%s)' % f)

ax.legend()
ax.set_ylim(0.5,0.8)
ax.set_xlim(0,500)
ax.set_title('conditional precession')

ax.set_xlabel('delay ($\mu$s)')
ax.set_ylabel(r'readout results $|0\rangle$')
