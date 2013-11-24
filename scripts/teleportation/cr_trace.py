import h5py
import numpy as np
from matplotlib import pyplot as plt


def get_pp_traces(fn, grp=None, prefix=None):
    f = h5py.File(fn, 'r')
    if grp == None:
        grp = f.keys()[0]

    if prefix == None:
        prefix = grp

    pump = f[grp+'/'+prefix+'_pump'].value
    probe = f[grp+'/'+prefix+'_probe'].value

    return pump, probe

def plot_pp_traces(pump, probe):
    fig, ax1 = plt.subplots(1,1, figsize=(8,4))
    x = np.arange(len(pump))

    ax1.plot(x, pump, '-b', drawstyle='steps', label='pump counts')
    ax1.set_xlabel('cycle')
    ax1.set_ylabel('pump counts', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    ax2 = ax1.twinx()
    ax2.plot(x, probe, '-r', drawstyle='steps', label='probe counts')
    ax2.set_ylabel('probe counts', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')


FN = r'D:\measuring\data\20130906\010233_CR_counter_trace_hans4_green_no_gate\010233_CR_counter_trace_hans4_green_no_gate.hdf5'
pump, probe = get_pp_traces(FN)

plot_pp_traces(pump, probe)


