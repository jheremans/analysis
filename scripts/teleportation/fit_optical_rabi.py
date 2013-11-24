import os, sys
import numpy as np
import h5py
import logging
from matplotlib import pyplot as plt

from analysis.lib import fitting
from analysis.lib.fitting import fit, rabi
from analysis.lib.tools import plot

from analysis.lib.tools import toolbox

timestamp = None#'20130827202900'

x_i = 0
x_f = 625-1
pts = 625

g_f = 1./10#8.593e-3#2.19290
g_A = 800
g_o = 800
g_x0 = 21
g_tau = 14.

f = fit.Parameter(g_f, 'f')
A = fit.Parameter(g_A, 'A')
o = fit.Parameter(g_o, 'o')
x0 = fit.Parameter(g_x0, 'x0')
tau = fit.Parameter(g_tau, 'c')


if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data()

empty,data = loadtxt(os.path.join(folder,'121500_th_30_3_CR_P_500_1.txt'), unpack = True)
#data = []
#for line in datafile:
#    li = line.strip()
#    if li.startswith(#)
#    data.append(line.strip().split())

#datafile.close


fig,ax = plt.subplots(1,1, figsize=(5,4))

x = np.linspace(x_i,x_f,pts)[76:200]*0.256
y = data[76:200]

plt.plot(x,y)

fit_result = fit.fit1d(x, y, rabi.fit_rabi_damped_exp_with_offset, 
		g_f, g_A, g_o, g_tau, g_x0, fixed=[], do_print=True, ret=True)
plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=True)

