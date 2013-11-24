import os, sys, time
import pickle
import pprint
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import ImageGrid

from analysis.lib.fitting import fit
from analysis import config

### configure stuff
config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'

rcParams['mathtext.default'] = 'regular'
rcParams['legend.numpoints'] = 1
rcParams['legend.frameon'] = False

fidspacing = 3 # if n, then every n-th datapoint from the fidelity data is used

### tools
def idx(arr, val):
    return argmin(abs(arr-val))

### get the data
fidsrcdir = os.path.join(config.outputdir, '20121102-ldefidelity')
savedir = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-ldefidelity')


fid = {}
f = np.load(os.path.join(fidsrcdir, 'fidelities.npz'))
for k in f.keys():
    fid[k] = f[k]
f.close()

corr = {}
f = np.load(os.path.join(fidsrcdir, 'correlations.npz'))
for k in f.keys():
    corr[k] = f[k]
f.close()

# extract the parity for some window slices, keep ch0, ch1 constant
# parity := P(00+11) - P(01+10)
parity = {
        'psi1' : {
            'ZZ' : [],
            'XX' : [],
            'XmX' : [],
            },
        'psi2' : {
            'ZZ' : [],
            'XX' : [],
            'XmX' : [],
            },
        }
u_parity = {
        'psi1' : {
            'ZZ' : [],
            'XX' : [],
            'XmX' : [],
            },
        'psi2' : {
            'ZZ' : [],
            'XX' : [],
            'XmX' : [],
            },
        }

ch0start = 639
ch1start = 668
windows = [60, 150]

ch0idx = idx(fid['ch0starts'], ch0start)
ch1idx = idx(fid['ch1starts'], ch1start)
winidxs = [ idx(fid['winvals'], w) for w in windows ]

ns = ['correctedpsi1correlations', 'correctedpsi2correlations']
uns = ['u_'+ns[0], 'u_'+ns[1]]
colors = ['k', 'r', 'b']

fig, axs = plt.subplots(2,2, figsize=(10,10))

for k,widx in enumerate(winidxs):
    for i,state in enumerate(['psi1', 'psi2']):
        for j,base in enumerate(['ZZ', 'XX', 'XmX']):
            parity[state][base].append( 2.*corr[ns[i]][:,widx,ch0idx,ch1idx,j,[0,3]].sum(-1)-1. )
            u_parity[state][base].append( np.sqrt( 4*corr[uns[i]][:,widx,ch0idx,ch1idx,j,0]**2 + \
                    4*corr[uns[i]][:,widx,ch0idx,ch1idx,j,3]**2) )

            axs[i,k].errorbar(fid['dtvals'], parity[state][base][-1], 
                    yerr=u_parity[state][base][-1], fmt='o', color=colors[j],
                    label=base)

            if i == 1:
                axs[i,k].set_xlabel('dt')
            if k == 0:
                axs[i,k].set_ylabel('P(even)-P(odd)')

            axs[i,k].set_xlim(0,150)
            axs[i,k].set_ylim(-1.1,1.1)
            axs[i,k].set_title(state+', window = %d' % windows[k])

axs[0,1].legend()




