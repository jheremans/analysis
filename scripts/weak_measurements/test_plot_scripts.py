import sys, os, time
import numpy as np
import matplotlib
#matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib import rcParams

import plots

LT1_clr = plots.colors['LT1']
LT2_clr = plots.colors['LT2']
common_clr = plots.colors['common']

rcParams['figure.subplot.hspace'] = 0.1
rcParams['figure.subplot.right'] = 0.95
rcParams['figure.subplot.top'] = 0.95
rcParams['figure.subplot.left'] = 0.15
rcParams['figure.subplot.bottom'] = 0.2

basepath = '../analysis_output/20121214_tpqi_manydts'
basepath = r'D:\machielblok\Desktop\PhD\QTlab\data\output'


fig,(ax1,ax2,ax3) = plt.subplots(1,3,sharey=True,figsize=(3,1.5)) 

ax1.plot([1,2,3],[1,2,3],color=common_clr['dark'])
ax2.plot([1,2,3],[1,4,9],color=common_clr['bright'])
ax3.plot([1,2,3],[1,2,4],color=common_clr['dark'])

ax1.set_ylabel(r'$g^(2)(\delta\tau)$')
ax2.set_xlabel(r'$\tau$ (ns)')
plt.subplots_adjust(wspace=0.15)

fig.savefig('tpqi.pdf', format='pdf')


