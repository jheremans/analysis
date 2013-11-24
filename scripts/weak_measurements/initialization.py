import os
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
from analysis.lib.tools import data_handling as dh

import plots

basepath=r'D:/machielblok/Desktop/PhD/QTlab/data/output'
name='desr'
fs=plots.fontsize


def analysis():
	fres_desr,data_desr,udata_desr = sc.plot_dark_esr(sc.get_latest_data('dark',date='20130311'),d='20130311')
	fres_zero,data_zero,udata_zero = sc.plot_dark_esr(sc.get_latest_data('init_',date='20130311'),d='20130311')
	fres_min1,data_min1,udata_min1 = sc.plot_dark_esr(sc.get_latest_data('init_',date='20130305'),d='20130305')
	esr_data={}
	xfit=np.linspace(fres_desr['x'].min(),fres_desr['x'].max(),501)
	esr_data['xfit']=xfit
	esr_data['desrfit'] = fres_desr['fitfunc'](xfit)
	esr_data['zerofit'] = fres_zero['fitfunc'](xfit)
	esr_data['min1fit'] = fres_min1['fitfunc'](xfit)
	esr_data['xdesr'] = fres_desr['x']
	esr_data['xzero'] = fres_zero['x']
	esr_data['xmin1'] = fres_min1['x']

	esr_data['datadesr'] = data_desr
	esr_data['datazero'] = data_zero
	esr_data['datamin1'] =data_min1

	esr_data['udatadesr'] = udata_desr
	esr_data['udatazero'] = udata_zero
	esr_data['udatamin1'] = udata_min1

	print os.path.join(basepath,name)
	print type(esr_data)
	np.savez(os.path.join(basepath,name), **esr_data)

def plot():
   esr_data=np.load(os.path.join(basepath,name)+'.npz')	
   fig=plt.figure(figsize=[2,2])
   fig.clf()
   ax=fig.add_subplot(111)
   ax.errorbar(esr_data['xdesr'],esr_data['datadesr'],yerr=esr_data['udatadesr'],
		   mfc=plots.colors['common'], mec=plots.colors['common'],ecolor=plots.colors['common'],
		   marker='o', linestyle='None',
		   label='  No initialization')
   ax.plot(esr_data['xfit'],esr_data['desrfit'],color=plots.colors['common'],ls='-',linewidth=1)
   
   ax.errorbar(esr_data['xzero'],esr_data['datazero'],yerr=esr_data['udatazero'],
		   mfc=plots.colors['N_zero'],mec=plots.colors['N_zero'],ecolor=plots.colors['N_zero'],
		   marker='o', linestyle='None',
		   label='  mI=0')

   ax.plot(esr_data['xfit'],esr_data['zerofit'],color=plots.colors['N_zero'],ls='-',linewidth=1)

   ax.errorbar(esr_data['xmin1'],esr_data['datamin1'],yerr=esr_data['udatamin1'],
		   mfc=plots.colors['N_one'],mec=plots.colors['N_one'],color=plots.colors['N_one'],
		   marker='o',linestyle='None',
		   label='  mI=-1')
   ax.plot(esr_data['xfit'],esr_data['min1fit'],color=plots.colors['N_one'],ls='-',linewidth=1)

   ax.set_xlabel('MW freq (GHz)',fontsize=fs)
   ax.set_ylabel('P(ms=0)',fontsize=fs)
   ax.set_ylim([0.4,1.15])
   ax.set_xlim([2.826,2.832])
   ax.set_yticks([0.5,0.75,1])
   ax.set_xticks([2.826,2.829,2.832])
   ax.set_yticklabels([0.5,0.75,1],fontsize=fs)
   ax.set_xticklabels([2.826,2.829,2.832],fontsize=fs)
   ax.xaxis.set_major_formatter(FormatStrFormatter('%0.3f'))
   #ax.set_title('Dark ESR')
   ax.legend(loc=4,prop={'size':fs})
   fig.savefig(os.path.join(basepath,name+'.pdf'),format='pdf')
#analysis()   
plot()
