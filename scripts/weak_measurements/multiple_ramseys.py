import os
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common,ramsey
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc
from matplotlib import rc

import plots

basepath=r'D:/machielblok/Desktop/PhD/QTlab/data/output'
name='ramsey'
fs=plots.fontsize

#datafolders=['1154','1252','1258','1303','1310','1347','1351','1316','1326', '1453']
#RO_time=[0,1,2,3,4,5,6,7,9, 11]
#
#
#datafolders=['164825','165126']
#date='20130111'

#datafolders=['111709','112023']
#date='20130115'

datafolders=['183118','183446']
date='20130116'
label=['mI=0','mI=+1']
meas_folder=r'D:\measuring\data'

datafolders=['190058','190909']
date='20130126'
label=['mI=-1','mI=0']
meas_folder=r'D:\measuring\data'

datafolders=['013420','014025']
date='20130228'
label=['mI=-1','mI=0']
meas_folder=r'D:\measuring\data'

datafolders=['120900','123020']
date='20130313'
label=['mI=-1','mI=0']
meas_folder=r'D:\machielblok\Desktop\PhD\QTlab\data'
def analysis():
    dp=os.path.join(meas_folder, date)
    amp=[]
    phase=[]
 
    j=0
    color=['Crimson','RoyalBlue']
    ramsey_data={}
    for i in datafolders:    
        tau_guess = 3000
        offset_guess=0.5   
        freq=1e-3
        amp=0.5
        result=sc.analyse_plot_results_vs_sweepparam(i,yname='P(mI=0)',Nuclcor=False,dataname='Spin_RO',d=date)
        resultcor=sc.analyse_plot_results_vs_sweepparam(i,yname='P(mI=0)',Nuclcor=True,dataname='Spin_RO',d=date)
	fit_result = fit.fit1d(result['x'], result['y'], ramsey.fit_ramsey_gaussian_decay, 
                    tau_guess, offset_guess, (freq,amp,90),
                    do_print = True , ret = True)
        fit_x=np.linspace(result['x'].min(),result['x'].max(),201)
	ramsey_data['fit_x']=fit_x
	if label[j]=='mI=-1':
	    ramsey_data['datamin1']=result['y']
	    ramsey_data['xmin1']=result['x']
	    ramsey_data['udatamin1']=result['uy']
	    ramsey_data['fitmin1']=fit_result['fitfunc'](fit_x)
	    ramsey_data['datamin1_cor']=resultcor['y']
        else:    
	    ramsey_data['datazero']=result['y']
	    ramsey_data['xzero']=result['x']
	    ramsey_data['udatazero']=result['uy']	
	    ramsey_data['fitzero']=fit_result['fitfunc'](fit_x)
	    ramsey_data['datazero1_cor']=resultcor['y']
	np.savez(os.path.join(basepath,name),**ramsey_data)
        j+=1 
	
def plot():	
    d=np.load(os.path.join(basepath,name)+'.npz')	
	#plot.plot_fit1d(fit_result,np.linspace(result['x'].min(),result['x'].max(),201))
        #measstrent=(90*result['x']/229 +0)/90
        #t=result['x']
        #result['x']=measstrent
        #result['x']=result['x']/191
        
    #plt.close('all')
    fig=plt.figure(figsize=[3,2])
    fig.clf()
    ax=fig.add_subplot(111)
    ax.errorbar(d['xzero'],d['datazero'],yerr=d['udatazero'],
		   mfc=plots.colors['N_zero'], mec=plots.colors['N_zero'],ecolor=plots.colors['N_zero'],
		   marker='o', linestyle='None',
		   label='  |0>')
    ax.plot(d['fit_x'],d['fitzero'],color=plots.colors['N_zero'],ls='-',linewidth=1)
   
    ax.errorbar(d['xmin1'],d['datamin1'],yerr=d['udatamin1'],
		   mfc=plots.colors['N_one'],mec=plots.colors['N_one'],ecolor=plots.colors['N_one'],
		   marker='o', linestyle='None',
		   label='  |1>')

    ax.plot(d['fit_x'],d['fitmin1'],color=plots.colors['N_one'],ls='-',linewidth=1)
    


    gcf().subplots_adjust(bottom=0.15,top=0.85,left=0.15,right=0.85)
    ax.set_xlabel (' free evolution time [us]', fontsize = fs)
    ax.set_ylabel ('P(ms=0)', fontsize = fs)   
    ax.set_ylim ([0, 1.1])
    ax.set_xlim ([0, 2000])
    yticks=[0,0.5,1]
    xticks=[0,500,1000,1500,2000]
    xticklabels=[0,0.5,1,1.5,2]
    ax3=ax.twiny()
    ax2=ax.twinx()
    ax.set_yticks(yticks)
    ax.set_xticks(xticks)
    ax.set_yticklabels(yticks,fontsize=fs)
    ax.set_xticklabels(xticklabels,fontsize=fs)
    
    
    min=d['datamin1'].min()
    max=d['datazero'].max()
    zero=0.5
    ax2.set_yticks([0,min,zero,max,1.1])
    ax2.set_yticklabels(['','-1','0','1',''],fontsize=fs)
    ax2.set_ylabel ('<Sz>N', fontsize = fs)

    
    theta=[90,180,270,360]
    t=[]
    t.append(0)
    for k in theta:
        t.append(int((k*229./90.)-12))
	print t
    t.append(2000)	
    ax3.set_xticks(t)
    ax3.set_xticklabels(['5','90','180','270','360',''],fontsize=fs)
    ax3.set_xlabel (r'$\theta$ (degrees)', fontsize = fs)
    
    ax.legend(loc=4,prop={'size':fs})
    fig.savefig(os.path.join(basepath,name+'.pdf'),format='pdf')
    
#analysis()
plot()
