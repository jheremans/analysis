
import numpy as np
import pylab as plt
from analysis.lib.spin import spin_control as sc
import os
from mpl_toolkits.mplot3d import Axes3D
from analysis.lib.tools import weaktools as tls
from analysis.lib.tools import plot as plot
from analysis.lib.math import tomography as tom

# 2 msmnts uncollapse very short
d='20130404'
foldernames=['192943','201204','212222','221836','233234','005033','015405','030325','035504']
zfoldernames=['183738','174132','180335','182237']
zd='20130403'
th=0
dir='down'

# 2 msmnts steer very short
d='20130404'
foldernames=['005342','012242','015106','022324','025507','033115','041030','044611','052016','055424']
zfoldernames=['154355']
zd='20130403'
th=''
dir='up'
'''
# 1 msmnt
d='20130405'
foldernames=['141755','144145','150814','153230','155539','161924','164248','171312']
zd=d
zfoldernames=['130443','173349']
th=''
'''
tau=50.
utau=1
for n in foldernames:
    datapath=sc.get_latest_data(n,date=d)
    
    files = os.listdir(datapath)
    f='Spin_RO'
    for k in files:
        if f in k:
            spin_ro_file = k

    data = np.load(datapath+'\\'+spin_ro_file)
    if (n==foldernames[0]):
        phase=np.zeros(len(data['FS']))
        rep=np.zeros(len(data['FS']))
    reps_array=data['SN']+data['FF']+data['FS']
    reps=reps_array[0]
    data_norm={}
    data_norm['sweep_par']=data['sweep_par']
    data_norm['sweep_par_name']=data['sweep_par_name']
    data_norm['SN']=data['SN']/(reps+0.)
    data_norm['FS']=data['FS']/(reps-data['SN']+0.)
    data_norm['FF']=data['FF']/(reps-data['SN']+0.)
    data_norm['uSN']=(data_norm['SN']*(1-data_norm['SN'])/reps)**.5
    data_norm['uFS']=(data_norm['FS']*(1-data_norm['FS'])/(reps-data['SN']+0.))**.5
    data_norm['uFF']=(data_norm['FF']*(1-data_norm['FF'])/(reps-data['SN']+0.))**.5

    data_norm['FinalRO_SN']=data['FinalRO_SN']/(data['SN']+0.)
    data_norm['FinalRO_FS']=data['FinalRO_FS']/(data['FS']+0.)
    data_norm['FinalRO_FF']=data['FinalRO_FF']/(data['FF']+0.)
    data_norm['FinalRO_Succes']=(data['FinalRO_SN']+data['FinalRO_FS'])/(data['SN']+data['FS']+0.)
    data_norm['FinalRO_All']=(data['FinalRO_SN']+data['FinalRO_FS']+data['FinalRO_FF'])/(reps+0.)

    print data['FS']
    for i in np.arange(len(data['SN'])):
        phase[i]+=data['FinalRO_FS'][i]
        rep[i]+=data['FS'][i]
    print rep
    phasenorm=phase/(rep+0.)
    x=data['sweep_par']
    data.close()
for n in zfoldernames:
    datapath=sc.get_latest_data(n,date=zd)
    
    files = os.listdir(datapath)
    f='Spin_RO'
    for k in files:
        if f in k:
            spin_ro_file = k

    data = np.load(datapath+'\\'+spin_ro_file)
    if (n==zfoldernames[0]):
        z=0
        zrep=0
    reps_array=data['SN']+data['FF']+data['FS']
    reps=reps_array[0]
    data_norm={}
    data_norm['sweep_par']=data['sweep_par']
    data_norm['sweep_par_name']=data['sweep_par_name']
    data_norm['SN']=data['SN']/(reps+0.)
    data_norm['FS']=data['FS']/(reps-data['SN']+0.)
    data_norm['FF']=data['FF']/(reps-data['SN']+0.)
    data_norm['uSN']=(data_norm['SN']*(1-data_norm['SN'])/reps)**.5
    data_norm['uFS']=(data_norm['FS']*(1-data_norm['FS'])/(reps-data['SN']+0.))**.5
    data_norm['uFF']=(data_norm['FF']*(1-data_norm['FF'])/(reps-data['SN']+0.))**.5

    data_norm['FinalRO_SN']=data['FinalRO_SN']/(data['SN']+0.)
    data_norm['FinalRO_FS']=data['FinalRO_FS']/(data['FS']+0.)
    data_norm['FinalRO_FF']=data['FinalRO_FF']/(data['FF']+0.)
    data_norm['FinalRO_Succes']=(data['FinalRO_SN']+data['FinalRO_FS'])/(data['SN']+data['FS']+0.)
    data_norm['FinalRO_All']=(data['FinalRO_SN']+data['FinalRO_FS']+data['FinalRO_FF'])/(reps+0.)

    print data['FS']
    for i in np.arange(len(data['SN'])):
        z+=data['FinalRO_FS'][i]
        zrep+=data['FS'][i]
    print zrep
    znorm=z/(zrep+0.)
    data.close()

phasecor,uphasecor=sc.get_nuclear_ROC(phasenorm,rep,sc.get_latest_data('SSRO',date=d))
zcor, uzcor=sc.get_nuclear_ROC(znorm,zrep,sc.get_latest_data('SSRO',date=zd))
figure1=plt.figure(2)
ax=figure1.add_subplot(111)
ax.errorbar(x,phasenorm,fmt='o',yerr=1/sqrt(rep),color='Crimson')
ax.errorbar(x,phasecor,fmt='o',yerr=uphasecor,color='RoyalBlue')
phase_fit=sc.fit_sin(x,phasecor,uphasecor)
xcor=np.abs(phase_fit['params'][1])+np.abs(phase_fit['params'][2])
uxcor=np.sqrt(phase_fit['error_dict']['A']**2+phase_fit['error_dict']['a']**2)

dm,f,uf,ideal=tls.calc_fidelity_psi(tau,1-zcor,xcor,utau,uzcor,uxcor,th=th,dir=dir)
idr=tls.make_rho(ideal[0]**2,ideal[1]*ideal[0])
print idr

tls.make_hist(dm[0],np.array([[0,0],[0,0]]))
print 'Fidelity',f,'  +-',uf
print 'Ideal state:', ideal
print 'uz: ',uzcor,'  ux: ',uxcor
