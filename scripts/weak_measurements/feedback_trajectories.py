import os
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D

def calc_meas_strength(x,t_zero,t_star):
    measstren=theta(x,t_zero,t_star)/90.
    return measstren

def theta(tau,t_zero,t_star):
    return 90-2*np.arccos(sqrt(S(tau,t_zero,t_star)))*180./np.pi

def S(tau,t_zero,t_star):
    return np.exp(-(tau/t_star)**2)*np.cos(np.pi/4-(tau+t_zero)*np.pi*.002185/2.)**2

data_norm,data_corr = sc.plot_feedback('172738')
y_once=data_corr['FinalRO_SN']
uy_once=data_corr['uFinalRO_SN']

y_twice=data_corr['FinalRO_FS']
uy_twice=data_corr['uFinalRO_FF']

tau=data_corr['sweep_par']

x=calc_meas_strength(tau,12,2400)
x_name=data_corr['sweep_par_name']

figure42=plt.figure(42)
plt.clf()
plt.errorbar(x,.9-y_once,fmt='o', yerr=uy_once,label='Collapse once',color='RoyalBlue')
plt.errorbar(x,y_twice,fmt='o', yerr=uy_twice,label='Collapse twice',color='Crimson')

plt.plot(np.linspace(0,x.max(),100),(np.sin(np.pi/4.+np.linspace(0,x.max(),100)*np.pi/4.))**2,'RoyalBlue')

N=(np.sin(np.pi/4.+np.linspace(0,x.max(),100)*np.pi/4.))**4+(np.sin(np.pi/4.-np.linspace(0,x.max(),100)*np.pi/4.))**4
plt.plot(np.linspace(0,x.max(),100),((np.sin(np.pi/4.+np.linspace(0,x.max(),100)*np.pi/4.))**4)/N,'Crimson')
plt.xlabel ('Measurement strength (a.u.)', fontsize = 16)
plt.ylabel ('|<sz>|', fontsize = 16)
plt.title('Collapse for one/two ROs',fontsize=16)
plt.ylim ([0, 1])
plt.xlim ([0,1 ])
plt.legend(loc=3)
plt.show()

