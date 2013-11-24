import os
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
import plots

fs=plots.fontsize

path=r'D:\machielblok\Desktop\PhD\QTlab\data\output\Weak_value'
title='weak_value_basisrot'

data_1ns=np.load(r'D:\machielblok\Desktop\PhD\QTlab\data\output\Weak_Value\data_1ns.npz')
data_30ns=np.load(r'D:\machielblok\Desktop\PhD\QTlab\data\output\Weak_Value\data_30ns.npz')
data_100ns=np.load(r'D:\machielblok\Desktop\PhD\QTlab\data\output\Weak_Value\data_100ns.npz')
data_sweep_tau=np.load(r'D:\machielblok\Desktop\PhD\QTlab\data\output\Weak_value\weak_value_sweepTau.npz')

fig=plt.figure(figsize=[2.5,1.5])
fig.clf()
gcf().subplots_adjust(bottom=0.2,left=0.1,top=0.96,right=0.96)

#sweep rot angle different tau
ax=fig.add_subplot(111)
print data_1ns
ax.errorbar(data_1ns['exp_angle'], data_1ns['exp_weakV'], yerr=data_1ns['err_weakV'],marker='.',mec='RoyalBlue',ecolor='RoyalBlue',mfc='None',elinewidth=0.5,ms=5,label=r'$\theta$ = $5^{\circ}$',linestyle='None')
ax.errorbar(data_1ns['model_angle'], data_1ns['model_weakV'],fmt='-',color='RoyalBlue')

ax.errorbar(data_30ns['exp_angle'], data_30ns['exp_weakV'], yerr=data_30ns['err_weakV'], marker='.',mec='Crimson',mfc='None',ecolor='Crimson',label=r'$\theta$ = $16 ^{\circ}$',ms=5,elinewidth=0.5,linestyle='None')
ax.errorbar(data_30ns['model_angle'], data_30ns['model_weakV'], fmt='-',color='Crimson')

ax.errorbar(data_100ns['exp_angle'], data_100ns['exp_weakV'], yerr=data_100ns['err_weakV'], marker='.',mec='DarkGreen',mfc='None',label=r'$\theta$ = $45^{\circ}$',ecolor='DarkGreen',elinewidth=0.5,ms=5,linestyle='None')
ax.errorbar(data_100ns['model_angle'], data_100ns['model_weakV'], fmt='-',color='DarkGreen')

#ax.plot([0,200],[1,1],color='Grey',label='Maximum value of <Sz>')
#ax.plot([0,200],[-1,-1],color='Grey')

#ax.set_title('Weak Value vs Basis rot length')
ax.set_xticks([0,45,90,135,180])
ax.set_xticklabels([0,r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$',r'$\frac{3 \pi}{4}$',r'$\pi$'],fontsize=4)
#ax.set_xticklabels([0,45,90,135,180],fontsize=4)<F5>
ax.set_yticks([-5,0,5,10])
ax.set_yticklabels([-5,0,5,10],fontsize=fs)
ax.set_ylim([-6,12])
ax.set_xlim([0,180])
ax.set_xlabel('Rotation angle $\phi$ (rad)',fontsize=fs)
ax.set_ylabel('Weak Value',fontsize=fs)
ax.legend(loc=2, prop={'size':6})
fig.savefig(os.path.join(path,title+'.pdf'),format='pdf')

#weakvalue vs theta
title='weak_value_theta'
fig=plt.figure(figsize=[1,0.5])
fig.clf()
bx=fig.add_subplot(111)
gcf().subplots_adjust(bottom=0.25,left=0.2,top=0.95,right=0.95)
bx.errorbar(data_sweep_tau['theta'], data_sweep_tau['weakV'], yerr=data_sweep_tau['err_weakV'],marker='.',mec='RoyalBlue',ecolor='RoyalBlue',mfc='None',elinewidth=0.5,ms=5,linestyle='None')
bx.errorbar(data_sweep_tau['th_theor'], data_sweep_tau['weakV_theor'],fmt='-',color='RoyalBlue')
bx.set_yticks([0,5,10])
bx.set_yticklabels([0,5,10],fontsize=4)
bx.set_ylim([0,12])
bx.set_xticks([0,1.57/2.,1.57])
bx.set_xticklabels([0,45,90],fontsize=4)
bx.set_xlim([0,1.57])
bx.set_xlabel(r'Measurement strength $\theta$ (degree)',fontsize=4)
bx.set_ylabel('Weak Value',fontsize=4)
fig.savefig(os.path.join(path,title+'.pdf'),format='pdf')



