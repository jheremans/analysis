import os
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib import rcParams

from mpl_toolkits.mplot3d import Axes3D
from analysis.lib.tools import weaktools as tls
from analysis.lib.math import tomography as tom

#datafolders=['1154','1252','1258','1303','1310','1347','1351','1316','1326', '1453']
#RO_time=[0,1,2,3,4,5,6,7,9, 11]
#
#
#datafolders=['164825','165126']
#date='20130111'

#datafolders=['111709','112023']
#date='20130115'
def make_rho(z,x):

    rho=np.array([[x,z],[1-z,x]])
    return rho
def make_rho_ideal(tau,input='X'):
    s=calc_meas_strength(tau,12.,2400.)
    x=np.cos(s*np.pi/2.)/2.
    z=np.sin(np.pi+s*np.pi/2.)/2.+0.5
    rho_cond=make_rho(z,x)
    rho_uncond=make_rho(0.5,x)
    return rho_cond, rho_uncond,s
def make_hist(data,ideal=np.array([[0.5,0.5],[0.5,0.5]]),title=''):
    column_names = ['','']
    row_names = ['','']

    fig = plt.figure()
    ax = Axes3D(fig)

    lx= len(data[0])            # Work out matrix dimensions
    ly= len(data[:,0])
    xpos = np.arange(0,lx,1)    # Set up a mesh of positions
    ypos = np.arange(0,ly,1)
    xpos, ypos = np.meshgrid(xpos+0.25, ypos+0.25)

    xpos = xpos.flatten()   # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx*ly)

    dx = 0.75 * np.ones_like(zpos)
    dy = dx.copy()
    dz = data.flatten()

    dzideal = ideal.flatten()

    ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color='RoyalBlue',alpha=0.65)
    ax.bar3d(xpos,ypos,zpos, dx, dy, dzideal, color='b',alpha=0)
    
    
#sh()
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_zlim3d([0,1])
    ax.set_zlabel('')
    ax.set_title(title)
    plt.show()
def calc_meas_strength(x,t_zero=12.,t_star=1400.):
    measstren=theta(x,t_zero,t_star)/90.
    return measstren

def theta(tau,t_zero,t_star):
    return 90-2*(np.arccos(sqrt(S(tau,t_zero,t_star))))*180./np.pi

def S(tau,t_zero,t_star):
    return np.exp(-(tau/t_star)**2)*np.cos(np.pi/4-(tau+t_zero)*np.pi*.002185/2.)**2



'''
corr=(.728,1-.088)
datafolders=['175338','175003']
date='20130225'
label=['mI=0','mI=+1']
meas_folder=r'D:\measuring\data'

dp=os.path.join(meas_folder, date)
dphase=[]
plt.close('all')
plt.figure(75)
j=0
for i in datafolders:    
    
    result=sc.plot_data_MBI(sc.get_latest_data(string=i,datapath=dp),fid=corr)
    plt.figure(75)
    measstrent=(90*result['x']/214 +2)/90
    t=result['x']
    result['x']=measstrent
    #result['x']=result['x']/191
    plt.errorbar(result['x'],2*(result['y']-0.5),fmt='o',yerr=result['yerr'],label=label[j])
    #plt.errorbar(result['x'],result['y'],fmt='o',yerr=result['yerr'],label=label[j])
    yerr_weak=np.sqrt(2*(1+np.cos(result['y']))/500.)/sin(result['x']*90*pi/180)
    #plt.errorbar(result['x'],2*(result['y']-0.5)/np.sin(result['x']*90*pi/180),fmt='o',yerr=yerr_weak,label=label[j])
    j=j+1
plt.plot(np.linspace(0,1,100),np.exp(-(np.linspace(0,0,100)/1300)**2)*np.sin(np.linspace(0,1,100)*np.pi/2),'r-')
plt.plot(np.linspace(0,1,100),np.exp(-(np.linspace(0,0,100)/1300)**2)*np.sin(np.pi+np.linspace(0,1,100)*np.pi/2),'r-')

plt.xlabel (' Measurement Strength (AU)', fontsize = 16)
plt.ylabel ('< $\sigma_{Z}$ > $_{nitrogen}$', fontsize = 16)   
plt.ylabel ('| $\uparrow$ > $_{nitrogen}$', fontsize = 16)   
plt.ylim ([-1, 1])
plt.xlim ([-0.01, 1])
plt.legend(loc=2)
plt.show()

datafolders=['_25ns','_125ns','_225ns']
date='20130118'
dp=os.path.join(meas_folder, date)
label=['mI=0','mI=+1']
meas_folder=r'D:\measuring\data'
'''

date='20130314'
#-y
state='Y'
zfolder='004118'
yfolder='233139'
xfolder='000101'

#state='-1'
#zfolder='020210'
#yfolder='010717'
#xfolder='013224'

#0
#zfolder='032201'
#yfolder='022809'
#xfolder='025425'

#x
#state='X'
#zfolder='050829'
#yfolder='040013'
#xfolder='044207'

#-x
#state='-X'
#zfolder='062515'
#yfolder='053455'
#xfolder='060013'

#y
#state='-Y'
#zfolder='074639'
#yfolder='065137'
#xfolder='071845'

#date='20130305'
#zfolder='005909'
#yfolder='232520'
#xfolder='024929'
#date='20130313'
# min y:
#zfolder='214008'
#yfolder='224959'
#xfolder='224131'

# mI-1:
#zfolder='193905'
#yfolder='185211'
#xfolder='191458'

# mI0:
#zfolder='210551'
#yfolder='201904'
#xfolder='204243'

# x:
#zfolder='224310'
#yfolder='213105'
#xfolder='221323'

#dp=os.path.join(meas_folder, date)
result_zmeas=sc.analyse_plot_results_vs_sweepparam(zfolder,yname='P(mI=0)',Nuclcor=True,dataname='Spin_RO',title='strong_meas_res_Z',d=date,save=False)
result_zcond=sc.analyse_weakcond_vs_sweepparam(zfolder,yname='P(mI=0)',Nuclcor=True,dataname='cond_Spin_RO',title='',d=date,save=False)

result_ymeas=sc.analyse_plot_results_vs_sweepparam(yfolder,yname='P(mI=0)',Nuclcor=True,dataname='Spin_RO',title='strong_meas_res_y',d=date,save=False)
result_ycond=sc.analyse_weakcond_vs_sweepparam(yfolder,yname='P(mI=0)',Nuclcor=True,dataname='cond_Spin_RO',title='',d=date,save=False)

result_xmeas=sc.analyse_plot_results_vs_sweepparam(xfolder,yname='P(mI=0)',Nuclcor=True,dataname='Spin_RO',title='strong_meas_res_x',d=date,save=False)
result_xcond=sc.analyse_weakcond_vs_sweepparam(xfolder,yname='P(mI=0)',Nuclcor=True,dataname='cond_Spin_RO',title='',d=date,save=False)

zmeas_x=calc_meas_strength(result_zmeas['x'],12,3400)
xmeas_x=calc_meas_strength(result_xmeas['x'],12,3400)

zmeas_x=calc_meas_strength(result_zmeas['x'],12,2000)
xmeas_x=calc_meas_strength(result_xmeas['x'],12,2000)

def plot_backaction():
    figure42=plt.figure(42,figsize=[0.8,0.8])
    #plt.figure(42)
    print 'plotting backaction'
    plt.clf()

    plt.errorbar(zmeas_x,result_zcond['y_cond'],fmt='o', yerr=result_zcond['uy_cond'],label='|mI=-1>',color='RoyalBlue')
    plt.errorbar(zmeas_x,1-result_zcond['y_cond'],fmt='o', yerr=result_zcond['uy_cond'],label='|mI=0>',color='Crimson')
    plt.errorbar(xmeas_x,result_xcond['y_cond']-0.5,fmt='o',yerr=result_xcond['uy_cond'],label='|y> (Im)',color='LimeGreen')
    plt.errorbar(xmeas_x,result_ycond['y_cond']-.5,fmt='o',yerr=result_ycond['uy_cond'],label='|x> (Re)',color='DarkOrange')
    if state =='Y':
    	plt.plot(np.linspace(0,1,100),(np.cos(np.linspace(0,1,100)*np.pi/2.)/2.),ls='-',color='LimeGreen')
    	plt.plot(np.linspace(0,1,100),(np.sin(np.pi+np.linspace(0,1,100)*np.pi/2.)/2.+0.5),ls='-',color='Crimson')
    	plt.plot(np.linspace(0,1,100),(np.sin(np.linspace(0,1,100)*np.pi/2.)/2.+0.5),'g-',color='RoyalBlue')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),'g-',color='DarkOrange')
	plt.text(0.2, 0.92, r'initial state |y>')
    	plt.ylim ([-0.5, 1])
    	plt.yticks([-0.5,0,0.5,1])

    if state =='-Y':
    	plt.plot(np.linspace(0,1,100),(np.cos(np.pi+np.linspace(0,1,100)*np.pi/2.)/2.),ls='-',color='LimeGreen')
    	plt.plot(np.linspace(0,1,100),(np.sin(np.pi+np.linspace(0,1,100)*np.pi/2.)/2.+0.5),ls='-',color='Crimson')
    	plt.plot(np.linspace(0,1,100),(np.sin(np.linspace(0,1,100)*np.pi/2.)/2.+0.5),'g-',color='RoyalBlue')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),'g-',color='DarkOrange')
	plt.text(0.2, 0.92, r'initial state |-y>')
    	plt.ylim ([-0.5, 1])
    	plt.yticks([-0.5,0,0.5,1])
	
    if state =='X':
    	plt.plot(np.linspace(0,1,100),(np.cos(np.linspace(0,1,100)*np.pi/2.)/2.),ls='-',color='DarkOrange')
    	plt.plot(np.linspace(0,1,100),(np.sin(np.pi+np.linspace(0,1,100)*np.pi/2.)/2.+0.5),ls='-',color='Crimson')
    	plt.plot(np.linspace(0,1,100),(np.sin(np.linspace(0,1,100)*np.pi/2.)/2.+0.5),'g-',color='RoyalBlue')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),'g-',color='LimeGreen')
	plt.text(0.2, 0.92, r'initial state |x>')
    	plt.ylim ([-0.5, 1])
    	plt.yticks([-0.5,0,0.5,1])
    if state =='-X':
    	plt.plot(np.linspace(0,1,100),(np.cos(np.pi+np.linspace(0,1,100)*np.pi/2.)/2.),ls='-',color='DarkOrange')
    	plt.plot(np.linspace(0,1,100),(np.sin(np.pi+np.linspace(0,1,100)*np.pi/2.)/2.+0.5),ls='-',color='Crimson')
    	plt.plot(np.linspace(0,1,100),(np.sin(np.linspace(0,1,100)*np.pi/2.)/2.+0.5),'g-',color='RoyalBlue')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),'g-',color='LimeGreen')
	plt.text(0.2, 0.92, r'initial state |-x>')
    	plt.ylim ([-0.5, 1])
    	plt.yticks([-0.5,0,0.5,1])
    if state =='-1':
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),ls='-',color='LimeGreen')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),ls='-',color='Crimson')
    	plt.plot(np.linspace(0,1,100),1+0*np.linspace(0,1,100),'g-',color='RoyalBlue')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),'g-',color='DarkOrange')
	plt.text(0.2, 1.2, r'initial state |mI=-1>')
	plt.ylim ([-0.2, 1.3])
    	plt.yticks([0,0.5,1])

    plt.xlabel (' Measurement Strength (a.u.)', fontsize = 16)
    plt.ylabel ('Density Matrix Element', fontsize = 16)
    plt.title('Backaction conditioned on weak measurement |mI=-1>',fontsize=16)

    plt.xlim ([0, 1])
    plt.xticks([0,0.5,1])
    plt.legend(loc=2,prop={'size':10})
    plt.show()

    figure43=plt.figure(43,figsize=[0.8,0.8])
    #plt.figure(42)
    print 'plotting backaction'
    plt.clf()

    plt.errorbar(zmeas_x,result_zmeas['y'],fmt='o', yerr=result_zmeas['uy'],label='|mI=-1>',color='RoyalBlue')
    plt.errorbar(zmeas_x,1-result_zmeas['y'],fmt='o', yerr=result_zmeas['uy'],label='|mI=0>',color='Crimson')
    plt.errorbar(xmeas_x,result_xmeas['y']-0.5,fmt='o',yerr=result_xmeas['uy'],label='|y> (Im)',color='LimeGreen')
    plt.errorbar(xmeas_x,result_ymeas['y']-0.5,fmt='o',yerr=result_ymeas['uy'],label='|x> (Re)',color='DarkOrange')
    
    if state =='Y':
    	plt.plot(np.linspace(0,1,100),(np.cos(np.linspace(0,1,100)*np.pi/2.)/2.),color='LimeGreen',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100)+0.5,color='RoyalBlue',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100)+0.5,color='Crimson',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),'g-',color='DarkOrange')
    	plt.text(0.2, 0.92, r'initial state |y>')
	plt.ylim ([-0.5, 1])
	plt.yticks([-0.5,0,0.5,1])

    if state =='-Y':
    	plt.plot(np.linspace(0,1,100),(np.cos(np.pi+np.linspace(0,1,100)*np.pi/2.)/2.),color='LimeGreen',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100)+0.5,color='RoyalBlue',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100)+0.5,color='Crimson',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),'g-',color='DarkOrange')
    	plt.text(0.2, 0.92, r'initial state |-y>')
	plt.ylim ([-0.5, 1])
	plt.yticks([-0.5,0,0.5,1])

    if state =='X':
    	plt.plot(np.linspace(0,1,100),(np.cos(np.linspace(0,1,100)*np.pi/2.)/2.),color='DarkOrange',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100)+0.5,color='RoyalBlue',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100)+0.5,color='Crimson',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),'g-',color='LimeGreen')
    	plt.text(0.2, 0.92, r'initial state |x>')
	plt.ylim ([-0.5, 1])
	plt.yticks([-0.5,0,0.5,1])
    if state =='-X':
    	plt.plot(np.linspace(0,1,100),(np.cos(np.pi+np.linspace(0,1,100)*np.pi/2.)/2.),color='DarkOrange',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100)+0.5,color='RoyalBlue',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100)+0.5,color='Crimson',ls='-')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),'g-',color='LimeGreen')
    	plt.text(0.2, 0.92, r'initial state |-x>')
	plt.ylim ([-0.5, 1])
	plt.yticks([-0.5,0,0.5,1])
    if state =='-1':
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),ls='-',color='LimeGreen')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),ls='-',color='Crimson')
    	plt.plot(np.linspace(0,1,100),1+0*np.linspace(0,1,100),'g-',color='RoyalBlue')
    	plt.plot(np.linspace(0,1,100),0*np.linspace(0,1,100),'g-',color='DarkOrange')
	plt.text(0.2, 1.2, r'initial state |mI=-1>')
	plt.ylim ([-0.2, 1.3])
    	plt.yticks([0,0.5,1])
	

    plt.xlabel (' Measurement Strength (a.u.)', fontsize = 16)
    plt.ylabel ('Density Matrix Element', fontsize = 16)
    plt.title('Backaction (unconditioned)',fontsize=16)
    
    
    plt.xlim ([0, 1])
    plt.xticks([0,0.5,1])
    plt.legend(loc=2,prop={'size':10})
    plt.show()

#plot_backaction()    

dir='up'
index=0
postselect='cond'

tau=result_xmeas['x'][index]
t=str(tau)+'ns_'+postselect
utau=1
th=''
if postselect!='cond':
    zcor=result_zmeas['y'][index]
    xcor=result_xmeas['y'][index]
    ycor=result_ymeas['y'][index]

    uzcor=result_zmeas['uy'][index]
    uxcor=result_xmeas['uy'][index]
    uycor=result_ymeas['uy'][index]
else:
    zcor=result_zcond['y_cond'][index]
    xcor=result_xcond['y_cond'][index]
    ycor=result_ycond['y_cond'][index]

    uzcor=result_zcond['uy_cond'][index]
    uxcor=result_xcond['uy_cond'][index]
    uycor=result_ycond['uy_cond'][index]
dm,f,uf,ideal=tls.calc_fidelity_psi(tau,zcor,xcor,utau,uzcor,uxcor,y=ycor,uy=uycor,th=th,dir=dir)
idr=tls.make_rho(ideal[0]**2,ideal[1]*ideal[0])
print idr

tls.make_hist(dm[0],np.array([[0,0],[0,0]]),path=r'D:\machielblok\Desktop\PhD\QTlab\data\output\backaction',title=t)
print 'Fidelity',f,'  +-',uf
print 'Ideal state:', ideal
print 'uz: ',uzcor,'  ux: ',uxcor
'''

result_25ns=sc.plot_rabi(sc.get_latest_data(datafolders[0],datapath=dp),fid=corr)
result_125ns=sc.plot_rabi(sc.get_latest_data(datafolders[1],datapath=dp),fid=corr)
result_225ns=sc.plot_rabi(sc.get_latest_data(datafolders[2],datapath=dp),fid=corr)

date='20130116'
dp=os.path.join(meas_folder, date)
result_zmeas=sc.plot_data_MBI(sc.get_latest_data('180000',datapath=dp),fid=corr)
rho_25 = make_rho(result_zmeas['y'][3],abs(result_25ns[0]['params'][1]))
rho_125 =make_rho(result_zmeas['y'][11],abs(result_125ns[0]['params'][1]))
rho_225 =make_rho(result_zmeas['y'][-1],abs(result_225ns[0]['params'][1]))

rho_25_cond = make_rho(result_zmeas['y_cond'][3],abs(result_25ns[0]['params'][1]))
rho_125_cond =make_rho(result_zmeas['y_cond'][11],abs(result_125ns[0]['params'][1]))
rho_225_cond =make_rho(result_zmeas['y_cond'][-1],abs(result_225ns[0]['params'][1]))

make_hist(rho_25_cond,title='tau = 25 ns')
make_hist(rho_125_cond,title='tau = 125 ns')
make_hist(rho_225_cond,title='tau = 225 ns')

make_hist(rho_25,title='tau = 25 ns')
make_hist(rho_125,title='tau = 125 ns')
make_hist(rho_225,title='tau = 225 ns')
'''
