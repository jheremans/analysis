
import numpy as np
import pylab as plt
from analysis.lib.spin import spin_control as sc
from analysis.lib.math import tomography as tom
import os
from mpl_toolkits.mplot3d import Axes3D

def make_rho(z,x):

    rho=np.array([[z,x],[x,1-z]])
    return rho
def make_rho_ideal(tau,input='X'):
    s=calc_meas_strength(tau,12.,2400.)
    x=np.cos(s*np.pi/2.)/2.
    z=np.sin(np.pi+s*np.pi/2.)/2.+0.5
    rho_cond=make_rho(1-z,x)
    rho_uncond=make_rho(0.5,x)
    return rho_cond, rho_uncond,s
def make_hist(data,ideal=np.array([[0.5,0.5],[0.5,0.5]]),title='',path=''):
    column_names = ['mI=0','mI=-1']
    row_names = ['mI=0','mI=-1']

    fig = plt.figure(figsize=[0.8,0.8])
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
    ax.plot([0,2],[0,0],[0.5,0.5],color='Grey',alpha=0.5)
    #ax.plot([0,0],[0,2],[0.5,0.5],color='Grey',alpha=0.5)
    #ax.plot([0,0.5],[2,2],[0.5,0.5],color='Grey',alpha=0.5)
    ax.plot([2,2],[0,2],[0.5,0.5],color='Grey',alpha=0.5)
    ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color='Crimson',alpha=1)
    ax.bar3d(xpos,ypos,zpos, dx, dy, dzideal, color='b',alpha=0)
 
    ax.w_xaxis.set_ticks([0.5,1.5])
    ax.w_yaxis.set_ticks([0.5,1.5])
    ax.w_xaxis.set_ticklabels(column_names,fontsize=4)
    ax.w_yaxis.set_ticklabels(row_names,fontsize=4)
    ax.set_zlim3d([0,1])
    ax.set_zticks([0,0.5,1])
    ax.set_zticklabels([0,0.5,1],fontsize=4)
    ax.set_zlabel('')
    ax.set_title(title)
    ax.view_init(elev=25.,azim=-40.)
    if path!='':
	print path+title+'.pdf'
	fig.savefig(os.path.join(path,title+'.pdf'),format='pdf')
    plt.show()
def calc_meas_strength(x,t_zero=12.,t_star=1400.):
    measstren=calc_theta(x,t_zero,t_star)/90.
    return measstren
def calc_theta(tau,t_zero,t_star):
    return 90-2*(np.arccos(np.sqrt(S(tau,t_zero,t_star))))*180./np.pi

def S(tau,t_zero,t_star):
    return np.exp(-(tau/t_star)**2)*np.cos(np.pi/4-(tau+t_zero)*np.pi*.002185/2.)**2

def projection(init,theta,direction):
    a=np.pi*theta/360.	
    if (direction=='up'):
        after=np.array([init[0]*(np.cos(a)+np.sin(a)),init[1]*(np.cos(a)-np.sin(a))])
    else:
	after=np.array([init[0]*(np.cos(a)-np.sin(a)),init[1]*(np.cos(a)+np.sin(a))])
    return after
def calc_fidelity_psi(tau,z,x,utau,uz,ux,y=0.5,uy=0,th='',dir='up',t_zero=12,t_star=1400):
    if (th==''):
        theta = calc_theta(tau,t_zero,t_star)
    else:
	theta=th
    print 'THETA'
    print theta
    init=np.array([1/np.sqrt(2.),1/np.sqrt(2.)])
    ideal=projection(init,theta,dir)
    dm=tom.measured_single_qubit_dm(x,y,z,ux,uy,uz)
    f,uf=tom.fidelity(ideal,dm[0],dm[1])
    return dm, f, uf, ideal
    
