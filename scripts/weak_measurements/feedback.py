import os
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
import plots
from analysis.lib.tools import weaktools as tls

from analysis.lib.math import tomography as tom

basepath=r'D:/machielblok/Desktop/PhD/QTlab/data/output'
name='feedback_state'
fs=plots.fontsize

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
def make_hist(data,title=''):

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

    ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color='RoyalBlue',alpha=0.75)

#sh()
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_zlim3d([0,1])
    ax.set_zlabel('')
    ax.set_title(title)
    plt.show()
def calc_meas_strength(x,t_zero,t_star):
    measstren=theta(x,t_zero,t_star)/90.
    return measstren

def theta(tau,t_zero,t_star):
    return 90-2*np.arccos(sqrt(S(tau,t_zero,t_star)))*180./np.pi

def S(tau,t_zero,t_star):
    return np.exp(-(tau/t_star)**2)*np.cos(np.pi/4-(tau+t_zero)*np.pi*.002185/2.)**2

#tau 50 ns, 4us RO 
date='20130401'
foldernamez_4 = '141914'
foldernamex_4 = '143456'

#tau 50 ns, 2 us RO
date='20130401'
foldernamez_2 = '140024'
foldernamex_2 = '133605'
tau=50.
utau=1.
dir='up'
th=''
#tau 50 ns, 6 us RO
date='20130401'
foldernamez_6 = '153222'
foldernamex_6 = '150233'

#tau 50 ns, 8 us RO
date='20130401'
foldernamez_8 = '180434'

#tau 50 ns, 15 us RO
date='20130401'
foldernamez_15 = '124303'
foldernamex_15 = ''

date2='20130404'
#tau 50 ns, 20 us RO
foldernamez_20 = '113803'
foldernamex_20 = ''
#tau 50 ns, 30 us RO
foldernamez_30 = '114334'
foldernamex_30 = ''
#tau 50 ns, 40 us RO
foldernamez_40 = '114951'
foldernamex_40 = ''
#tau 50 ns, 50 us RO
foldernamez_50 = '115605'
foldernamex_50 = ''

succesprob_fn=[foldernamez_2,foldernamez_4,foldernamez_6,foldernamez_8,foldernamez_15]
succesprob_fn2=['111620','113803','114334','110834','114951','115605','120436','120919','121458','122045','131359']
RO_times=[2,4,6,8,15,10,20,30,6,40,50,1,2,4,8,50]

def targetstate_analysis(foldernamez,foldernamex, filename='Spin_RO',date='',RO_time=''):
    zdata_norm,zdata_corr=sc.plot_feedback(foldernamez, filename='Spin_RO',d=date)
    print RO_time
    if ((RO_time=='2us') or (RO_time=='4us') or (RO_time=='6us')): 
        xdata_norm,xdata_corr=sc.plot_feedback(foldernamex, filename='Spin_RO',d=date)
        SNfit=sc.fit_sin(xdata_norm['sweep_par'],xdata_corr['FinalRO_SN'],xdata_corr['uFinalRO_SN'])
        FSfit=sc.fit_sin(xdata_norm['sweep_par'],xdata_corr['FinalRO_FS'],xdata_corr['uFinalRO_FS'])
        allfit=sc.fit_sin(xdata_norm['sweep_par'],xdata_corr['FinalRO_All'],xdata_corr['uFinalRO_All'])
        N=(xdata_norm['SN'][2]+xdata_norm['FS'][2])
        xsucces= (xdata_norm['SN'][2]*(abs(SNfit['params'][1])*2)+
			xdata_norm['FS'][2]*(abs(FSfit['params'][1])*2))/N
	uxsucces=np.sqrt(((xdata_norm['SN'][2]*abs(SNfit['error_dict']['a'])*2)/N)**2+
	                 ((xdata_norm['FS'][2]*abs(FSfit['error_dict']['a'])*2)/N)**2)
        Sx=[abs(FSfit['params'][1])*2,abs(SNfit['params'][1])*2,xsucces]
        uSx=[FSfit['error_dict']['a']*2,SNfit['error_dict']['a']*2,uxsucces,0]
	i=2
	Sz=[zdata_corr['FinalRO_FS'][i],zdata_corr['FinalRO_SN'][i],zdata_corr['FinalRO_Succes'][i]]
    else:
	Sx=[0,0,0]
	uSx=[0,0,0]
	i=1
	Sz=[1-zdata_corr['FinalRO_FS'][i],1-zdata_corr['FinalRO_SN'][i],1-zdata_corr['FinalRO_Succes'][i]]
    if (RO_time=='15us'):	
	i=2
	Sz=[zdata_corr['FinalRO_FS'][i],zdata_corr['FinalRO_SN'][i],zdata_corr['FinalRO_Succes'][i]]


    

    Sy=[0,0,0,0]
    fdata={}
    fdata['zdata_norm']=zdata_norm
    fdata['zdata_corr']=zdata_corr
    fdata['res_FS']=[2*(1-Sz[0])-1,Sx[0],Sy[0],0]
    fdata['res_SN']=[2*(1-Sz[1])-1,Sx[1],Sy[1],0]
    fdata['res_Succes']=[2*(1-Sz[2])-1,Sx[2],Sy[2],0]
    fdata['ures_FS']=[zdata_corr['uFinalRO_FS'][i]*2,uSx[0],0,0]
    fdata['ures_SN']=[zdata_corr['uFinalRO_SN'][i]*2,uSx[1],0,0]
    fdata['ures_Succes']=[zdata_corr['uFinalRO_Succes'][i]*2,uSx[2],0,0]
    #fdata['SNfit']=SNfit
    #fdata['FSfit']=FSfit
    #fdata['allfit']=allfit
    meas_strength = calc_meas_strength(50,12,1400)
    fdata['meas_strength'] = meas_strength
    fdata['res_ideal']=[np.sin(meas_strength*np.pi/2.),np.cos(meas_strength*np.pi/2.),0,0]
    fdata['dm'],fdata['f'],fdata['uf'],fdata['ideal']=tls.calc_fidelity_psi(tau,(fdata['res_Succes'][0]+1)/2.,fdata['res_Succes'][1]/2.+0.5,utau,fdata['ures_Succes'][0]/2.,fdata['ures_Succes'][1],th=th,dir=dir)
    tls.make_hist(fdata['dm'][0],np.array([[0,0],[0,0]]))

    #print 'Fidelity',f,'  +-',uf
    #print 'Ideal state:', ideal
    #print 'uz: ',uzcor,'  ux: ',uxcor
    np.savez(os.path.join(basepath,name+RO_time),**fdata)

def targetstate_plot(RO_time):
    d=np.load(os.path.join(basepath,name+RO_time)+'.npz')
    fig=plt.figure(figsize=[1.7,1])
    gcf().subplots_adjust(bottom=0.15,top=0.975,left=0.15,right=0.975)
    ax=fig.add_subplot(111)
    width=0.2
    x=np.arange(4)

    barall=ax.bar(x, d['res_FS'],width, yerr=d['ures_FS'], color='RoyalBlue',ecolor='RoyalBlue')
    barfee=ax.bar(x+width, d['res_SN'], width, yerr=d['ures_SN'],color='Crimson',ecolor='Crimson')
    barher=ax.bar(x+2*width, d['res_Succes'], width, yerr=d['ures_Succes'],color='DarkGreen',ecolor='DarkGreen')
    baride=ax.bar(x+3*width, d['res_ideal'], width, color='White')
    ax.plot([0,3],[0,0],color='Black')
   
    #ax.set_title('Final State (tau=50ns,RO 2 us)')
    ax.set_xticks(x+width*2)
    ax.set_ylim([0,1])
    yticks=[0,0.5,1]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks,fontsize=fs)
    ax.set_xticklabels( ('<Sz>','<Sx>', '<Sy>'),fontsize=fs )
    ax.legend( (barall[0], barfee[0],barher[0],baride[0]), ('2 measurements', '1 measurement','All Succes','Target'), prop={'size':5})
    #plt.show()
    fig.savefig(os.path.join(basepath,name+'.pdf'),format='pdf')


def succesprob_analysis(folders,folderstwo,RO_times,date,date2):
    s_heralded=[]
    s_feedback=[]
    us_heralded=[]
    us_feedback=[]
    for fn in folders:
        zdata_norm,zdata_corr=sc.plot_feedback(fn, filename='Spin_RO',d=date)
        s_heralded.append(zdata_norm['SN'][2])
	us_heralded.append(zdata_norm['uSN'][2])
        s_feedback.append(zdata_norm['SN'][2]+(1-zdata_norm['SN'][2])*zdata_norm['FS'][2])
	us_feedback.append(np.sqrt(((1-zdata_norm['FS'][2])*zdata_norm['uSN'][2])**2+((1-zdata_norm['SN'][2])*zdata_norm['uFS'][2])**2))
    for fn in folderstwo:
        zdata_norm,zdata_corr=sc.plot_feedback(fn, filename='Spin_RO',d=date2)
        s_heralded.append(zdata_norm['SN'][1])
	us_heralded.append(zdata_norm['uSN'][1])
        s_feedback.append(zdata_norm['SN'][1]+(1-zdata_norm['SN'][1])*zdata_norm['FS'][1])
        us_feedback.append(np.sqrt(((1-zdata_norm['FS'][2])*zdata_norm['uSN'][2])**2+((1-zdata_norm['SN'][2])*zdata_norm['uFS'][2])**2))
    sdata={}
    sdata['s_heralded']=s_heralded
    sdata['us_heralded']=us_heralded
    sdata['s_feedback']=s_feedback
    sdata['us_feedback']=us_feedback
    sdata['RO_times']=RO_times
    np.savez(os.path.join(basepath,'Psucces'),**sdata)
def fit_exp(d,yname,xname='RO_times'):
    offset_guess=d[yname].max()
    amp_guess=offset_guess
    decay_guess=10
    feedback_fit_result = fit.fit1d(d[xname], d[yname], common.fit_exp_decay_with_offset, 
		    offset_guess, amp_guess,decay_guess,
		    do_print = False , ret = True)
    x=np.linspace(d[xname].min(),d[xname].max(),501)
    y=feedback_fit_result['fitfunc'](x)
    return x,y
def succesprob_plot():
    d=np.load(os.path.join(basepath,'Psucces')+'.npz')

    #fit succes probabilities
    feedback_fit_x,feedback_fit_y = fit_exp(d,'s_feedback','RO_times')
    heralded_fit_x,heralded_fit_y = fit_exp(d,'s_heralded','RO_times')
    fig=plt.figure(figsize=[1.7,2])

    ax=fig.add_subplot(212)
    gcf().subplots_adjust(bottom=0.15,top=0.975,left=0.15,right=0.975)

    ax.plot(feedback_fit_x,feedback_fit_y,'-',color='Crimson')
    ax.plot(heralded_fit_x,heralded_fit_y,'-',color='RoyalBlue')
    ax.errorbar(d['RO_times'],d['s_feedback'],yerr=d['us_feedback'],marker='o',
		    mfc='Crimson', mec='Crimson', ecolor='Crimson',linestyle='None',
		    label='2 measurements')
    ax.errorbar(d['RO_times'],d['s_heralded'],yerr=d['us_heralded'],marker='o',
		    mfc='RoyalBlue',mec='RoyalBlue',ecolor='RoyalBlue',linestyle='None',
		    label='1 measurement')

    yticks=[0,0.5,1]
    xticks=[0,25,50]
    ax.set_ylim([0,1.05])
    ax.set_xlim([0,50])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks,fontsize=fs)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks,fontsize=fs)
    ax.set_xlabel('RO time [us]',fontsize=fs)
    ax.set_ylabel('P(Succes)',fontsize=fs)
    ax.legend(loc=1, prop={'size':5})
    
   
   #collect fidelity data
    d2=np.load(os.path.join(basepath,name+'2us')+'.npz')
    d4=np.load(os.path.join(basepath,name+'4us')+'.npz')
    d6=np.load(os.path.join(basepath,name+'6us')+'.npz')
    d15=np.load(os.path.join(basepath,name+'15us')+'.npz')
    d20=np.load(os.path.join(basepath,name+'20us')+'.npz')
    d30=np.load(os.path.join(basepath,name+'30us')+'.npz')
    d40=np.load(os.path.join(basepath,name+'40us')+'.npz')
    d50=np.load(os.path.join(basepath,name+'50us')+'.npz')
    fRO=[2,4,6,15,20,30,40,50]
    f=[d2['f'],d4['f'],d6['f'],d15['f'],d20['f'],d30['f'],d40['f'],d50['f']]
    uf=[d2['uf'],d4['uf'],d6['uf'],d15['uf'],d20['uf'],d30['uf'],d40['uf'],d50['uf']]

    z=[2*float(np.real(d2['dm'][0][0,0]))-1,2*float(np.real(d4['dm'][0][0,0]))-1,
       2*float(np.real(d6['dm'][0][0,0]))-1,2*float(np.real(d15['dm'][0][0,0]))-1,
       2*float(np.real(d20['dm'][0][0,0]))-1,2*float(np.real(d30['dm'][0][0,0]))-1,
       2*float(np.real(d40['dm'][0][0,0]))-1,2*float(np.real(d50['dm'][0][0,0]))-1]
    uz=[2*float(np.real(d2['dm'][1][0,0])),2*float(np.real(d4['dm'][1][0,0])),
        2*float(np.real(d6['dm'][1][0,0])),2*float(np.real(d15['dm'][1][0,0])),
	2*float(np.real(d20['dm'][1][0,0])),2*float(np.real(d30['dm'][1][0,0])),
        2*float(np.real(d40['dm'][1][0,0])),2*float(np.real(d50['dm'][1][0,0]))]
    
    x=[2*float(np.real(d2['dm'][0][1,0])),2*float(np.real(d4['dm'][0][1,0])),
       2*float(np.real(d6['dm'][0][1,0])),2*float(np.real(d15['dm'][0][1,0])),
       2*float(np.real(d20['dm'][0][1,0])),2*float(np.real(d30['dm'][0][1,0])),
       2*float(np.real(d40['dm'][0][1,0])),2*float(np.real(d50['dm'][0][1,0]))]
    ux=[2*float(np.real(d2['dm'][1][1,0])),2*float(np.real(d4['dm'][1][1,0])),
        2*float(np.real(d6['dm'][1][1,0])),2*float(np.real(d15['dm'][1][1,0])),
	2*float(np.real(d20['dm'][1][1,0])),2*float(np.real(d30['dm'][1][1,0])),
        2*float(np.real(d40['dm'][1][1,0])),2*float(np.real(d50['dm'][1][1,0]))]
    print x
    #fit fidelity
    fid={}
    fid['f']=np.array(f)
    fid['RO_times']=np.array(fRO)
    fid_fit_x,fid_fit_y = fit_exp(fid,'f','RO_times')

    sx={}
    sx['x']=np.array(x)
    sx['RO_times']=np.array(fRO)
    sx_fit_x,sx_fit_y = fit_exp(sx,'x','RO_times')
    #plot fidelity
    bx=fig.add_subplot(211)
    bx.plot(fid_fit_x,fid_fit_y,'-',color='DarkGreen')
    bx.plot(sx_fit_x,sx_fit_y,'-',color='LimeGreen')
    bx.errorbar(fRO,f,yerr=uf,marker='o',
		    mfc='DarkGreen',mec='DarkGreen',ecolor='DarkGreen',linestyle='None',
		    label='Fidelity')
    bx.errorbar(fRO,x,yerr=ux,marker='s',
		    mfc='LimeGreen',mec='LimeGreen',ecolor='LimeGreen',linestyle='None',
		    label='<Sx>')
    bx.errorbar(fRO,z,yerr=uz,marker='^',
		    mfc='LightGreen',mec='LightGreen',ecolor='LightGreen',linestyle='None',
		    label='<Sz>')
    zideal=(d50['ideal'][0]**2)*2-1
    xideal=(d50['ideal'][0]*d50['ideal'][1])*2
    bx.plot([0,50],[zideal,zideal],'--',ms=0.5,
		    color='LightGreen')
    bx.plot([0,50],[xideal,xideal],'--',ms=0.5,
		    color='LimeGreen')
    #ax.set_ylabel('P(Succes)',fontsize=fs)
    #ax.set_xlabel('RO duration (us)',fontsize=fs)
    #ax.set_title('Probability to reach Target state')
    subplots_adjust(hspace=0,left=0.2)
    yticks=[0,0.5,1]
    xticks=[]
    bx.set_ylim([-0.05,1])
    bx.set_xlim([0,50])
    bx.set_xticks(xticks)
    bx.set_xticklabels(xticks,fontsize=fs)
    bx.set_yticks(yticks)
    bx.set_yticklabels(yticks,fontsize=fs)
    bx.legend(loc=1, prop={'size':5})
    fig.savefig(os.path.join(basepath,'Psucces.pdf'),format='pdf')

#res_FS,res_SN,res_Succes,ures_FS,ures_SN,ures_Succes=targetstate_analysis(foldernamez_2,foldernamex_2,date=date) 
'''
RO_time='2us'
targetstate_analysis(foldernamez_2,foldernamex_2,date=date,RO_time=RO_time)

RO_time='4us'
targetstate_analysis(foldernamez_4,foldernamex_4,date=date,RO_time=RO_time)
RO_time='6us'
targetstate_analysis(foldernamez_6,foldernamex_6,date=date,RO_time=RO_time)
RO_time='15us'
targetstate_analysis(foldernamez_15,foldernamex_15,date=date,RO_time=RO_time)


RO_time='20us'
targetstate_analysis(foldernamez_20,foldernamex_20,date=date2,RO_time=RO_time)

RO_time='30us'
targetstate_analysis(foldernamez_30,foldernamex_30,date=date2,RO_time=RO_time)
RO_time='40us'
targetstate_analysis(foldernamez_40,foldernamex_40,date=date2,RO_time=RO_time)
RO_time='50us'
targetstate_analysis(foldernamez_50,foldernamex_50,date=date2,RO_time=RO_time)
'''
#targetstate_plot('2us')
#succesprob_analysis(succesprob_fn,succesprob_fn2,RO_times,date,date2)
succesprob_plot()





