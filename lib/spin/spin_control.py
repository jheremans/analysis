import os, time
import numpy as np
from numpy import *
import sympy
import matplotlib
from matplotlib import pyplot as plt
from analysis.lib.fitting import fit, rabi, common, esr, ramsey, SE
from analysis.lib.tools import plot
from analysis.lib.m2.ssro import ssro, readout_correction
from analysis.lib.ssro import ssro_adwin as ssrosc
from analysis.lib.math import error
from analysis.lib.spin import Nspin_twolevel_correction as Ncor
from measurement.lib.config import experiment_lt2 as explt2


def num2str(num, precision): 
    return "%0.*f" % (precision, num)

def dec2str(num): 
    return "%d" % (num)

def find_nearest(array,value):
    idx=(abs(array-value)).argmin()
    return idx

def get_latest_data(string = 'ADwin_SSRO', datapath = '',date=''):
    meas_folder = r'D:\measuring\data'
    if date=='':
        currdate = time.strftime('%Y%m%d')
    else:
        currdate=date
    
    if datapath == '':
        df = os.path.join(meas_folder, currdate)
    else:
        df = datapath
    
    right_dirs = list()

    if os.path.isdir(df):
        for k in os.listdir(df):
            if string in k:
                right_dirs.append(k)
        
        if len(right_dirs) > 0:
            latest_dir = os.path.join(df,right_dirs[len(right_dirs)-1])
        else:
            print 'No measurements containing %s in %s'%(string, df)
            latest_dir=''
        
        print '\nAnalyzing data in %s'%latest_dir

    else:
        print 'Folder %s does not exist'%df
        latest_dir = False

    return latest_dir

def SSRO_correct(SSRO_meas, F0, F1, F0err = 0.01, F1err = 0.01):
    w1 = SSRO_meas[:2].copy()
    w1 = w1[::-1]
    
    w1total = np.sum(w1)
    w1 /= float(w1total)

    norm = 1./(F0+F1-1.)
    w1staterr = np.sqrt(w1*(1.-w1)/float(w1total))
    
    w1p0 = w1[1]
    w1p1 = w1[0]

    w1ms0err = np.sqrt( 
            (norm**2 * (-w1p1 + F0*(w1p0+w1p1)) * F1err)**2 +\
            (norm**2 * (w1p0 - F1*(w1p0+w1p1)) * F0err)**2 +\
            (norm * (F1 - 1) * w1staterr[1])**2 +\
            (norm * F1 * w1staterr[0])**2
            )
    w1ms1err = w1ms0err
    w1err = np.array([w1ms0err, w1ms1err])


    corrmat = 1./(F0*F1 - (1.-F1)*(1.-F0)) * np.array([[F1, F1-1.],[F0-1., F0]])
    return np.dot(corrmat, w1.reshape(2,1)).reshape(-1), \
            w1err

def get_MBI_readout_result(datapath,fname,key=''):

    
    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    f=fname
    for k in files:
        if 'statics_and_parameters.npz' in k:
            stats_params_file = k
        if f in k:
            spin_ro_file = k

    e = np.load(datapath+'\\'+stats_params_file)
    par = e['sweep_par']
    noof_datapoints = len(par)
    noof_reps = e['completed_repetitions']
    if 'sweep_par_name' in e.keys():
        x_name=e['sweep_par_name']
    else:
        x_name='Total free evolution time [ns]'

    e.close()
    print fname

    
    ###########################################
    ######## SPIN RO normalization  ###########
    ###########################################
    print datapath
    print spin_ro_file
    f = np.load(datapath+'\\'+spin_ro_file)
  
    if key=='':
        key=fname


    SSRO_readout = f[key]
    print SSRO_readout
    
    ret_dict={'x':par,'y':SSRO_readout,'xname':x_name,'reps': noof_reps}
    
    return ret_dict


def get_electron_ROC(normalized_ssro,reps,ssro_calib_folder=get_latest_data('SSRO'),RO_dur=''):
    p0 = np.zeros(normalized_ssro.shape)
    u_p0 = np.zeros(normalized_ssro.shape)
    if RO_dur=='':
        ro_duration = explt2.ssroprotocol['RO_duration']
    else:
        ro_duration=RO_dur
    roc = error.SingleQubitROC()
    [[roc.F0,roc.F1],[roc.uF0,roc.uF1]]=ssrosc.get_readout_corr(ro_duration,ssro_calib_folder)
    u_normalized_ssro = (normalized_ssro*(1.-normalized_ssro)/reps)**0.5
    p0, u_p0 = roc.num_eval(normalized_ssro,u_normalized_ssro)
    print p0
    print u_p0
    return p0,u_p0

def get_nuclear_ROC(normalized_ssro,reps,ssro_calib_folder=get_latest_data('SSRO')):
    p0 = np.zeros(normalized_ssro.shape)
    u_p0 = np.zeros(normalized_ssro.shape)
    
    ro_duration = explt2.ssroprotocol['RO_duration']
    roc = Ncor.NuclearSpinROC()
    [[roc.F0_ssro,roc.F1_ssro],[roc.uF0_ssro,roc.uF1_ssro]]=ssrosc.get_readout_corr(ro_duration,ssro_calib_folder)
    print roc.F0_ssro
    u_normalized_ssro = (normalized_ssro*(1.-normalized_ssro)/reps)**0.5
    p0, u_p0 = roc.num_eval(normalized_ssro,u_normalized_ssro)
    print p0
    print u_p0
    print 'F1 and F0'
    print roc.F1
    print roc.F0
    return p0,u_p0


def plot_ssro_vs_sweep(x,y,yerr,ylim=[0,1],xname='',yname='',label='',title='',datapath='',nr='',splot=111,save=True):
    if nr=='':
            if splot==111:
                figure1 = plt.figure()
            else: 
                figure1 = plt.figure(figsize=(32,24))
    else:
            if splot==111:
                figure1 = plt.figure(nr)
            else:
                figure1 = plt.figure(nr,figsize=(32,24))
    matplotlib.rc('xtick',labelsize=14)
    matplotlib.rc('ytick',labelsize=14)
    ax=figure1.add_subplot(splot)
    ax.errorbar(x,y, fmt='o',yerr=yerr,label=label)
    
    ax.set_xlabel(xname,fontsize=14)
    ax.set_ylabel(yname,fontsize=14)
    ax.set_ylim(ylim)
    ax.set_title(title)

    if save:
        print datapath
        print title
        print'where it should be saved'
        print os.path.join(datapath,title+'.png')
        figure1.savefig(os.path.join(datapath,title+'.png'))
        figure1.savefig(os.path.join(datapath,title+'.pdf'),format='pdf')
    print 'look for file'
    print datapath
    figure1.subplots_adjust(wspace=.3,hspace=.3)
    if label!='':
        #ax.legend(loc=1)
        a=2
    return ax

def analyse_plot_results_vs_sweepparam(fname,yname,ylim=[0,1],Nuclcor=False,title='result',dataname='',datapath='',key='SSRO_counts',d='',save=True):
    if d == '':
        currdate = time.strftime('%Y%m%d')
    else:
        currdate=d

    dp=get_latest_data(fname,datapath,d)

    ssro_result=get_MBI_readout_result(dp,dataname,key=key)

    
    reps=ssro_result['reps']
    x=ssro_result['x']
    y_norm=ssro_result['y']/float(reps)
    if Nuclcor == False:
        [ssro_ro_cor,ussro_ro_cor]=get_electron_ROC(y_norm,reps,ssro_calib_folder=get_latest_data('SSRO',datapath,date=currdate))
    else:
        [ssro_ro_cor,ussro_ro_cor]=get_nuclear_ROC(y_norm,reps,ssro_calib_folder=get_latest_data('SSRO',datapath,date=currdate))
    plot_ssro_vs_sweep(x,ssro_ro_cor,ussro_ro_cor,ylim,ssro_result['xname'],yname,datapath=dp,title=currdate+'-'+fname+title,nr=1,save=save)
    
    res={'y':ssro_ro_cor,'uy':ussro_ro_cor,'x':x}

    return res
def analyse_weakcond_vs_sweepparam(fname,yname,ylim=[0,1],Nuclcor=False,title='',dataname='',datapath='',d='',save=True):
    if d == '':
        currdate = time.strftime('%Y%m%d')
    else:
        currdate=d
    dp=get_latest_data(fname,datapath,d)
    ssro_result_cond=get_MBI_readout_result(dp,dataname,key='SSRO_cond_counts')
    ssro_result_weak=get_MBI_readout_result(dp,dataname,key='SSRO_weak_counts')
    
    reps=ssro_result_weak['reps']

    x=ssro_result_weak['x']
    y_norm_weak=ssro_result_weak['y']/float(reps)
    y_norm_cond=ssro_result_cond['y']/map(float,ssro_result_weak['y'])
    [ssro_ro_cor_weak,ussro_ro_cor_weak]=get_electron_ROC(y_norm_weak,reps,ssro_calib_folder=get_latest_data('SSRO',date=currdate),RO_dur=explt2.MBIprotocol['weak_RO_duration'])
    print 'do cond'
    if Nuclcor == False:
        [ssro_ro_cor_cond,ussro_ro_cor_cond]=get_electron_ROC(y_norm_cond,np.mean(ssro_result_weak['y'])
                                                        ,ssro_calib_folder=get_latest_data('SSRO',date=currdate),RO_dur=explt2.MBIprotocol['RO_duration'])
    else:
        [ssro_ro_cor_cond,ussro_ro_cor_cond]=get_nuclear_ROC(y_norm_cond,np.mean(ssro_result_weak['y'])
                                                        ,ssro_calib_folder=get_latest_data('SSRO',date=currdate))
    plot_ssro_vs_sweep(x,ssro_ro_cor_weak,ussro_ro_cor_weak,ylim,ssro_result_weak['xname'],yname,
            currdate+'-'+fname+'weak msmsnt r175003esult',dp,nr=1,save=save)
    plot_ssro_vs_sweep(x,ssro_ro_cor_cond,ussro_ro_cor_cond,ylim,ssro_result_cond['xname'],yname,
            currdate+'-'+fname+'cond msmsnt result',dp,nr=2,save=save)
    res={'y_cond':ssro_ro_cor_cond,'y_weak':ssro_ro_cor_weak,'uy_cond':ussro_ro_cor_cond,'uy_weak':ussro_ro_cor_weak,'x':x}

    return res
def analyse_correlations_weakstrong(fname,yname,Nuclcor=False,title='',dataname='',datapath='',d='',save=True):
    if d == '':
        currdate = time.strftime('%Y%m%d')
    else:
        currdate=d

    dp=get_latest_data(fname,datapath,d)
    ssro_result_00=get_MBI_readout_result(dp,fname='correlations_',key='00')
    ssro_result_01=get_MBI_readout_result(dp,fname='correlations_',key='01')
    ssro_result_10=get_MBI_readout_result(dp,fname='correlations_',key='10')
    ssro_result_11=get_MBI_readout_result(dp,fname='correlations_',key='11')

    reps=ssro_result_01['y']+ssro_result_11['y']
    x=ssro_result_11['x']
    print ssro_result_11['y']
    y_norm=ssro_result_11['y']/float(reps)
    if Nuclcor == False:
        [ssro_ro_cor,ussro_ro_cor]=get_electron_ROC(y_norm,reps,ssro_calib_folder=get_latest_data('SSRO',datapath,date=currdate))
    else:
        [ssro_ro_cor,ussro_ro_cor]=get_nuclear_ROC(y_norm,reps,ssro_calib_folder=get_latest_data('SSRO',datapath,date=currdate))
    plot_ssro_vs_sweep(x,ssro_ro_cor,ussro_ro_cor,ssro_result['xname'],yname,currdate+'-'+fname+title,dp,save=save)

def get_feedback_data(foldername,filename='Spin_RO',d=''):
    datapath=get_latest_data(foldername,date=d)
    
    files = os.listdir(datapath)
    f=filename
    for k in files:
        if f in k:
            spin_ro_file = k

    data = np.load(datapath+'\\'+spin_ro_file)
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

    data_corr={}
    data_corr['sweep_par']=data['sweep_par']
    data_corr['sweep_par_name']=data['sweep_par_name']
    data_corr['FinalRO_SN'],data_corr['uFinalRO_SN']= get_nuclear_ROC(data_norm['FinalRO_SN'],data['SN'],get_latest_data('SSRO',date=d))
    data_corr['FinalRO_FS'],data_corr['uFinalRO_FS']= get_nuclear_ROC(data_norm['FinalRO_FS'],data['FS'],get_latest_data('SSRO',date=d))
    data_corr['FinalRO_FF'],data_corr['uFinalRO_FF']= get_nuclear_ROC(data_norm['FinalRO_FF'],data['FF'],get_latest_data('SSRO',date=d))
    data_corr['FinalRO_Succes'],data_corr['uFinalRO_Succes']= get_nuclear_ROC(data_norm['FinalRO_Succes'],
            data['FS']+data['SN'],get_latest_data('SSRO',date=d))
    data_corr['FinalRO_All'],data_corr['uFinalRO_All']= get_nuclear_ROC(data_norm['FinalRO_All'],reps,get_latest_data('SSRO',date=d))
    return data_norm,data_corr,datapath

def plot_feedback(foldername, filename='Spin_RO',d=''):
    data_norm, data_corr,dp = get_feedback_data (foldername, filename,d=d)
    x = data_norm['sweep_par']
    fig=plt.figure(2)
    fig.clf()
    plot_ssro_vs_sweep(x,data_norm['SN'],data_norm['uSN'],ylim=[0,1],xname=data_norm['sweep_par_name'],yname='P(click)',title='',label='1st msmnt',datapath=dp,nr=2,splot=221,save=False)
    plot_ssro_vs_sweep(x,data_norm['FS'],data_norm['uFS'],ylim=[0,1],xname=data_norm['sweep_par_name'],yname='P(click)',title='',label='2nd msmnt',datapath=dp,nr=2,splot=221,save=False)
    plot_ssro_vs_sweep(x,data_norm['SN']+(1-data_norm['SN'])*data_norm['FS'],(data_norm['uSN']+data_norm['uFS'])/2.,ylim=[0,0.5],xname=data_norm['sweep_par_name'],yname='P(click)',title='',label='P(succes) total',datapath=dp,nr=2,splot=221,save=False)
    plot_ssro_vs_sweep(x,data_corr['FinalRO_SN'],data_corr['uFinalRO_SN'],ylim=[0,1],xname=data_norm['sweep_par_name'],yname='P(mI=-1)',title='Final RO SN',datapath=dp,nr=2,splot=222,save=False)
    plot_ssro_vs_sweep(x,data_corr['FinalRO_FS'],data_corr['uFinalRO_FS'],ylim=[0,1],xname=data_norm['sweep_par_name'],yname='P(mI=-1)',title='Final RO FS',datapath=dp,nr=2,splot=223,save=False)
    plot_ssro_vs_sweep(x,data_corr['FinalRO_FF'],data_corr['uFinalRO_FF'],ylim=[0,1],xname=data_norm['sweep_par_name'],yname='P(mI=-1)',title='Final RO FF',datapath=dp,nr=2,splot=224,save=True)

    plot_ssro_vs_sweep(x,data_corr['FinalRO_All'],data_corr['uFinalRO_All'],ylim=[0,1],xname=data_norm['sweep_par_name'],yname='P(mI=-1)',title='Final RO Compare',label='All',datapath=dp,nr=3,splot=111,save=False)
    plot_ssro_vs_sweep(x,data_corr['FinalRO_Succes'],data_corr['uFinalRO_Succes'],ylim=[0,1],xname=data_norm['sweep_par_name'],yname='P(mI=-1)',title='Final RO Compare',label='Feedback',datapath=dp,nr=3,splot=111,save=False)
    plot_ssro_vs_sweep(x,data_corr['FinalRO_SN'],data_corr['uFinalRO_SN'],ylim=[0,1],xname=data_norm['sweep_par_name'],yname='P(mI=-1)',title='Final RO Compare',label='Heralded',datapath=dp,nr=3,splot=111,save=True)

    return data_norm,data_corr


##############################################
def fit_sin(x,y,uy,savepath='',fixed=[],fix_param=[],plot_fft=False,do_plot=True):
    FFT = fft.fft(y)
    N = int(len(x))
    timestep = (x[2]-x[1])
    freq = fft.fftfreq(N,d = timestep)

    #Remove offset:
    FFT[freq == 0] = 0
    if plot_fft:
        
        figure1 = plt.figure(1)
        plt.bar(freq*1E3,abs(FFT), width = 0.4*1E3*(freq[1]-freq[0]), align = 'center')
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Amplitude')
        plt.title('FFT (offset removed)')
        
        figure1.savefig(savepath+'\\fft_signal_rabi.png')

    freq_guess = freq[find_nearest(abs(FFT),abs(FFT).max())]
    amp_guess = (y.max()+y.min())/2.0
    offset_guess = y.min()+(y.max()+y.min())/2.0
    phase_guess=0
    if 0 in fixed:
        freq_guess=fix_param[0]
    figure3 = plt.figure(3)
    plt.clf()
    ax1 = figure3.add_subplot(111)

    
    fit_result = fit.fit1d(x, y, rabi.fit_rabi_simple, 
               freq_guess, amp_guess, offset_guess, phase_guess,fixed=fixed,
            do_print = do_plot , ret = True)
    
    if fit_result:
        fit_result['yerr']=uy  
        if do_plot:
            plot.plot_fit1d(fit_result,np.linspace(x.min(),x.max(),201),plot_data=do_plot)
    return fit_result



def plot_data_MBI(datapath,fid=(0.829,0.984),fiderr=(5.05e-03,1.82e-03), fit_data = True, title='',with_detuning = False, save = True):

    
    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    do_weak=False
    do_cond=False
    for k in files:
        if 'statics_and_parameters.npz' in k:
            stats_params_file = k
        if '0_Spin_RO.npz' in k:
            spin_ro_file = k
        if 'SP_histogram.npz' in k:
            sp_file = k
        if 'weak_Spin_RO.npz' in k:
            weak_spin_ro_file = k
            do_weak=True
        if 'cond_Spin_RO.npz' in k:
            cond_spin_ro_file = k
            do_cond=True

    e = np.load(datapath+'\\'+stats_params_file)
    par = e['sweep_par']


    noof_datapoints = len(par)
    if 'sweep_par_name' in e.keys():
        print e['sweep_par_name']
    else:
        print "Sweep parameter name not defined."
    
 
    noof_reps = e['completed_repetitions']
    if 'sweep_par_name' in e.keys():
        sp_name=e['sweep_par_name']
    else:
        sp_name='Total free evolution time [ns]'

    e.close()

    
    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    f = np.load(datapath+'\\'+spin_ro_file)
    
  

    t = f['time']

    idx = 0

    SSRO_readout = f['SSRO_counts']/float(noof_reps)
    SSRO_readout_corr = zeros(len(SSRO_readout))
    retdata = zeros(len(SSRO_readout))
    readout_error = zeros(len(SSRO_readout))
    print len(SSRO_readout)
    if len(par) == 1:
        par = np.linspace(1,75,75)
        sp_name='NR of RO steps'
    par_min = par.min()
    par_max = par.max()
    print len(par)
    for i in arange(len(SSRO_readout)):
        ms0_events = SSRO_readout[i]*noof_reps
        ms1_events = noof_reps*(1-SSRO_readout[i])
        corr = SSRO_correct(array([ms1_events,ms0_events]),F0=fid[0],F1=fid[1],F0err=fiderr[0],F1err=fiderr[1])
        SSRO_readout_corr[i]=corr[0][0]
        retdata[i]=corr[0][0]
        readout_error[i] = corr[1][0]  


    figure2 = plt.figure(2)
    plt.clf()
    plt.plot(par,SSRO_readout, 'sk')

    plt.ylim([0,1])
    plt.xlabel(sp_name)
    plt.ylabel('Fraction of events with > 0 counts')
    plt.title(datapath)
    plt.text(1.01*par_min,1.1*max(SSRO_readout),'')
    if save:
        figure2.savefig(datapath+'\\histogram_integrated_SSRO.png')

    
    figure3 = plt.figure(3)
    plt.clf()
    plt.errorbar(par,SSRO_readout_corr, fmt='o',yerr=readout_error)

    plt.ylim([0,1.1])
    plt.xlabel(sp_name)
    plt.ylabel('P(ms=0), corrected for readout error')
    plt.title(os.path.split(datapath)[1])
    plt.text(1.01*par_min,1.1*max(SSRO_readout_corr),'')
    if save:
        figure3.savefig(datapath+'\\SSRO_corrected.png')    

    x = 6.0
    y = 8.0

    if do_weak:
        
        f = np.load(datapath+'\\'+weak_spin_ro_file)
        
      

        t = f['time']

        idx = 0

        SSRO_weak_readout = f['SSRO_weak_counts']/float(noof_reps)
        SSRO_weak_readout_corr = zeros(len(SSRO_weak_readout))
        readout_error = zeros(len(SSRO_weak_readout))
        print len(SSRO_weak_readout)
        if len(par) == 1:
            par = np.linspace(1,75,75)
            sp_name='NR of RO steps'
        par_min = par.min()
        par_max = par.max()
        print len(par)
        for i in arange(len(SSRO_weak_readout)):
            ms0_events = SSRO_weak_readout[i]*noof_reps
            ms1_events = noof_reps*(1-SSRO_weak_readout[i])
            corr = SSRO_correct(array([ms1_events,ms0_events]),F0=fid[0],F1=fid[1],F0err=fiderr[0],F1err=fiderr[1])
            SSRO_weak_readout_corr[i]=corr[0][0]
            readout_error[i] = corr[1][0]  


        figure2 = plt.figure(2)
        plt.clf()
        plt.plot(par,SSRO_weak_readout, 'sk')

        plt.ylim([0,1])
        plt.xlabel(sp_name)
        plt.ylabel('Fraction of events with > 0 counts')
        plt.title(datapath)
        plt.text(1.01*par_min,1.1*max(SSRO_weak_readout),'')
        if save:
            figure2.savefig(datapath+'\\histogram_integrated_SSRO_weak.png')

        
        figure3 = plt.figure(3)
        plt.clf()
        plt.errorbar(par,SSRO_weak_readout_corr, fmt='o',yerr=readout_error)

        plt.ylim([0,1.1])
        plt.xlabel(sp_name)
        plt.ylabel('P(ms=0) weak meas, corrected for readout error')
        plt.title(datapath)
        plt.text(1.01*par_min,1.1*max(SSRO_weak_readout_corr),'')
        if save:
            figure3.savefig(datapath+'\\SSRO_corrected_weak.png')    

        x = 6.0
        y = 8.0
    if (do_cond) and (sum(SSRO_weak_readout)>.1):
        
        f = np.load(datapath+'\\'+cond_spin_ro_file)
        
      

        t = f['time']

        idx = 0

        SSRO_cond_readout = f['SSRO_cond_counts']
        SSRO_cond_readout_corr = zeros(len(SSRO_cond_readout))
        correlation = zeros(len(SSRO_cond_readout))
        readout_error = zeros(len(SSRO_cond_readout))
        print len(SSRO_cond_readout)
        if len(par) == 1:
            par = np.linspace(1,75,75)
            sp_name='NR of RO steps'
        par_min = par.min()
        par_max = par.max()
        print len(par)
        for i in arange(len(SSRO_cond_readout)):
            ms0_events = SSRO_cond_readout[i]
            
            ms1_events = float(noof_reps)*SSRO_weak_readout[i]*(1-SSRO_cond_readout[i]/SSRO_weak_readout[i]/float(noof_reps))
            correlation[i]= ((float(SSRO_cond_readout[i])/float(SSRO_weak_readout[i])/float(noof_reps))*2)-1
            corr = SSRO_correct(array([ms1_events,ms0_events]),F0=fid[0],F1=fid[1],F0err=fiderr[0],F1err=fiderr[1])
            SSRO_cond_readout_corr[i]=corr[0][0]
            readout_error[i] = corr[1][0]  
        print 'correlation'
        print correlation

        figure2 = plt.figure(2)
        plt.clf()
        plt.plot(par,correlation, 'sk')

        plt.ylim([-1,1])
        plt.xlabel(sp_name)
        plt.ylabel('Correlation with weak meas result')
        plt.title(datapath)
        plt.text(1.01*par_min,1.1*max(correlation),'')
        if save:
            figure2.savefig(datapath+'\\histogram_integrated_SSRO_cond.png')

        
        figure3 = plt.figure(3)
        plt.clf()
        plt.errorbar(par,SSRO_cond_readout_corr, fmt='o',yerr=readout_error)

        plt.ylim([0,1])
        plt.xlabel(sp_name)
        plt.ylabel('P(ms=0), corrected for readout error')
        plt.title(datapath)
        plt.text(1.01*par_min,1.1*max(SSRO_cond_readout_corr),'')
        if save:
            figure3.savefig(datapath+'\\SSRO_corrected_cond.png')    

        x = 6.0
        y = 8.0
    ###########################################
    ######## SPIN PUMPING #####################
    ###########################################
    v = np.load(datapath+'\\'+sp_file)
    sp_counts = v['counts']
    sp_time = v['time']

    offset_guess = sp_counts[len(sp_counts)-1]
    init_amp_guess = sp_counts[2]
    decay_guess = 10
    print 'reloaded'
    figure5 = plt.figure(10)
    fit.fit1d(sp_time/1E3, sp_counts, common.fit_exp_decay_with_offset, 
            offset_guess, init_amp_guess, decay_guess,
            do_plot = True, do_print = True, newfig = False,
            plot_fitparams_xy = (0.5,0.5),ret=False)
    
    plt.plot(sp_time/1E3,sp_counts,'sg')
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Integrated counts')
    plt.title('Spin pumping')
    v.close()
    print par
    if save:
        figure5.savefig(datapath+'\\spin_pumping.png')

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: MW length (ns)\tCol1: Integrated counts\tCol2: SSRO Readout corrected\tCol3:error SSRO Readout Cor\n'
        col_vals = str()
        for k in arange(len(par)):
            col_vals += num2str(par[k],2)+'\t'+str(SSRO_readout[k])+'\t' + str(SSRO_readout_corr[k]) +'\t'+str(readout_error[k]) +'\n'
        fo = open(datapath+'\\SSRO_readout.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        #curr_date = '#'+time.ctime()+'\n'
        #col_names = '#Col0: MW length (ns)\tCol1: SSRO corrected\n'
        #col_vals = str()
        #for k in arange(noof_datapoints):
        #    col_vals += num2str(par[k],2)+'\t'+str(SSRO_readout[k])+'\n'
        #fo = open(datapath+'\\SSRO_readout.dat', "w")
        #for item in [curr_date, col_names, col_vals]:
        #    fo.writelines(item)
        #fo.close()
    data={}    
    data['x']=par
    data['y']=retdata
    data['yerr']=readout_error

    if (do_cond) and (sum(SSRO_weak_readout)>.1):
        data['y_cond']=SSRO_cond_readout_corr

    if (do_cond) and (sum(SSRO_weak_readout)>.1):
        data['y_weak']=SSRO_weak_readout_corr
    return data

def plot_data_MBI(datapath,fid=(0.829,0.984),fiderr=(5.05e-03,1.82e-03), fit_data = True, title='',with_detuning = False, save = True):
    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    do_weak=False
    do_cond=False
    for k in files:
        if 'statics_and_parameters.npz' in k:
            stats_params_file = k
        if '0_Spin_RO.npz' in k:
            spin_ro_file = k
        if 'SP_histogram.npz' in k:
            sp_file = k
        if 'weak_Spin_RO.npz' in k:
            weak_spin_ro_file = k
            do_weak=True
        if 'cond_Spin_RO.npz' in k:
            cond_spin_ro_file = k
            do_cond=True

    e = np.load(datapath+'\\'+stats_params_file)
    par = e['sweep_par']
    par_min = par.min()
    par_max=par.max()
    
    
    print par_min*1e-6
    noof_datapoints = e['noof_datapoints']
    noof_reps = e['completed_repetitions']
    if 'sweep_par_name' in e.keys():
        sp_name=e['sweep_par_name']
    else:
        sp_name='Total free evolution time [ns]'

    

    
    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    f = np.load(datapath+'\\'+spin_ro_file)
    #raw_counts = f['counts']
    SSRO_counts = f['SSRO_counts']
    repetitions = f['sweep_axis']
    t = f['time']

    tot_size = len(repetitions)
    reps_per_point = tot_size/float(noof_datapoints)

    idx = 0
    counts_during_readout = zeros(noof_datapoints)
    par = linspace(par_min,par_max,noof_datapoints)
    print par
    #counts_during_readout = sum(raw_counts, axis = 1)
    #counts_during
    #SSRO_readout = sum(SSRO_counts, axis = 1)/float(noof_reps)
    SSRO_readout=SSRO_counts
    SSRO_readout_corr = zeros(len(SSRO_readout))
    readout_error = zeros(len(SSRO_readout))

    if 'sweep_par_name' in e.keys():
        print e['sweep_par_name']
    else:
        print "Sweep parameter name not defined."
    
 
    noof_reps = e['completed_repetitions']
    if 'sweep_par_name' in e.keys():
        sp_name=e['sweep_par_name']
    else:
        sp_name='Total free evolution time [ns]'

    e.close()

    
    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    f = np.load(datapath+'\\'+spin_ro_file)
    
  

    t = f['time']

    idx = 0

    SSRO_readout = f['SSRO_counts']/float(noof_reps)
    SSRO_readout_corr = zeros(len(SSRO_readout))
    retdata = zeros(len(SSRO_readout))
    readout_error = zeros(len(SSRO_readout))
    print len(SSRO_readout)
    if len(par) == 1:
        par = np.linspace(1,75,75)
        sp_name='NR of RO steps'
    par_min = par.min()
    par_max = par.max()
    print len(par)
    for i in arange(len(SSRO_readout)):
        ms0_events = SSRO_readout[i]*noof_reps
        ms1_events = noof_reps*(1-SSRO_readout[i])
        corr = SSRO_correct(array([ms1_events,ms0_events]),F0=fid[0],F1=fid[1],F0err=fiderr[0],F1err=fiderr[1])
        SSRO_readout_corr[i]=corr[0][0]
        retdata[i]=corr[0][0]
        readout_error[i] = corr[1][0]  


    figure2 = plt.figure(2)
    plt.clf()
    plt.plot(par,SSRO_readout, 'sk')

    plt.ylim([0,1])
    plt.xlabel(sp_name)
    plt.ylabel('Fraction of events with > 0 counts')
    plt.title(datapath)
    plt.text(1.01*par_min,1.1*max(SSRO_readout),'')
    if save:
        figure2.savefig(datapath+'\\histogram_integrated_SSRO.png')

    
    figure3 = plt.figure(3)
    plt.clf()
    plt.errorbar(par,SSRO_readout_corr, fmt='o',yerr=readout_error)

    plt.ylim([0,1.1])
    plt.xlabel(sp_name)
    plt.ylabel('P(ms=0), corrected for readout error')
    plt.title(datapath)
    plt.text(1.01*par_min,1.1*max(SSRO_readout_corr),'')
    if save:
        figure3.savefig(datapath+'\\SSRO_corrected.png')    

    x = 6.0
    y = 8.0

    if do_weak:
        
        f = np.load(datapath+'\\'+weak_spin_ro_file)
        
      

        t = f['time']

        idx = 0

        SSRO_weak_readout = f['SSRO_weak_counts']/float(noof_reps)
        SSRO_weak_readout_corr = zeros(len(SSRO_weak_readout))
        readout_error = zeros(len(SSRO_weak_readout))
        print len(SSRO_weak_readout)
        if len(par) == 1:
            par = np.linspace(1,75,75)
            sp_name='NR of RO steps'
        par_min = par.min()
        par_max = par.max()
        print len(par)
        for i in arange(len(SSRO_weak_readout)):
            ms0_events = SSRO_weak_readout[i]*noof_reps
            ms1_events = noof_reps*(1-SSRO_weak_readout[i])
            corr = SSRO_correct(array([ms1_events,ms0_events]),F0=fid[0],F1=fid[1],F0err=fiderr[0],F1err=fiderr[1])
            SSRO_weak_readout_corr[i]=corr[0][0]
            readout_error[i] = corr[1][0]  


        figure2 = plt.figure(2)
        plt.clf()
        plt.plot(par,SSRO_weak_readout, 'sk')

        plt.ylim([0,1])
        plt.xlabel(sp_name)
        plt.ylabel('Fraction of events with > 0 counts')
        plt.title(datapath)
        plt.text(1.01*par_min,1.1*max(SSRO_weak_readout),'')
        if save:
            figure2.savefig(datapath+'\\histogram_integrated_SSRO_weak.png')

        
        figure3 = plt.figure(3)
        plt.clf()
        plt.errorbar(par,SSRO_weak_readout_corr, fmt='o',yerr=readout_error)

        plt.ylim([0,1.1])
        plt.xlabel(sp_name)
        plt.ylabel('P(ms=0) weak meas, corrected for readout error')
        plt.title(datapath)
        plt.text(1.01*par_min,1.1*max(SSRO_weak_readout_corr),'')
        if save:
            figure3.savefig(datapath+'\\SSRO_corrected_weak.png')    

        x = 6.0
        y = 8.0
    if (do_cond) and (sum(SSRO_weak_readout)>.1):
        
        f = np.load(datapath+'\\'+cond_spin_ro_file)
        
      

        t = f['time']

        idx = 0

        SSRO_cond_readout = f['SSRO_cond_counts']
        SSRO_cond_readout_corr = zeros(len(SSRO_cond_readout))
        correlation = zeros(len(SSRO_cond_readout))
        readout_error = zeros(len(SSRO_cond_readout))
        print len(SSRO_cond_readout)
        if len(par) == 1:
            par = np.linspace(1,75,75)
            sp_name='NR of RO steps'
        par_min = par.min()
        par_max = par.max()
        print len(par)
        for i in arange(len(SSRO_cond_readout)):
            ms0_events = SSRO_cond_readout[i]
            
            ms1_events = float(noof_reps)*SSRO_weak_readout[i]*(1-SSRO_cond_readout[i]/SSRO_weak_readout[i]/float(noof_reps))
            correlation[i]= ((float(SSRO_cond_readout[i])/float(SSRO_weak_readout[i])/float(noof_reps))*2)-1
            corr = SSRO_correct(array([ms1_events,ms0_events]),F0=fid[0],F1=fid[1],F0err=fiderr[0],F1err=fiderr[1])
            SSRO_cond_readout_corr[i]=corr[0][0]
            readout_error[i] = corr[1][0]  
        print 'correlation'
        print correlation

        figure2 = plt.figure(2)
        plt.clf()
        plt.plot(par,correlation, 'sk')

        plt.ylim([-1,1])
        plt.xlabel(sp_name)
        plt.ylabel('Correlation with weak meas result')
        plt.title(datapath)
        plt.text(1.01*par_min,1.1*max(correlation),'')
        if save:
            figure2.savefig(datapath+'\\histogram_integrated_SSRO_cond.png')

        
        figure3 = plt.figure(3)
        plt.clf()
        plt.errorbar(par,SSRO_cond_readout_corr, fmt='o',yerr=readout_error)

        plt.ylim([0,1])
        plt.xlabel(sp_name)
        plt.ylabel('P(ms=0), corrected for readout error')
        plt.title(datapath)
        plt.text(1.01*par_min,1.1*max(SSRO_cond_readout_corr),'')
        if save:
            figure3.savefig(datapath+'\\SSRO_corrected_cond.png')    

        x = 6.0
        y = 8.0
    ###########################################
    ######## SPIN PUMPING #####################
    ###########################################
    v = np.load(datapath+'\\'+sp_file)
    sp_counts = v['counts']
    sp_time = v['time']

    offset_guess = sp_counts[len(sp_counts)-1]
    init_amp_guess = sp_counts[2]
    decay_guess = 10
    print 'reloaded'
    figure5 = plt.figure(10)
    fit.fit1d(sp_time/1E3, sp_counts, common.fit_exp_decay_with_offset, 
            offset_guess, init_amp_guess, decay_guess, do_print = True,ret=False,
            plot_fitparams_xy = (0.5,0.5))
    
    plt.plot(sp_time/1E3,sp_counts,'sg')
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Integrated counts')
    plt.title('Spin pumping')
    v.close()
    print par
    if save:
        figure5.savefig(datapath+'\\spin_pumping.png')

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: MW length (ns)\tCol1: Integrated counts\tCol2: SSRO Readout corrected\tCol3:error SSRO Readout Cor\n'
        col_vals = str()
        for k in arange(len(par)):
            col_vals += num2str(par[k],2)+'\t'+str(SSRO_readout[k])+'\t' + str(SSRO_readout_corr[k]) +'\t'+str(readout_error[k]) +'\n'
        fo = open(datapath+'\\SSRO_readout.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        #curr_date = '#'+time.ctime()+'\n'
        #col_names = '#Col0: MW length (ns)\tCol1: SSRO corrected\n'
        #col_vals = str()
        #for k in arange(noof_datapoints):
        #    col_vals += num2str(par[k],2)+'\t'+str(SSRO_readout[k])+'\n'
        #fo = open(datapath+'\\SSRO_readout.dat', "w")
        #for item in [curr_date, col_names, col_vals]:
        #    fo.writelines(item)
        #fo.close()
    data={}    
    data['x']=par
    data['y']=retdata
    data['yerr']=readout_error

    if (do_cond) and (sum(SSRO_weak_readout)>.1):
        data['y_cond']=SSRO_cond_readout_corr

    if (do_cond) and (sum(SSRO_weak_readout)>.1):
        data['y_weak']=SSRO_weak_readout_corr
    return data



def plot_rabi(datapath, fit_data = True, with_detuning = False, save = True, ro_correct=False,fid=(0.8554,0.984),fiderr=(4.97e-03,1.73e-03)):

    plt.close('all')
    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    
    for k in files:
        if 'statics_and_parameters.npz' in k:
            stats_params_file = k
        if 'Spin_RO.npz' in k:
            spin_ro_file = k
        if 'SP_histogram.npz' in k:
            sp_file = k

    e = np.load(datapath+'\\'+stats_params_file)
    f_drive = e['mw_drive_freq']
    mwpower = e['mw_power']
    mw_min_len = e['sweep_par'].min()
    mw_max_len = e['sweep_par'].max()
    noof_datapoints = e['noof_datapoints']
    noof_reps = e['completed_repetitions']
    mw_len = e['sweep_par'] 
    e.close()

    
    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    f = np.load(datapath+'\\'+spin_ro_file)
    #raw_counts = f['counts']
    SSRO_counts = f['SSRO_counts']
    #repetitions = f['sweep_axis']
    t = f['time']


    idx = 0
    counts_during_readout = zeros(noof_datapoints)
    #mw_len = linspace(mw_min_len,mw_max_len,noof_datapoints)
    #counts_during_readout = sum(raw_counts, axis = 1)
    #SSRO_readout = sum(SSRO_counts, axis = 1)/float(noof_reps)
    counts_during_readout=SSRO_counts/float(noof_reps)
    SSRO_readout=SSRO_counts/float(noof_reps)

    SSRO_readout_corr = zeros(len(SSRO_readout))
    readout_error = zeros(len(SSRO_readout))
    for i in arange(len(SSRO_readout)):
        ms0_events = SSRO_readout[i]*noof_reps
        ms1_events = noof_reps*(1-SSRO_readout[i])
        corr = SSRO_correct(array([ms1_events,ms0_events]),F0=fid[0],F1=fid[1],F0err=fiderr[0],F1err=fiderr[1])
        SSRO_readout_corr[i]=corr[0][0]
        readout_error[i] = corr[1][0] 

    #########################################
    ############ FITTING ####################
    #########################################
    
    if fit_data:
        FFT = fft.fft(SSRO_readout_corr)
        N = int(noof_datapoints)
        timestep = (mw_max_len-mw_min_len)/float(noof_datapoints-1)
        freq = fft.fftfreq(N,d = timestep)

        #Remove offset:
        FFT[freq == 0] = 0

        figure1 = plt.figure(1)
        plt.bar(freq*1E3,abs(FFT), width = 0.4*1E3*(freq[1]-freq[0]), align = 'center')
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Amplitude')
        plt.title('FFT (offset removed)')
        
        if save:
            figure1.savefig(datapath+'\\fft_signal_rabi.png')

        freq_guess = freq[find_nearest(abs(FFT),abs(FFT).max())]
        amp_guess = (SSRO_readout_corr.max()+SSRO_readout_corr.min())/2.0
        offset_guess = SSRO_readout_corr.min()+(SSRO_readout_corr.max()+\
                SSRO_readout_corr.min())/2.0
        phase_guess = 0

    figure3 = plt.figure(3)
    plt.clf()
    ax1 = figure3.add_subplot(111)
   
    if fit_data and with_detuning:
        tau_guess = 200
        fit.fit1d(mw_len, counts_during_readout, 
                rabi.fit_rabi_multiple_detunings, 
                amp_guess, offset_guess, freq_guess, tau_guess,
                (0,0),(2.2E-3,0),(2*2.2E-3,0),
                do_print = True,ret=False)
    elif fit_data:
        fit_result = fit.fit1d(mw_len, SSRO_readout_corr, rabi.fit_rabi_simple, 
                freq_guess, amp_guess, offset_guess, phase_guess,
                do_print = True , ret = True)
        fit_result['yerr']=readout_error
        plot.plot_fit1d(fit_result,np.linspace(mw_len.min(),mw_len.max(),201))
        #pi_pulse_len = 0.5*1/fit_result['params_dict']['f']
        #pi_2_pulse_len = 0.5*pi_pulse_len

        #print 'pi pulse length = ', pi_pulse_len
        #print 'pi/2 pulse length = ',pi_2_pulse_len
    

    #ax2 = figure2.add_subplot(111, sharex=ax1, frameon=False)
    #ax2.plot(mw_len,counts_during_readout, 'sk')

    plt.xlabel('MW length (ns)')
    plt.ylabel('Integrated counts')
    plt.title('MW length sweep, driving $f$ ='+num2str(f_drive/1E6,1)+\
            ' MHz, power = '+num2str(mwpower,0)+' dBm')
    plt.text(0.1*(mw_max_len+mw_min_len),max(counts_during_readout),datapath)
    if save:
        figure3.savefig(datapath+'\\histogram_integrated.png')


    figure2 = plt.figure(2)
    plt.clf()
    ax1 = figure2.add_subplot(111)
    #plt.plot(mw_len,SSRO_readout, 'sk')


    #ax2 = figure2.add_subplot(111, sharex=ax1, frameon=False)
    #ax2.plot(mw_len,counts_during_readout, 'sk')
    plt.ylim([0,1])
    plt.xlabel('MW length (ns)')
    plt.ylabel('Fraction of events N >0')
    plt.title('MW length sweep, driving $f$ ='+num2str(f_drive/1E6,1)+\
            ' MHz, power = '+num2str(mwpower,0)+' dBm')
    plt.text(0.1*(mw_max_len+mw_min_len),max(counts_during_readout),datapath)
    if save:
        figure2.savefig(datapath+'\\histogram_integrated_SSRO.png')

    x = 6.0
    y = 8.0


    ###########################################
    ######## RO correct   #####################
    ###########################################

    if ro_correct:

       
        
        figure5=plt.figure()

        amp_guess_corr = (SSRO_readout_corr.max()+SSRO_readout_corr.min())/2.0
        offset_guess_corr =SSRO_readout_corr.min()+(SSRO_readout_corr.max()+\
                SSRO_readout_corr.min())/2.0
        phase_guess = 0
        if fit_data:
            [fit_result, p] = fit.fit1d(mw_len, SSRO_readout_corr, rabi.fit_rabi_simple, 
                    freq_guess, amp_guess_corr, offset_guess_corr, phase_guess,
                    do_plot = True, do_print = True, newfig = False, ret = True)

        else: plt.plot(mw_len,SSRO_readout_corr,'bo')
        plt.title('MW length sweep, driving $f$ ='+num2str(f_drive/1E6,1)+\
                ' MHz, power = '+num2str(mwpower,0)+' dBm')
        plt.text(0.1*(mw_max_len+mw_min_len),0.9,datapath)
        plt.xlabel('MW length (ns)')
        plt.ylabel('Read out corrected P_ms=0')
        plt.ylim(ymin=0)

        if save:
            figure5.savefig(datapath+'\\Rabi_SSRO_corrected.png')


    ###########################################
    ######## SPIN PUMPING #####################
    ###########################################
    v = np.load(datapath+'\\'+sp_file)
    sp_counts = v['counts']
    sp_time = v['time']

    offset_guess = sp_counts[len(sp_counts)-1]
    init_amp_guess = sp_counts[2]
    decay_guess = 10

    figure6 = plt.figure()
    fit.fit1d(sp_time/1E3, sp_counts, common.fit_exp_decay_with_offset, 
            offset_guess, init_amp_guess, decay_guess,
            do_plot = True, do_print = True, newfig = False,
            plot_fitparams_xy = (0.5,0.5),ret=False)
    
    plt.plot(sp_time/1E3,sp_counts,'sg')
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Integrated counts')
    plt.title('Spin pumping')
    v.close()
    if save:
        figure6.savefig(datapath+'\\spin_pumping.png')

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: MW length (ns)\tCol1: Integrated counts\n'
        col_vals = str()
        for k in arange(noof_datapoints):
            col_vals += num2str(mw_len[k],2)+'\t'+num2str(counts_during_readout[k],0)+'\n'
        fo = open(datapath+'\\integrated_histogram.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()

    return  [fit_result]


def plot_dark_esr(datapath, fit_data = True, save = True, f_dip = 2.8295E9,d=''):
  
    if d=='':
       date=time.strftime('%Y%m%d')
    else:
       date=d
    plt.close('all')


    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    
    for k in files:
        if 'statics_and_parameters.npz' in k:
            stats_params_file = k
        if 'Spin_RO.npz' in k:
            spin_ro_file = k
        if 'SP_histogram.npz' in k:
            sp_file = k

    e = np.load(datapath+'\\'+stats_params_file)
    f_center = e['mw_center_freq']
    mwpower = e['mw_power']
    #mw_min_freq = e['min_mw_freq']
    #mw_max_freq = e['max_mw_freq']
    noof_datapoints = e['noof_datapoints']
    noof_reps = e['completed_repetitions']
    repetitions = e['completed_repetitions']
    if 'sweep_par' in e.keys():
        mw_freq=e['sweep_par']
    else:
        mw_min_freq = e['min_mw_freq']
        mw_max_freq = e['max_mw_freq']
        mw_freq = np.linspace(mw_min_freq,mw_max_freq,noof_datapoints)*1e-9
        
    e.close()

    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    f = np.load(datapath+'\\'+spin_ro_file)
    #raw_counts = f['counts']
    SSRO_counts = f['SSRO_counts']
    
    t = f['time']
    if 'sweep_par' in e.keys():
        SSRO_readout = SSRO_counts/float(noof_reps)
        counts_during_readout = SSRO_counts
    else:
        SSRO_readout = sum(SSRO_counts, axis = 1)/float(noof_reps)

    idx = 0




    [ssro_ro_cor,ussro_ro_cor]=get_electron_ROC(SSRO_readout,noof_reps,ssro_calib_folder=get_latest_data('SSRO',date=date))
    #########################################
    ############ FITTING ####################
    #########################################

    guess_offset = ssro_ro_cor.max()
    guess_A_min1 = 0.04
    guess_A_plus1 = 0.04
    guess_A_0 = 0.00
    guess_x0 = 2.8289
    guess_sigma = 0.00045

    splitting = 2.19e-3
    
    A_min1 = fit.Parameter(guess_A_min1, 'A_min1')
    A_plus1 = fit.Parameter(guess_A_plus1, 'A_plus1')
    A_0 = fit.Parameter(guess_A_0, 'A_0')
    o = fit.Parameter(guess_offset, 'o')
    x0 = fit.Parameter(guess_x0, 'x0')
    sigma = fit.Parameter(guess_sigma, 'sigma')
    def fitfunc(x):
        return o() - A_min1()*np.exp(-((x-(x0()-splitting))/sigma())**2) \
            - A_plus1()*np.exp(-((x-(x0()+splitting))/sigma())**2) \
            - A_0()*np.exp(-((x-x0())/sigma())**2) 
    figure2 = plt.figure(2)
    dp=datapath
    print dp

    ax=plot_ssro_vs_sweep(mw_freq,ssro_ro_cor,ussro_ro_cor,[0.4,1.05],'MW freq (GHz)','P(ms=0)',title='Dark_ESR',datapath=dp,nr=1)

    if fit_data:
        fit_result = fit.fit1d(mw_freq, ssro_ro_cor, None, p0 = [A_min1, A_plus1, A_0, o, x0, sigma],fitfunc = fitfunc, do_print=True, ret=True, fixed=[])
        plot.plot_fit1d(fit_result,np.linspace(mw_freq.min(),mw_freq.max(),751),
                ax=ax, plot_data=False,info_xy=(0.128,0.328))
        center_peak = fit_result['params_dict']['x0']




        print '-1: f = ', (center_peak - splitting), ' GHz'
        print '0: f = ', center_peak, ' GHz'
        print '1: f = ', (center_peak + splitting), ' GHz'
        print 'Populations: ',

    pol = error.Formula()
    a0, am1, ap1 = sympy.symbols('a0, am1, ap1')
    pol.formula = max(0.,am1) / (max(0.,a0) + max(0.,ap1) + max(0.,am1))
    pol.values[a0] = A_0()
    pol.values[am1] = A_min1()
    pol.values[ap1] = A_plus1()
    pol.uncertainties[a0] = fit_result['error_dict']['A_0']
    pol.uncertainties[am1] = fit_result['error_dict']['A_min1']
    pol.uncertainties[ap1] = fit_result['error_dict']['A_plus1']

    print 'Spin polarization mI=-1 = %.3f +/- %.3f' \
            % (float(pol.value()), float(pol.uncertainty()))

    plt.xlabel('MW frequency (GHz)')
    plt.ylabel('Fraction of events N >0')
    plt.title('MW frequency sweep, MW parking spot $f$ ='+num2str(f_center/1E6,1)+\
            ' MHz, power = '+num2str(mwpower,0)+' dBm')
   
    if save:
        figure2.savefig(datapath+'\\histogram_integrated.png')

    x = 6.0
    y = 8.0

   

    ###########################################
    ######## SPIN PUMPING #####################
    ###########################################

    v = np.load(datapath+'\\'+sp_file)
    sp_counts = v['counts']
    sp_time = v['time']

    offset_guess = sp_counts[len(sp_counts)-1]
    init_amp_guess = sp_counts[2]
    decay_guess = 10

    figure4 = plt.figure(4)
    fit.fit1d(sp_time/1E3, sp_counts, common.fit_exp_decay_with_offset, 
            offset_guess, init_amp_guess, decay_guess, do_print = True, ret=False,
            plot_fitparams_xy = (0.5,0.5))
    
    plt.plot(sp_time/1E3,sp_counts,'sg')
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Integrated counts')
    plt.title('Spin pumping')
    v.close()
    if save:
        figure4.savefig(datapath+'\\spin_pumping.png')

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: MW freq (GHz)\tCol1: Integrated counts\n'
        col_vals = str()
        for k in arange(noof_datapoints):
            col_vals += num2str(mw_freq[k]/1E9,10)+'\t'+num2str(SSRO_readout[k],0)+'\n'
        fo = open(datapath+'\\integrated_histogram.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()

        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: MW length (ns)\tCol1: Integrated counts\n'
        col_vals = str()
        for k in arange(noof_datapoints):
            col_vals += num2str(mw_freq[k]/1E9,10)+'\t'+str(SSRO_readout[k]) +'\n'
        fo = open(datapath+'\\SSRO_readout.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()

    return fit_result,ssro_ro_cor,ussro_ro_cor

def plot_ramsey(datapath, fid=(0.7611,0.9895),fiderr=(4.264e-03,1.019e-03),fit_data = True, save = True):

    plt.close('all')
    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    
    for k in files:
        if 'statics_and_parameters.npz' in k:
            stats_params_file = k
        if 'Spin_RO.npz' in k:
            spin_ro_file = k
        if 'SP_histogram.npz' in k:
            sp_file = k

    e = np.load(datapath+'\\'+stats_params_file)
    f_drive = e['mw_drive_freq']
    detuning = e['detuning']
    mwpower = e['mw_power']
    tau_min = e['min_wait_time']
    tau_max = e['max_wait_time']
    #pi_2_duration = e['pi_2_duration']
    pi_2_duration=25
    noof_datapoints = e['noof_datapoints']
    noof_reps = e['completed_repetitions']
    e.close()

    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    f = np.load(datapath+'\\'+spin_ro_file)
    raw_counts = f['counts']
    SSRO_counts = f['SSRO_counts']
    repetitions = f['sweep_axis']
    t = f['time']

    tot_size = len(repetitions)
    reps_per_point = tot_size/float(noof_datapoints)

    idx = 0
    counts_during_readout = zeros(noof_datapoints)
    tau = linspace(tau_min,tau_max,noof_datapoints)
    counts_during_readout = sum(raw_counts, axis = 1)
    SSRO_readout = sum(SSRO_counts, axis = 1)/float(noof_reps)
    SSRO_readout_corr = zeros(len(SSRO_readout))
    readout_error = zeros(len(SSRO_readout))
    SSRO_uncor=SSRO_readout

    for i in arange(len(SSRO_readout)):
        ms0_events = SSRO_readout[i]*noof_reps
        ms1_events = noof_reps*(1-SSRO_readout[i])
        corr = SSRO_correct(array([ms1_events,ms0_events]),F0=fid[0],F1=fid[1],F0err=fiderr[0],F1err=fiderr[1])
        SSRO_readout_corr[i]=corr[0][0]
        readout_error[i] = corr[1][0] 
    SSRO_readout=SSRO_readout_corr
    #########################################
    ############ FITTING ####################
    #########################################
    print fit_data
    if fit_data:
        FFT = fft.fft(SSRO_readout)
        N = int(noof_datapoints)
        timestep = (tau_max-tau_min)/float(noof_datapoints-1)
        freq = fft.fftfreq(N,d = timestep)

        #Remove offset:
        FFT[freq == 0] = 0

        figure1 = plt.figure(1)
        plt.bar(freq*1E3,abs(FFT), width = 0.4*1E3*(freq[1]-freq[0]), align = 'center')
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Amplitude')
        plt.title('FFT (offset removed)')
        
        if save:
            figure1.savefig(datapath+'\\fft_signal_rabi.png')

        offset_guess = SSRO_readout.min()+(SSRO_readout.max()+\
                SSRO_readout.min())/2.0
        tau_guess = 1000
        
        #modulations:
        mod_freq_guess = freq[find_nearest(abs(FFT),abs(FFT).max())]
        mod_amp_guess = (SSRO_readout.max()+SSRO_readout.min())/2.0
        mod_phase_guess = 0

        print 'modulation frequency guess = ',mod_freq_guess 
        print 'modulation amplitude guess = ',mod_amp_guess

    figure2 = plt.figure(2)
        
    if fit_data:
        tau_guess = 200
        fit.fit1d(tau, SSRO_readout, 
                ramsey.fit_ramsey_gaussian_decay, 
                tau_guess, offset_guess, 
                (mod_freq_guess,mod_amp_guess),
                do_plot = True, do_print = True, newfig = False)
    print SSRO_readout
    plt.plot(tau,SSRO_readout, 'sk')
    plt.xlabel('\Delta t$ (ns)')
    plt.ylabel('P(ms=0)')
    plt.title('Ramsey, driving $f$ ='+num2str(f_drive/1E6,1)+\
            ' MHz, power = '+num2str(mwpower,0)+' dBm,\n Detuning = '+\
            num2str(detuning/1E6,0)+' MHz, $\pi/2$ length = '+\
            num2str(pi_2_duration,0)+' ns')
    plt.text(0.1*(tau_max+tau_min),max(counts_during_readout),datapath)
    if save:
        figure2.savefig(datapath+'\\histogram_integrated.png')

    x = 6.0
    y = 8.0

    figure3 = plt.figure(figsize=(x,y))
    plt.pcolor(raw_counts, cmap = 'hot', antialiased=False)
    plt.xlabel('Readout time (us)')
    plt.ylabel('MW repetition number')
    plt.title('Total histogram, integrated over repetitions')
    plt.colorbar()
    if save:
        figure3.savefig(datapath+'\\histogram_counts_2d.png')

    f.close()
    

    ###########################################
    ######## SPIN PUMPING #####################
    ###########################################
    v = np.load(datapath+'\\'+sp_file)
    sp_counts = v['counts']
    sp_time = v['time']

    offset_guess = sp_counts[len(sp_counts)-1]
    init_amp_guess = sp_counts[2]
    decay_guess = 10

    figure4 = plt.figure(4)
    fit.fit1d(sp_time/1E3, sp_counts, common.fit_exp_decay_with_offset, 
            offset_guess, init_amp_guess, decay_guess,
            do_plot = True, do_print = True, newfig = False,
            plot_fitparams_xy = (0.5,0.5))
    
    plt.plot(sp_time/1E3,sp_counts,'sg')
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Integrated counts')
    plt.title('Spin pumping')
    v.close()
    if save:
        figure4.savefig(datapath+'\\spin_pumping.self.nr_of_shelving_pulsespng')

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: tau (ns)\tCol1: Integrated counts\tCol2: SSRO Readout corrected\tCol3:error SSRO Readout Cor\n'
        col_vals = str()
        for k in arange(noof_datapoints):
            col_vals += num2str(tau[k],2)+'\t'+str(SSRO_readout[k])+'\t' + str(SSRO_readout_corr[k]) +'\t'+str(readout_error[k]) +'\n'
        fo = open(datapath+'\\SSRO_readout.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()

    return True

def plot_SE(datapath,fid=(0.788,0.991),fiderr=(2.885e-03,6.659e-04),fit_data = True, exp_guess=2.8, nr_of_pulses=1,save = True):

    plt.close('all')
    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    
    for k in files:
        if 'statics_and_parameters.npz' in k:
            stats_params_file = k
        if 'Spin_RO.npz' in k:
            spin_ro_file = k
        if 'SP_histogram.npz' in k:
            sp_file = k

    e = np.load(datapath+'\\'+stats_params_file)
    f_drive = e['mw_drive_freq']
    mwpower = e['mw_power']
    noof_datapoints = e['noof_datapoints']
    if 'fe_min' in e.keys():
        tau_max = e['fe_max']/2.
        tau_min = e['fe_min']/2.
        if 'free_evol' in e.keys():
            
            tau_free_evol = e['free_evol']*1e-3
        else:
            tau_free_evol = np.linspace(tau_min,tau_max,noof_datapoints)*2*1e-6
    else:    
        if 'max_sweep_par' in e.keys():
            tau_max = e['max_sweep_par']
            tau_min = e['min_sweep_par']
        else:
            tau_max = e['tau_max']
            tau_min = e['tau_min']

        tau_free_evol = nr_of_pulses*linspace(tau_min,tau_max,noof_datapoints)*2*1e-6
    #tau_free_evol=np.linspace(0,4260,11)
    #tau_free_evol[0]=4
    noof_reps = e['completed_repetitions']
    e.close()

    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    f = np.load(os.path.join(datapath,spin_ro_file))
    raw_counts = f['counts']
    SSRO_counts = f['SSRO_counts']
    repetitions = f['sweep_axis']
    t = f['time']

    tot_size = len(repetitions)
    reps_per_point = tot_size/float(noof_datapoints)

    idx = 0
    counts_during_readout = zeros(noof_datapoints)
    counts_during_readout = sum(raw_counts, axis = 1)
    SSRO_readout = sum(SSRO_counts, axis = 1)/float(noof_reps)
    SSRO_readout_corr = zeros(len(SSRO_readout))
    readout_error = zeros(len(SSRO_readout))

    for i in arange(len(SSRO_readout)):
        ms0_events = SSRO_readout[i]*noof_reps
        ms1_events = noof_reps*(1-SSRO_readout[i])
        corr = SSRO_correct(array([ms1_events,ms0_events]),F0=fid[0],F1=fid[1],F0err=fiderr[0],F1err=fiderr[1])
        SSRO_readout_corr[i]=corr[0][0]
        readout_error[i] = corr[1][0] 

    #########################################
    ############ FITTING ####################
    #########################################
    #FIXME: write nice fitting routine for SE
    
    figure2 = plt.figure(2)
        

    plt.plot(tau_free_evol,counts_during_readout, 'sk')
    plt.xlabel('Free evolution time [us]')
    plt.ylabel('Integrated counts')
    plt.title('Spin Echo, driving $f$ ='+num2str(f_drive/1E6,1)+\
            ' MHz')
    plt.text(0.1*(max(tau_free_evol)+min(tau_free_evol)),max(counts_during_readout),datapath)
    if save:
        figure2.savefig(datapath+'\\histogram_integrated.png')

    figure4 = plt.figure(4)
  
    plt.plot(tau_free_evol,SSRO_readout, 'sk')
    plt.xlabel('Free evolution time [us]')
    plt.ylabel('Fraction of events >0 counts')
    plt.title('Spin Echo SSRO cor, driving $f$ ='+num2str(f_drive/1E6,1)+\
            ' MHz')
    plt.text(0.1*(max(tau_free_evol)+min(tau_free_evol)),max(SSRO_readout),datapath)
    if save:
        figure4.savefig(datapath+'\\SSRO_readout.png')    

    figure6 = plt.figure(6)
        
    if fit_data:
        tau_guess = 500
        offset_guess = SSRO_readout_corr.min()
        amp_guess=SSRO_readout_corr.max()-offset_guess

        print tau_guess
        print offset_guess
        print amp_guess
        print tau_guess
        print exp_guess
        print tau_free_evol

        fitresult=fit.fit1d(tau_free_evol, SSRO_readout_corr, 
                SE.fit_echo, tau_guess, amp_guess,offset_guess, 
                exp_guess, do_plot = True, do_print = True, newfig = False)
        print fitresult
    plt.errorbar(tau_free_evol,SSRO_readout_corr,fmt='sk',yerr=readout_error)
    plt.xlabel('Free evolution time [us]')
    plt.ylabel('P ms=0')
    plt.title('Spin Echo SSRO cor, driving $f$ ='+num2str(f_drive/1E6,1)+\
            ' MHz')
    plt.text(0.1*(max(tau_free_evol)+min(tau_free_evol)),max(SSRO_readout_corr),datapath)
    if save:
        figure6.savefig(datapath+'\\SSRO_readout_corr.png')  


    x = 6.0
    y = 8.0

    figure3 = plt.figure(figsize=(x,y))
    plt.pcolor(raw_counts, cmap = 'hot', antialiased=False)
    plt.xlabel('Readout time (us)')
    plt.ylabel('MW repetition number')
    plt.title('Total histogram, integrated over repetitions')
    plt.colorbar()
    if save:
        figure3.savefig(datapath+'\\histogram_counts_2d.png')

    f.close()
    

    ###########################################
    ######## SPIN PUMPING #####################
    ###########################################
    v = np.load(datapath+'\\'+sp_file)
    sp_counts = v['counts']
    sp_time = v['time']

    offset_guess = sp_counts[len(sp_counts)-1]
    init_amp_guess = sp_counts[2]
    decay_guess = 10

    figure5 = plt.figure(5)
    fit.fit1d(sp_time/1E3, sp_counts, common.fit_exp_decay_with_offset, 
            offset_guess, init_amp_guess, decay_guess,
            do_plot = True, do_print = True, newfig = False,
            plot_fitparams_xy = (0.5,0.5))
    
    plt.plot(sp_time/1E3,sp_counts,'sg')
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Integrated counts')
    plt.title('Spin pumping')
    v.close()
    if save:

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: MW length (ns)\tCol1: Integrated counts\tCol2: SSRO Readout corrected\tCol3: error SSRO Readout corrected\n'
        col_vals = str()
        for k in arange(noof_datapoints):
            col_vals += num2str(tau_free_evol[k],2)+'\t'+str(SSRO_readout[k])+'\t' + str(SSRO_readout_corr[k]) +'\t'+str(readout_error[k]) +'\n'
        fo = open(datapath+'\\SSRO_readout.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()



        figure5.savefig(datapath+'\\spin_pumping.png')

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: Interpulse delay (ns)\tCol1: Integrated counts\tCol2: SSRO corr\n'
        col_vals = str()
        for k in arange(noof_datapoints):
            col_vals += num2str(tau_free_evol[k],2)+'\t'+num2str(counts_during_readout[k],0)+'\t'+num2str(SSRO_readout_corr[k],2)+'\n'
        fo = open(datapath+'\\integrated_histogram.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()

    return True

def plot_esmw(datapath, fit_data = True, save = True):
    plt.close('all')

    suffix = 'esmw'

    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    
    e = np.load(datapath+'\\'+suffix+'-1_statics_and_parameters.npz')
    f_drive = e['mw_drive_freq']
    mwpower = e['mw_power']
    mw_min_freq = e['min_esmw_freq']
    mw_max_freq = e['max_esmw_freq']
    noof_datapoints = e['noof_datapoints']
    e.close()
    
    spin_ro_file = list()
    sp_file = list()
    
    for idx in arange(noof_datapoints):
        for k in files:
            if k == suffix+'-'+num2str(idx,0)+'_Spin_RO.npz':
                spin_ro_file.append(k)
            if k == suffix+'-'+num2str(idx,0)+'_SP_histogram.npz':
                sp_file.append(k)



    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    mw_freq = linspace(mw_min_freq,mw_max_freq,noof_datapoints)
    fit_par=[]
    for idx in arange(noof_datapoints):
        f = np.load(datapath+'\\'+spin_ro_file[idx])
        raw_counts = f['counts']
        repetitions = f['sweep_axis']
        t = f['time']

        tot_size = len(repetitions)
        reps_per_point = tot_size/float(noof_datapoints)

        counts_during_readout = sum(raw_counts, axis = 0)
        ro_time = arange(0,shape(raw_counts)[1])
        figure1 = plt.figure(1)
        
        offset_guess = counts_during_readout.min()
        init_amp_guess = counts_during_readout.max()
        decay_guess = 10

        if fit_data:
            fit_result,p=fit.fit1d(ro_time[arange(4,len(ro_time)-1)], 
                    counts_during_readout[arange(4,len(ro_time)-1)], 
                    common.fit_exp_decay_with_offset, 
                    offset_guess, init_amp_guess, decay_guess,
                    do_plot = True, do_print = True, newfig = False,
                    plot_fitparams_xy = (0.5,0.5),
                    ret=True)
            if fit_result != False:
                fit_par.append(fit_result['params'][2])
            else:
                fit_par.append(0)
        plt.plot(ro_time, counts_during_readout, 'or')
        plt.xlabel('Read-out duration ($\mu$s)')
        plt.ylabel('Integrated counts')
        plt.title('Read-out with MW, driving $f$ ='+num2str(f_drive/1E6,1)+\
                ' MHz, power = '+num2str(mwpower,0)+' dBm')
        if save:
            figure1.savefig(datapath+'\\histogram_integrated'+num2str(idx,0)+'.png')
        plt.clf()

        x = 6.0
        y = 8.0

        figure2 = plt.figure(figsize=(x,y))
        plt.pcolor(raw_counts, cmap = 'hot', antialiased=False)
        plt.xlabel('Readout time (us)')
        plt.ylabel('MW repetition number')
        plt.title('Total histogram, integrated over repetitions')
        plt.colorbar()
        if save:
            figure2.savefig(datapath+'\\histogram_counts_2d'+num2str(idx,0)+'.png')
        plt.clf()

        f.close()
        
        ###########################################
        ######## SPIN PUMPING #####################
        ###########################################
        v = np.load(datapath+'\\'+sp_file[idx])
        sp_counts = v['counts']
        sp_time = v['time']

        offset_guess = sp_counts[len(sp_counts)-1]
        init_amp_guess = sp_counts[2]
        decay_guess = 10

        figure3 = plt.figure(3)
        if fit_data:
            fit.fit1d(sp_time/1E3, sp_counts, common.fit_exp_decay_with_offset, 
                    offset_guess, init_amp_guess, decay_guess,
                    do_plot = True, do_print = True, newfig = False,
                    plot_fitparams_xy = (0.5,0.5))
        
        plt.plot(sp_time/1E3,sp_counts,'sg')
        plt.xlabel('Time ($\mu$s)')
        plt.ylabel('Integrated counts')
        plt.title('Spin pumping')
        v.close()
        if save:
            figure3.savefig(datapath+'\\spin_pumping'+num2str(idx,0)+'.png')
        plt.clf()
    
    if save:
        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: MW Freq (ns)\tCol1: Integrated counts\n'
        col_vals = str()
        for k in arange(len(mw_freq)):
            col_vals += num2str(mw_freq[k],4)+'\t'+num2str(counts_during_readout[k],0)+'\n'
        fo = open(datapath+'\\integrated_histogram.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()
    
    if fit_data:
        plt.close('all')
        print len(mw_freq)
        print len(fit_par)
        figure4 = plt.figure(4)
        plt.plot(mw_freq/1E6,fit_par, '-r')
        plt.xlabel('MW frequency (MHz)')
        plt.ylabel('Decay constant')
        
        if save:
            figure4.savefig(datapath+'\\decay_constant_vs_mw_freq.png')
        plt.clf()



    return fit_par

def plot_esr(datapath, fit_data = True, save = True, f_dip = 2.828E9,msplusone=False,Zsplitting=35e6):
    plt.close('all')
    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    
    for k in files:
        if '.npz' in k:
            data_file = k


    data = np.load(datapath+'\\'+data_file)
    mw_freq = data['freq']
    counts = data['counts']
    data.close()

    f_dip_guess=f_dip
    offset_guess = counts.max()
    dip_depth_guess = offset_guess - counts.min()
    width_guess = 5e-3
    

    if msplusone:
        noof_dips = 2
        dip_separation = Zsplitting
    else:
        noof_dips = 1
        dip_separation=0

    if fit_data:
        figure2 = plt.figure(2)
        figure2.clf()
        fit_res=fit.fit1d(mw_freq/1E9, counts, common.fit_gauss, 
            offset_guess, dip_depth_guess, f_dip_guess/1E9,width_guess,
            do_plot = True, do_print = True, newfig = False)
    
    
    print fit_res
    plt.plot(mw_freq/1E9,counts, '-k')
    plt.xlabel('MW frequency (GHz)')
    plt.ylabel('Integrated counts')
    plt.title('MW frequency sweep')
#, power = '+num2str(mwpower,0)+' dBm)
    #plt.text(0.1*(mw_max_freq+mw_min_freq),max(counts_during_readout),datapath)
    if save:
        figure2.savefig(datapath+'\\esr_data.png')

def plot_Pulse_cal(datapath, fit_data = True, save = True):

    plt.close('all')
    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    
    for k in files:
        if 'statics_and_parameters.npz' in k:
            stats_params_file = k
        if 'Spin_RO.npz' in k:
            spin_ro_file = k
        if 'SP_histogram.npz' in k:
            sp_file = k

    e = np.load(datapath+'\\'+stats_params_file)
    f_drive = e['mw_drive_freq']
    mwpower = e['mw_power']
    min_pulse_nr = e['min_pulse_nr']
    max_pulse_nr = e['max_pulse_nr']
    noof_datapoints = e['noof_datapoints']
    noof_reps = e['completed_repetitions']
    e.close()

    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    f = np.load(datapath+'\\'+spin_ro_file)
    raw_counts = f['counts']
    repetitions = f['sweep_axis']
    SSRO_counts = f['SSRO_counts']
    t = f['time']

    tot_size = len(repetitions)
    reps_per_point = tot_size/float(noof_datapoints)

    idx = 0
    counts_during_readout = zeros(noof_datapoints)
    pulse_nr = linspace(min_pulse_nr,max_pulse_nr,noof_datapoints)
    counts_during_readout = sum(raw_counts, axis = 1)
    SSRO_readout = sum(SSRO_counts, axis = 1)/float(noof_reps)
    

    #########################################
    ############ FITTING ####################
    #########################################
    
    #FIXME to be implemented
    figure2=plt.figure(2)
    figure2.clf()
    plt.plot(pulse_nr,SSRO_readout, 'sk')
    plt.xlabel('Pulse nr')
    plt.ylabel('P ms=0')
    plt.title('MW length sweep, driving $f$ ='+num2str(f_drive/1E6,1)+\
            ' MHz, power = '+num2str(mwpower,0)+' dBm')
    #plt.text(0.1*(mw_max_len+mw_min_len),max(counts_during_readout),datapath)
    if save:
        figure2.savefig(datapath+'\\histogram_integrated.png')

    x = 6.0
    y = 8.0

    figure3 = plt.figure(figsize=(x,y))
    plt.pcolor(raw_counts, cmap = 'hot', antialiased=False)
    plt.xlabel('Readout time (us)')
    plt.ylabel('MW repetition number')
    plt.title('Total histogram, integrated over repetitions')
    plt.colorbar()
    if save:
        figure3.savefig(datapath+'\\histogram_counts_2d.png')

    f.close()
    

    ###########################################
    ######## SPIN PUMPING #####################
    ###########################################
    v = np.load(datapath+'\\'+sp_file)
    sp_counts = v['counts']
    sp_time = v['time']

    offset_guess = sp_counts[len(sp_counts)-1]
    init_amp_guess = sp_counts[2]
    decay_guess = 10

    figure4 = plt.figure(4)
    fit.fit1d(sp_time/1E3, sp_counts, common.fit_exp_decay_with_offset, 
            offset_guess, init_amp_guess, decay_guess,
            do_plot = True, do_print = True, newfig = False,
            plot_fitparams_xy = (0.5,0.5))
    
    plt.plot(sp_time/1E3,sp_counts,'sg')
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Integrated counts')
    plt.title('Spin pumping')
    v.close()
    if save:
        figure4.savefig(datapath+'\\spin_pumping.png')

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: MW length (ns)\tCol1: Integrated counts\n'
        col_vals = str()
        for k in arange(noof_datapoints):
            col_vals += num2str(pulse_nr[k],2)+'\t'+num2str(counts_during_readout[k],0)+'\n'
        fo = open(datapath+'\\integrated_histogram.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()

    return True

def plot_Pulse_cal_amp(datapath, fit_data = True, save = True):

    plt.close('all')
    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    
    for k in files:
        if 'statics_and_parameters.npz' in k:
            stats_params_file = k
        if 'Spin_RO.npz' in k:
            spin_ro_file = k
        if 'SP_histogram.npz' in k:
            sp_file = k

    e = np.load(datapath+'\\'+stats_params_file)
    f_drive = e['mw_drive_freq']
    mwpower = e['mw_power']
    min_amp = e['min_pulse_amp']
    max_amp = e['max_pulse_amp']
    noof_datapoints = e['noof_datapoints']
    noof_reps = e['completed_repetitions']
    e.close()

    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    f = np.load(datapath+'\\'+spin_ro_file)
    raw_counts = f['counts']
    repetitions = f['sweep_axis']
    SSRO_counts = f['SSRO_counts']
    t = f['time']

    tot_size = len(repetitions)
    reps_per_point = tot_size/float(noof_datapoints)

    idx = 0
    counts_during_readout = zeros(noof_datapoints)
    amp = linspace(min_amp,max_amp,noof_datapoints)
    counts_during_readout = sum(raw_counts, axis = 1)
    SSRO_readout = sum(SSRO_counts, axis = 1)/float(noof_reps)

    

    #########################################
    ############ FITTING ####################
    #########################################
    
    #FIXME to be implemented
    figure2=plt.figure(2)
    plt.plot(amp,SSRO_readout, 'sk')
    plt.xlabel('Pulse amp')
    plt.ylabel('P ms=0')
    plt.title('MW length sweep, driving $f$ ='+num2str(f_drive/1E6,1)+\
            ' MHz, power = '+num2str(mwpower,0)+' dBm')
    #plt.text(0.1*(mw_max_len+mw_min_len),max(counts_during_readout),datapath)
    if save:
        figure2.savefig(datapath+'\\histogram_integrated.png')

    x = 6.0
    y = 8.0

    figure3 = plt.figure(figsize=(x,y))
    plt.pcolor(raw_counts, cmap = 'hot', antialiased=False)
    plt.xlabel('Readout time (us)')
    plt.ylabel('MW repetition number')
    plt.title('Total histogram, integrated over repetitions')
    plt.colorbar()
    if save:
        figure3.savefig(datapath+'\\histogram_counts_2d.png')

    f.close()
    

    ###########################################
    ######## SPIN PUMPING #####################
    ###########################################
    v = np.load(datapath+'\\'+sp_file)
    sp_counts = v['counts']
    sp_time = v['time']

    offset_guess = sp_counts[len(sp_counts)-1]
    init_amp_guess = sp_counts[2]
    decay_guess = 10

    figure4 = plt.figure(4)
    fit.fit1d(sp_time/1E3, sp_counts, common.fit_exp_decay_with_offset, 
            offset_guess, init_amp_guess, decay_guess,
            do_plot = True, do_print = True, newfig = False,
            plot_fitparams_xy = (0.5,0.5))
    
    plt.plot(sp_time/1E3,sp_counts,'sg')
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Integrated counts')
    plt.title('Spin pumping')
    v.close()
    if save:
        figure4.savefig(datapath+'\\spin_pumping.png')

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: MW length (ns)\tCol1: Integrated counts\n'
        col_vals = str()
        for k in arange(noof_datapoints):
            col_vals += num2str(amp[k],2)+'\t'+num2str(counts_during_readout[k],0)+'\n'
        fo = open(datapath+'\\integrated_histogram.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()

    return True

def plot_Pulse_cal_time(datapath, fit_data = True, save = True):

    plt.close('all')
    ###########################################
    ######## MEASUREMENT SPECS ################
    ###########################################
    files = os.listdir(datapath)
    
    for k in files:
        if 'statics_and_parameters.npz' in k:
            stats_params_file = k
        if 'Spin_RO.npz' in k:
            spin_ro_file = k
        if 'SP_histogram.npz' in k:
            sp_file = k

    e = np.load(datapath+'\\'+stats_params_file)
    f_drive = e['mw_drive_freq']
    mwpower = e['mw_power']
    min_amp = e['min_pulse_amp']
    max_amp = e['max_pulse_amp']
    min_time = e['min_sweep_par']
    max_time = e['max_sweep_par']
    noof_datapoints = e['noof_datapoints']
    noof_reps = e['completed_repetitions']
    e.close()

    ###########################################
    ######## SPIN RO  #########################
    ###########################################
    
    f = np.load(datapath+'\\'+spin_ro_file)
    raw_counts = f['counts']
    repetitions = f['sweep_axis']
    SSRO_counts = f['SSRO_counts']
    t = f['time']

    tot_size = len(repetitions)
    reps_per_point = tot_size/float(noof_datapoints)

    idx = 0
    counts_during_readout = zeros(noof_datapoints)
    delay_time = linspace(min_time,max_time,noof_datapoints)
    amp = linspace(min_amp,max_amp,noof_datapoints)
    counts_during_readout = sum(raw_counts, axis = 1)
    SSRO_readout = sum(SSRO_counts, axis = 1)/float(noof_reps)

    

    #########################################
    ############ FITTING ####################
    #########################################
    
    #FIXME to be implemented
    figure2=plt.figure(2)
    plt.plot(delay_time,SSRO_readout, 'sk')
    plt.xlabel('time between CORPSE pulses [ns]')
    plt.ylabel('P ms=0')
    plt.title('MW length sweep, driving $f$ ='+num2str(f_drive/1E6,1)+\
            ' MHz, power = '+num2str(mwpower,0)+' dBm')
    #plt.text(0.1*(mw_max_len+mw_min_len),max(counts_during_readout),datapath)
    if save:
        figure2.savefig(datapath+'\\histogram_integrated.png')

    x = 6.0
    y = 8.0

    figure3 = plt.figure(figsize=(x,y))
    plt.pcolor(raw_counts, cmap = 'hot', antialiased=False)
    plt.xlabel('Readout time (us)')
    plt.ylabel('MW repetition number')
    plt.title('Total histogram, integrated over repetitions')
    plt.colorbar()
    if save:
        figure3.savefig(datapath+'\\histogram_counts_2d.png')

    f.close()
    

    ###########################################
    ######## SPIN PUMPING #####################
    ###########################################
    v = np.load(datapath+'\\'+sp_file)
    sp_counts = v['counts']
    sp_time = v['time']

    offset_guess = sp_counts[len(sp_counts)-1]
    init_amp_guess = sp_counts[2]
    decay_guess = 10

    figure4 = plt.figure(4)
    fit.fit1d(sp_time/1E3, sp_counts, common.fit_exp_decay_with_offset, 
            offset_guess, init_amp_guess, decay_guess,
            do_plot = True, do_print = True, newfig = False,
            plot_fitparams_xy = (0.5,0.5))
    
    plt.plot(sp_time/1E3,sp_counts,'sg')
    plt.xlabel('Time ($\mu$s)')
    plt.ylabel('Integrated counts')
    plt.title('Spin pumping')
    v.close()
    if save:
        figure4.savefig(datapath+'\\spin_pumping.png')

        #Save a dat file for use in e.g. Origin with the rabi oscillation.
        curr_date = '#'+time.ctime()+'\n'
        col_names = '#Col0: MW length (ns)\tCol1: Integrated counts\n'
        col_vals = str()
        for k in arange(noof_datapoints):
            col_vals += num2str(amp[k],2)+'\t'+num2str(counts_during_readout[k],0)+'\n'
        fo = open(datapath+'\\integrated_histogram.dat', "w")
        for item in [curr_date, col_names, col_vals]:
            fo.writelines(item)
        fo.close()

    return True
