
import numpy as np
import pylab as plt
from  matplotlib import rc
from analysis.lib.spin import spin_control as sc

def conditioned_data(s_array,f_array):
    sarray=np.load(s_array)
    seg_array=sarray['segmented_RO_data']
    farray=np.load(f_array)
    full_RO_array=farray['full_RO_data']
    bins = max(seg_array)+1
    print "Number of segments: " + str(bins)
    j=0
    corr_array=np.zeros(bins)
    corr_norm_array=np.zeros(bins)
    corrnormalized=np.zeros(bins)
    for i in seg_array:
        corr_array[i]+=full_RO_array[j]
        corr_norm_array[i]+=1
        j+=1
    corrnormalized = corr_array/(corr_norm_array+0.)
    return corrnormalized,corr_array,corr_norm_array
'''

for n in np.arange(6):
        sdata=r'D:\measuring\data\20130628\140659_Adwin_SSRO__LT2_Segmented_RO_sweep_power_conditioned_data_0_2nW_noClickRecorded\Adwin_SSRO-00'+str(n)+'_segmented_RO_data.npz'
        fdata=r'D:\measuring\data\20130628\140659_Adwin_SSRO__LT2_Segmented_RO_sweep_power_conditioned_data_0_2nW_noClickRecorded\Adwin_SSRO-00'+str(n)+'_full_RO_data.npz' 
   
        d,corr,cornorm=conditioned_data(sdata,fdata)
        fname = 'C:\Documents and Settings\localadmin\Desktop\correlData\\' +str(n)+'.npz' 
        np.savez (fname, corr= corr, norm = cornorm, norm_correlation = d)
 
plt.figure()
plt.plot (d)
plt.ylim([0,1])
plt.ylabel ('read-out fidelity|photon detected')
plt.xlabel ('Segment number')
plt.show()
'''

def average_RO_fid(folder):

    ms_zero_seg=np.load (folder+'Adwin_SSRO-000_segmented_RO_data.npz')
    ms_zero_seg_fid=1-list(ms_zero_seg['segmented_RO_data']).count(max(ms_zero_seg['segmented_RO_data']))/float(len(ms_zero_seg['segmented_RO_data']))
    
    ms_zero_seg.close()
    ms_one_seg=np.load (folder+'Adwin_SSRO-001_segmented_RO_data.npz')
    ms_one_seg_fid=list(ms_one_seg['segmented_RO_data']).count(max(ms_one_seg['segmented_RO_data']))/float(len(ms_one_seg['segmented_RO_data']))
    print len(ms_one_seg['segmented_RO_data'])
    print max(ms_one_seg['segmented_RO_data'])
    ms_one_seg.close()
 
    fig=plt.figure()
    
  
    plt.bar(0,ms_one_seg_fid,width=1,color='Crimson',alpha=0.8,linewidth=1.5,edgecolor='DarkRed')
    plt.bar(0,1-ms_zero_seg_fid,width=1,color='RoyalBlue',alpha=0.8,linewidth=1.5,edgecolor='b')
    plt.bar(1,ms_zero_seg_fid,width=1,color='RoyalBlue',alpha=0.8,linewidth=1.5,edgecolor='b')
    plt.bar(1,1-ms_one_seg_fid,width=1,color='Crimson',alpha=0.8,linewidth=1.5,edgecolor='DarkRed')
    s='Average Fid: \n' + str((ms_zero_seg_fid+ms_one_seg_fid)/2.)+'\nms=0 :   '+str(int(ms_zero_seg_fid*1000)/1000.)+'\nms=-1:   '+str(int(ms_one_seg_fid*1000)/1000.)
    plt.text(2.2,0.5,s)
    plt.ylim([0,1])
    plt.xlim([0,3])
    plt.xticks([0.5,1.5],['0','1'])
    plt.ylabel ('Probability')
    plt.xlabel ('nr of photons')
    plt.legend(['ms=-1','ms=0'])
    plt.show()
    ms_zero=np.load (folder+'Adwin_SSRO-000_segment_number.npz')
    ms_one=np.load (folder+'Adwin_SSRO-001_segment_number.npz')
    print sum(ms_zero['segment_number'])/5000.
    print 1-sum(ms_one['segment_number'])/5000.
    ms_zero.close()
    ms_one.close()

def segmented_N_ramsey(nr_of_datapoints=1,folder = r'D:\measuring\data\20130709\123046_Seg_RO__LT2_N_Ramsey_seg_RO_100us_750pW/'):
    a= np.load (folder+'Seg_RO-000_segment_number.npz')
    seg_nr=a['segment_number']
    a.close()
    
    a= np.load (folder+'Seg_RO-000_segment_number.npz')
    seg_nr=a['segment_number']
    a.close()
    
    c= np.load (folder+'Seg_RO-000_Spin_RO.npz')
    SSRO_counts=c['SSRO_counts']
    cond_RO_data=c['cond_RO_data']
    x=c['sweep_axis']
    c.close()
    
    nr_of_datapoints=len(SSRO_counts)
    reps=int(sum(seg_nr)/float(nr_of_datapoints))
    total_segments=len(seg_nr)
    #total_segments=45
    
    plt.figure()
    phase=[]
    segnr_phase=[]
    reconstr_sum=np.zeros(201)
    corrected=np.zeros(201)
    accum_seg_nr=np.zeros(len(seg_nr))
    for i in np.arange(total_segments):
        if i in [4,5,6,7,8]:
            do_plot=True
            print i
        else:
            do_plot=False
        y=cond_RO_data[nr_of_datapoints*i:nr_of_datapoints*(i+1)]
        #FIX: This normalization is very rough, we should devide each point with the number of repetitions (get it from segmented_RO_data.npz)
        y_norm=y/float(max(y))
        
        dict=sc.fit_sin(x,y,sqrt(max(y)),fixed=[0],fix_param=[1/360.],do_plot=do_plot)
        if dict:
            phase.append(dict['params'][2])
            segnr_phase.append(i)
            reconstr_sum=reconstr_sum+dict['fitfunc'](np.linspace(0,360,201))    
            A=dict['params'][0]
            a=dict['params'][1]
            phi=dict['params'][2]
            f=1/360.
            corrected=corrected+fit_func(abs(a),abs(A),f,100,np.linspace(0,360,201))
            accum_seg_nr[i]=sum(seg_nr[0:i])
   
    plt.show()    
    plt.figure()
    plt.errorbar(segnr_phase,phase,dict['error_dict']['phi'],fmt='o')
    plt.xlabel('Segment number of photon click',fontsize=14)
    plt.ylabel('Phase second RF pulse',fontsize=14)
    plt.yticks([-90,-45,0,45,90])
    plt.xticks([0,25,50,75,100])
    #plt.ylim([-180,180])
    plt.show()
    
    plt.clf()
    plt.show()    
    plt.figure()
    plt.plot(np.linspace(0,360,201),reconstr_sum/float(reps))
    plt.plot(np.linspace(0,360,201),corrected/float(reps))
    plt.xlabel('Phase',fontsize=14)
    plt.ylabel('P(ms=0) reconstructed from separate fits',fontsize=14)
    plt.legend(['Total','Corrected for phase shift'],loc=4)
    #plt.yticks([0,45,90,135,180])
    #plt.xticks([0,10,20,30,40])
    #plt.ylim([0.2,0.5])
    plt.show()
    
    plt.figure()
    plt.plot(np.arange(total_segments-1)+1,seg_nr[0:total_segments-1]/float(sum(seg_nr)),'o')
    plt.xlabel('Segment number',fontsize=14)
    plt.ylabel('P(click)',fontsize=14)
    plt.show()
      
    plt.figure()
    plt.plot(np.arange(total_segments-1)+1,accum_seg_nr[0:total_segments-1]/float(sum(seg_nr)),'o')
    plt.xlabel('Segment number',fontsize=14)
    plt.ylabel('P(click) accummulated',fontsize=14)
    plt.show()
    
    print 'Contrast for total:'
    print (max(reconstr_sum/float(reps))-min(reconstr_sum/float(reps)))/0.5
    print 'Contrast for corrected:'
    print (max(corrected/float(reps))-min(corrected/float(reps)))/0.5
    print dict['error_dict']['phi']
    
    return SSRO_counts, cond_RO_data, seg_nr
def fit_func(a,A,f,phi,x):
    return (a + A * cos(2.*np.pi*(f*x + phi/360.)))
def segmented_RO_analysis(nr_of_datapoints=1,folder = r'D:\measuring\data\20130703\144857_Adwin_SSRO__LT2_SIL10_Segmented_RO_sweep_power/'):
    analysis = True

    Ndet = np.zeros(nr_of_datapoints)
    fid1 = np.zeros(nr_of_datapoints)
    fid2 = np.zeros(nr_of_datapoints)
    power = np.zeros(nr_of_datapoints)
    QND = np.zeros(nr_of_datapoints)
    for n in (np.arange(nr_of_datapoints)):
        if (n<10):
            nn = '0'+str(n)
        else:
            nn=str(n)
        #if analysis:
        #	ssro.run_single (folder=fold, index = n)
            
        a = np.load (folder+'Seg_RO-0'+nn+'_segment_number.npz')
        par=np.load (folder+'Seg_RO-0'+nn+'_parameters_dict.npz')
        Ndet [n] = np.sum(a['segment_number'])
        b = np.load (folder+'Seg_RO-0'+nn+'_segmented_RO_data.npz')
        c = np.load (folder+'Seg_RO-0'+nn+'_full_RO_data.npz')
        
        fid1[n]=1-list(b['segmented_RO_data']).count(max(b['segmented_RO_data']))/float(len(b['segmented_RO_data']))
        fid2[n]=sum(c['full_RO_data'])/float(len(b['segmented_RO_data']))
        QND[n]=fid1[n]*fid2[n]/0.854
        power[n]=par['segmented_Ex_RO_amplitude']*1e9
        a.close()
        b.close()
        c.close()
        
        
    #print max(fid1)
    #print fid1*fid2
    plt.figure()
    plt.plot (power, Ndet/50000., 'o')
    plt.ylim([0,1])
    plt.ylabel ('read-out fidelity')
    plt.xlabel ('read-out power [nW]')
    plt.show()
    '''
    plt.figure()
    
    plt.ylim([0,1])
    plt.ylabel ('fidelity segmented RO')
    plt.xlabel ('read-out power [nW]')
    plt.show()
    '''
    print QND
    plt.figure()
    plt.plot (power, QND, 'o')
    plt.ylim([0,1])
    plt.ylabel ('QND fid')
    plt.xlabel ('read-out power [nW]')
    plt.show()
    
    plt.figure()
    plt.plot (power, fid1,'o')
    plt.plot (power, fid2,'o')
    plt.ylim([0,1])
    plt.ylabel ('P(click)')
    plt.xlabel ('read-out power [nW]')
    plt.legend(['Segmented RO','Second RO'],loc=4)
    plt.show()
    
    
def plot_dephasing ():
    
    pRO_time = [0,4, 6, 8, 10, 15, 25, 50, 75, 100]
    segmRO_A = np.array([0.346, 0.244136, 0.19500, 0.160282, 0.13989, 0.11728, 0.09146, 0.09296, 0.09089, 0.09178])
    segmRO_dA = np.array([0.0085, 0.0066, 0.0026, 0.0051, 0.0054, 0.0048, 0.0044, 0.0026, 0.0053, 0.0049])
    fullRO_A = np.array([0.0085, 0.2514, 0.1686, 0.1163, 0.0854, 0.0431, 0.027, 0.010, 0.008, 0.005])
    fullRO_dA = np.array ([0.0085, 0.0064, 0.0024, 0.004, 0.0065, 0.0036, 0.006, 0.005, 0.004, 0.006]) 
    RO_time = [4, 6, 8, 10, 15, 25, 50, 75, 100]
    sRO_ampl = [0.53, 0.52, 0.52, 0.54, 0.51, 0.5, 0.51,0.5, .48]
    fRO_ampl = [0.52,0.51, 0.52, 0.51, 0.53, 0.53, 0.53, 0.53, 0.52]
    usRO_ampl=np.ones(9)*0.013
    ufRO_ampl=np.ones(9)*0.013
    colors={
        'fs':'RoyalBlue',
        'fn':'Crimson'
        }
    plt.errorbar (pRO_time, segmRO_A/0.346, segmRO_dA, label='Segmented RO',marker='o', linestyle='None',mfc=colors['fn'], mec=colors['fn'],ecolor=colors['fn'],color=colors['fn'])
    plt.errorbar (pRO_time, fullRO_A/0.346, fullRO_dA, label='Normal RO',marker='o', linestyle='None',mfc=colors['fs'], mec=colors['fs'],ecolor=colors['fs'],color=colors['fs'])
    plt.ylim ([0, 1.05])
    plt.xlim ([0, 110])
    plt.ylabel ('normalized phase contrast', fontsize=16)
    plt.xlabel ('read-out time [us]', fontsize=16)
    plt.legend()
    plt.show()
   

    plt.errorbar (RO_time, sRO_ampl,usRO_ampl, label='Segmented RO',marker='o', linestyle='None',mfc=colors['fn'], mec=colors['fn'],ecolor=colors['fn'],color=colors['fn'])
    plt.errorbar (RO_time, fRO_ampl,ufRO_ampl, label='Normal RO',marker='o', linestyle='None',mfc=colors['fs'], mec=colors['fs'],ecolor=colors['fs'],color=colors['fs'])
    plt.ylim ([0, 1])
    plt.xlim ([0, 110])
    plt.ylabel ('Z', fontsize=16)
    plt.xlabel ('read-out time [us]', fontsize=16)
    plt.legend()
    plt.show()

plot_dephasing()    