<<<<<<< HEAD
import os
from lde import fidelities
        
npz_file_ZZ=r'D:\measuring\analysis\data\lde\ZZ\all_ro_cor_fast.npz'
npz_file_XX=r'D:\measuring\analysis\data\lde\XX\all_ro_cor_fast.npz'
npz_file_XmX=r'D:\measuring\analysis\data\lde\X-X\all_ro_cor_fast.npz'
dZZ=load(npz_file_ZZ)
dXX=load(npz_file_XX)
dXmX=load(npz_file_XmX)
ZZ_corr=dZZ['corr']
ZZ_corr_00=dZZ['corr_00']
ZZ_corr_01=dZZ['corr_01']
ZZ_corr_10=dZZ['corr_10']
ZZ_corr_11=dZZ['corr_11']
XX_corr=dXX['corr']
XX_corr_00=dXX['corr_00']
XX_corr_01=dXX['corr_01']
XX_corr_10=dXX['corr_10']
XX_corr_11=dXX['corr_11']
XmX_corr=dXmX['corr']
XmX_corr_00=dXmX['corr_00']
XmX_corr_01=dXmX['corr_01']
XmX_corr_10=dXmX['corr_10']
XmX_corr_11=dXmX['corr_11']

w_lengths=dZZ['w_lengths']
w_dts=dZZ['w_dts']
        
F_psi1=zeros((len(w_lengths),len(w_dts)))
dF_psi1=zeros((len(w_lengths),len(w_dts)))
F_psi2=zeros((len(w_lengths),len(w_dts)))
dF_psi2=zeros((len(w_lengths),len(w_dts)))

dZZ_psi1=zeros((len(w_lengths),len(w_dts)))
dXX_psi1=zeros((len(w_lengths),len(w_dts)))
dZZ_psi2=zeros((len(w_lengths),len(w_dts)))
dXX_psi2=zeros((len(w_lengths),len(w_dts)))

ro_correct=True

for i,w_length in enumerate(w_lengths):
    for j,w_dt in enumerate(w_dts):
        
        ZZ_corr_psi1=ZZ_corr_00[i,j] + ZZ_corr_11[i,j]
        XX_corr_psi1=XX_corr_00[i,j] + XX_corr_11[i,j]
        XmX_corr_psi1=XmX_corr_00[i,j]+XmX_corr_11[i,j]
        
        F_psi1[i,j],dF_psi1[i,j],_dict = \
                fidelities.get_fidelity(ZZ_corr_psi1,XX_corr_psi1,XmX_corr_psi1, ro_correct=ro_correct, psi1=True)
        dZZ_psi1[i,j]=1/4.*_dict['dZZ']**2 + _dict['dZZS']**2 
        dXX_psi1[i,j]=_dict['dXXavg']**2

        ZZ_corr_psi2=ZZ_corr_01[i,j]+ ZZ_corr_10[i,j]
        XX_corr_psi2=XX_corr_01[i,j]+ XX_corr_10[i,j]
        XmX_corr_psi2=XmX_corr_01[i,j]+XmX_corr_10[i,j]
        
        F_psi2[i,j],dF_psi2[i,j],_dict = \
                fidelities.get_fidelity(ZZ_corr_psi2,XX_corr_psi2,XmX_corr_psi2, ro_correct=ro_correct, psi1=False)
        dZZ_psi2[i,j]=1/4.*_dict['dZZ']**2 + _dict['dZZS']**2 
        dXX_psi2[i,j]=_dict['dXXavg']**2

if ro_correct:
    F2stdev_psi1=(F_psi1-0.5*ones((len(w_lengths),len(w_dts))))/dF_psi1
    F2stdev_psi2=(F_psi2-0.5*ones((len(w_lengths),len(w_dts))))/dF_psi2

do_plot=True
save_plots=True
do_plot_errors=True
save_path=r'H:\My Documents\processed data\lde'
w_dt_bins=[0,1,2,4,8]
w_length_bins=[0,1,2,3,4]
if do_plot:
     fig=plt.figure()
     plt.suptitle('Psi 1 RO corrected Fidelity from (637,666)')
     ax=plt.subplot(221)
     #plt.plot(w_lengths,F_psi1[:,w_dt_bins])
     for j in w_dt_bins:
        #print j
        plt.errorbar(w_lengths,F_psi1[:,j], yerr=dF_psi1[:,j])
     plt.legend((w_dts[w_dt_bins]), loc=4)
     #ax.set_xlabel('window length [bins=0.256 ns]')
     ax.set_ylabel('Fidelity')
     ax.set_ylim([0,1])
     ax.set_title('Vs w_length for different dts')

     ax=plt.subplot(222)
     #plt.plot(w_dts,transpose(F_psi1[w_length_bins,:]))
     for i in w_length_bins:
        #print j
        plt.errorbar(w_dts,F_psi1[i,:], yerr=dF_psi1[i,:])
     plt.legend((w_lengths[w_length_bins]), loc=4)
     #ax.set_xlabel('dt window [bins=0.256 ns]')
     ax.set_ylabel('Fidelity')
     ax.set_ylim([0,1])
     ax.set_title('vs dts for different w_length')
     
     if do_plot_errors and ro_correct:

        ax=plt.subplot(223)
        plt.plot(w_lengths,F2stdev_psi1[:,w_dt_bins])
        plt.legend((w_dts[w_dt_bins]), loc=4)
        ax.set_xlabel('window length [bins=0.256 ns]')
        ax.set_ylabel('(Fidelity-0.5)/sigma ')
        #ax.set_title('vs w_length for different dts')
          
        ax=plt.subplot(224)
        plt.plot(w_dts,transpose(F2stdev_psi1[w_length_bins,:]))
        plt.legend((w_lengths[w_length_bins]), loc=4)
        ax.set_xlabel('dt window [bins=0.256 ns]')
        ax.set_ylabel('(Fidelity-0.5)/sigma ')
        #ax.set_title('vs dts for different w_length')
     
     if save_plots:
        fig.savefig(save_path+'_psi1.pdf')     
if do_plot:
     fig=plt.figure()
     plt.suptitle('Psi 2 RO corrected Fidelity from (637,666) ' )
     
     ax=plt.subplot(221)
     #plt.plot(w_lengths,F_psi2[:,w_dt_bins])
     for j in w_dt_bins:
        #print j
        plt.errorbar(w_lengths,F_psi2[:,j], yerr=dF_psi2[:,j])
     plt.legend((w_dts[w_dt_bins]), loc=4)
     #ax.set_xlabel('window length [bins=0.256 ns]')
     ax.set_ylabel('Fidelity')
     ax.set_ylim([0,1])
     ax.set_title('vs w_length for different dts')
     
     ax=plt.subplot(222)
     #plt.plot(w_dts,transpose(F_psi2[w_length_bins,:]))
     for i in w_length_bins:
        #print j
        plt.errorbar(w_dts,F_psi2[i,:], yerr=dF_psi2[i,:])
     plt.legend((w_lengths[w_length_bins]), loc=4)
     #ax.set_xlabel(' dt window  [bins=0.256 ns]')
     ax.set_ylabel('Fidelity')
     ax.set_ylim([0,1])
     ax.set_title('vs dts vor different w_length')
     
     
     #-----------------------------------------------
     
     if do_plot_errors and ro_correct:
     
         ax=plt.subplot(223)
         plt.plot(w_lengths,F2stdev_psi2[:,w_dt_bins])
         plt.legend((w_dts[w_dt_bins]), loc=4)
         ax.set_xlabel('window length [bins=0.256 ns]')
         ax.set_ylabel('(Fidelity-0.5)/sigma ')
         #ax.set_title('vs w_length for different dts')
         
         ax=plt.subplot(224)
         plt.plot(w_dts,transpose(F2stdev_psi2[w_length_bins,:]))
         plt.legend((w_lengths[w_length_bins]), loc=4)
         ax.set_xlabel(' dt window  [bins=0.256 ns]')
         ax.set_ylabel('(Fidelity-0.5)/sigma ')
         #ax.set_title('vs dts for different w_length')     
     
         

        
         #######
         
     if save_plots:
        fig.savefig(save_path+'_psi2.pdf')     
     
plot_XX_XX_errors=False
if plot_XX_XX_errors:
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_lengths,dZZ_psi2[:,:])
     plt.legend((w_dts))
     ax.set_xlabel('window length [bins=0.256 ns]')
     ax.set_ylabel('Fidelity-2*sigma ')
     ax.set_title('Psi 1 RO corrected Fidelity-2*sigma vs w_length from (637,666) for different dts')
     
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_lengths,dXX_psi2[:,:])
     plt.legend((w_dts))
     ax.set_xlabel('window length [bins=0.256 ns]')
     ax.set_ylabel('Fidelity-2*sigma ')
     ax.set_title('Psi 2 RO corrected Fidelity-2*sigma  vs w_length from (637,666) for different dts')
     
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_dts,transpose(dZZ_psi2[:,:]))
     plt.legend((w_lengths))
     ax.set_xlabel('dt window [bins=0.256 ns]')
     ax.set_ylabel('Fidelity-2*sigma ')
     ax.set_title('Psi 1 RO corrected Fidelity-2*sigma  vs dts from (637,666) for different w_length')
     
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_dts,transpose(dXX_psi2[:,:]))
     plt.legend((w_lengths))
     ax.set_xlabel(' dt window  [bins=0.256 ns]')
     ax.set_ylabel('Fidelity-2*sigma ')
     ax.set_title('Psi 2 RO corrected Fidelity-2*sigma  vs dts from (637,666) vor different w_length')     
     
         
=======
import os
from analysis.lib.lde import fidelities
        
npz_file_ZZ=r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde/ZZ/all_ro_cor_fast.npz'
npz_file_XX=r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde/XX/all_ro_cor_fast.npz'
npz_file_XmX=r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde/X-X/all_ro_cor_fast.npz'
dZZ=load(npz_file_ZZ)
dXX=load(npz_file_XX)
dXmX=load(npz_file_XmX)
ZZ_corr=dZZ['corr']
ZZ_corr_00=dZZ['corr_00']
ZZ_corr_01=dZZ['corr_01']
ZZ_corr_10=dZZ['corr_10']
ZZ_corr_11=dZZ['corr_11']
XX_corr=dXX['corr']
XX_corr_00=dXX['corr_00']
XX_corr_01=dXX['corr_01']
XX_corr_10=dXX['corr_10']
XX_corr_11=dXX['corr_11']
XmX_corr=dXmX['corr']
XmX_corr_00=dXmX['corr_00']
XmX_corr_01=dXmX['corr_01']
XmX_corr_10=dXmX['corr_10']
XmX_corr_11=dXmX['corr_11']

w_lengths=dZZ['w_lengths']
w_dts=dZZ['w_dts']
        
F_psi1=zeros((len(w_lengths),len(w_dts)))
dF_psi1=zeros((len(w_lengths),len(w_dts)))
F_psi2=zeros((len(w_lengths),len(w_dts)))
dF_psi2=zeros((len(w_lengths),len(w_dts)))

dZZ_psi1=zeros((len(w_lengths),len(w_dts)))
dXX_psi1=zeros((len(w_lengths),len(w_dts)))
dZZ_psi2=zeros((len(w_lengths),len(w_dts)))
dXX_psi2=zeros((len(w_lengths),len(w_dts)))

ro_correct=True

for i,w_length in enumerate(w_lengths):
    for j,w_dt in enumerate(w_dts):
        
        ZZ_corr_psi1=ZZ_corr_00[i,j] + ZZ_corr_11[i,j]
        XX_corr_psi1=XX_corr_00[i,j] + XX_corr_11[i,j]
        XmX_corr_psi1=XmX_corr_00[i,j]+XmX_corr_11[i,j]
        
        F_psi1[i,j],dF_psi1[i,j],_dict = \
                fidelities.get_fidelity(ZZ_corr_psi1,XX_corr_psi1,XmX_corr_psi1, ro_correct=ro_correct, psi1=True)
        dZZ_psi1[i,j]=1/4.*_dict['dZZ']**2 + _dict['dZZS']**2 
        dXX_psi1[i,j]=_dict['dXXavg']**2

        ZZ_corr_psi2=ZZ_corr_01[i,j]+ ZZ_corr_10[i,j]
        XX_corr_psi2=XX_corr_01[i,j]+ XX_corr_10[i,j]
        XmX_corr_psi2=XmX_corr_01[i,j]+XmX_corr_10[i,j]
        
        F_psi2[i,j],dF_psi2[i,j],_dict = \
                fidelities.get_fidelity(ZZ_corr_psi2,XX_corr_psi2,XmX_corr_psi2, ro_correct=ro_correct, psi1=False)
        dZZ_psi2[i,j]=1/4.*_dict['dZZ']**2 + _dict['dZZS']**2 
        dXX_psi2[i,j]=_dict['dXXavg']**2

if ro_correct:
    F2stdev_psi1=(F_psi1-0.5*ones((len(w_lengths),len(w_dts))))/dF_psi1
    F2stdev_psi2=(F_psi2-0.5*ones((len(w_lengths),len(w_dts))))/dF_psi2

do_plot=True
save_plots=True
do_plot_errors=True
save_path=r'H:\My Documents\processed data\lde'
w_dt_bins=[0,1,2,4,8]
w_length_bins=[0,1,2,3,4]
if do_plot:
     fig=plt.figure()
     plt.suptitle('Psi 1 RO corrected Fidelity from (637,666)')
     ax=plt.subplot(221)
     #plt.plot(w_lengths,F_psi1[:,w_dt_bins])
     for j in w_dt_bins:
        #print j
        plt.errorbar(w_lengths,F_psi1[:,j], yerr=dF_psi1[:,j])
     plt.legend((w_dts[w_dt_bins]), loc=4)
     #ax.set_xlabel('window length [bins=0.256 ns]')
     ax.set_ylabel('Fidelity')
     ax.set_ylim([0,1])
     ax.set_title('Vs w_length for different dts')

     ax=plt.subplot(222)
     #plt.plot(w_dts,transpose(F_psi1[w_length_bins,:]))
     for i in w_length_bins:
        #print j
        plt.errorbar(w_dts,F_psi1[i,:], yerr=dF_psi1[i,:])
     plt.legend((w_lengths[w_length_bins]), loc=4)
     #ax.set_xlabel('dt window [bins=0.256 ns]')
     ax.set_ylabel('Fidelity')
     ax.set_ylim([0,1])
     ax.set_title('vs dts for different w_length')
     
     if do_plot_errors and ro_correct:

        ax=plt.subplot(223)
        plt.plot(w_lengths,F2stdev_psi1[:,w_dt_bins])
        plt.legend((w_dts[w_dt_bins]), loc=4)
        ax.set_xlabel('window length [bins=0.256 ns]')
        ax.set_ylabel('(Fidelity-0.5)/sigma ')
        #ax.set_title('vs w_length for different dts')
          
        ax=plt.subplot(224)
        plt.plot(w_dts,transpose(F2stdev_psi1[w_length_bins,:]))
        plt.legend((w_lengths[w_length_bins]), loc=4)
        ax.set_xlabel('dt window [bins=0.256 ns]')
        ax.set_ylabel('(Fidelity-0.5)/sigma ')
        #ax.set_title('vs dts for different w_length')
     
     if save_plots:
        fig.savefig(save_path+'_psi1.pdf')     
if do_plot:
     fig=plt.figure()
     plt.suptitle('Psi 2 RO corrected Fidelity from (637,666) ' )
     
     ax=plt.subplot(221)
     #plt.plot(w_lengths,F_psi2[:,w_dt_bins])
     for j in w_dt_bins:
        #print j
        plt.errorbar(w_lengths,F_psi2[:,j], yerr=dF_psi2[:,j])
     plt.legend((w_dts[w_dt_bins]), loc=4)
     #ax.set_xlabel('window length [bins=0.256 ns]')
     ax.set_ylabel('Fidelity')
     ax.set_ylim([0,1])
     ax.set_title('vs w_length for different dts')
     
     ax=plt.subplot(222)
     #plt.plot(w_dts,transpose(F_psi2[w_length_bins,:]))
     for i in w_length_bins:
        #print j
        plt.errorbar(w_dts,F_psi2[i,:], yerr=dF_psi2[i,:])
     plt.legend((w_lengths[w_length_bins]), loc=4)
     #ax.set_xlabel(' dt window  [bins=0.256 ns]')
     ax.set_ylabel('Fidelity')
     ax.set_ylim([0,1])
     ax.set_title('vs dts vor different w_length')
     
     
     #-----------------------------------------------
     
     if do_plot_errors and ro_correct:
     
         ax=plt.subplot(223)
         plt.plot(w_lengths,F2stdev_psi2[:,w_dt_bins])
         plt.legend((w_dts[w_dt_bins]), loc=4)
         ax.set_xlabel('window length [bins=0.256 ns]')
         ax.set_ylabel('(Fidelity-0.5)/sigma ')
         #ax.set_title('vs w_length for different dts')
         
         ax=plt.subplot(224)
         plt.plot(w_dts,transpose(F2stdev_psi2[w_length_bins,:]))
         plt.legend((w_lengths[w_length_bins]), loc=4)
         ax.set_xlabel(' dt window  [bins=0.256 ns]')
         ax.set_ylabel('(Fidelity-0.5)/sigma ')
         #ax.set_title('vs dts for different w_length')     
     
         

        
         #######
         
     if save_plots:
        fig.savefig(save_path+'_psi2.pdf')     
     
plot_XX_XX_errors=False
if plot_XX_XX_errors:
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_lengths,dZZ_psi2[:,:])
     plt.legend((w_dts))
     ax.set_xlabel('window length [bins=0.256 ns]')
     ax.set_ylabel('Fidelity-2*sigma ')
     ax.set_title('Psi 1 RO corrected Fidelity-2*sigma vs w_length from (637,666) for different dts')
     
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_lengths,dXX_psi2[:,:])
     plt.legend((w_dts))
     ax.set_xlabel('window length [bins=0.256 ns]')
     ax.set_ylabel('Fidelity-2*sigma ')
     ax.set_title('Psi 2 RO corrected Fidelity-2*sigma  vs w_length from (637,666) for different dts')
     
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_dts,transpose(dZZ_psi2[:,:]))
     plt.legend((w_lengths))
     ax.set_xlabel('dt window [bins=0.256 ns]')
     ax.set_ylabel('Fidelity-2*sigma ')
     ax.set_title('Psi 1 RO corrected Fidelity-2*sigma  vs dts from (637,666) for different w_length')
     
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_dts,transpose(dXX_psi2[:,:]))
     plt.legend((w_lengths))
     ax.set_xlabel(' dt window  [bins=0.256 ns]')
     ax.set_ylabel('Fidelity-2*sigma ')
     ax.set_title('Psi 2 RO corrected Fidelity-2*sigma  vs dts from (637,666) vor different w_length')     
     
         
>>>>>>> 998bfb9d754ee59f58c639b71f91f58a0a5b6921
