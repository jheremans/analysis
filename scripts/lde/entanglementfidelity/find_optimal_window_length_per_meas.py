import os

def ro_c_v(corr):
    roc=sscorr.ssro_correct_twoqubit_state_photon_numbers(corr,0.905,0.805, 0.9937,0.998, verbose=False)
    return (roc[0][0]+roc[3][0]),sscorr.get_fidelity_error(corr,0.905,0.805, 0.9937,0.998, verbose=False)
    
def ro_c_v2(corr):
    roc=sscorr.ssro_correct_twoqubit_state_photon_numbers(corr,0.905,0.805, 0.9937,0.998, verbose=False)
    return (roc[1][0]+roc[2][0]),sscorr.get_fidelity_error(corr,0.905,0.805, 0.9937,0.998, verbose=False)

vcorr1=zeros((len(w_lengths),len(w_dts)))
vcorr_psi1=zeros((len(w_lengths),len(w_dts)))
vcorr_psi2=zeros((len(w_lengths),len(w_dts)))
 
for i,w_length in enumerate(w_lengths):
    for j,w_dt in enumerate(w_dts):
        vcorr[i,j]=ro_c_v2(corr[i,j]) if xx else ro_c_v(corr[i,j])
        vcorr_psi1[i,j]=ro_c_v2(corr_00[i,j]+corr_00[i,j]) if xx else ro_c_v(corr_00[i,j]+corr_11[i,j])
        vcorr_psi2[i,j]=ro_c_v2(corr_01[i,j]+corr_10[i,j]) if xx else ro_c_v(corr_01[i,j]+corr_10[i,j])
    
    
if len(w_lengths) == 1 or len(w_dts)==1:
    x= w_lengths if len(w_dts) == 1 else w_dts
    plt.figure()
    ax=plt.subplot(111)
    plt.plot(x,vcorr)
    plt.plot(x,vcorr_psi1, color='g')
    plt.plot(x,vcorr_psi2, color='r')
    ax.set_xlabel('window length [bins=0.256 ns]')
    ax.set_ylabel('01 + 10 probability')
    ax.set_title('RO corrected correlation vs w_length from (637,666)')
else:
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_lengths,vcorr_psi1[:,:])
     plt.legend((w_dts))
     ax.set_xlabel('window length [bins=0.256 ns]')
     ax.set_ylabel('00 + 11 probability')
     ax.set_title('Psi 1 RO corrected correlation vs w_length from (637,666) for different dts')
     
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_lengths,vcorr_psi2[:,:])
     plt.legend((w_dts))
     ax.set_xlabel('window length [bins=0.256 ns]')
     ax.set_ylabel('01 + 10 probability')
     ax.set_title('Psi 2 RO corrected correlation vs w_length from (637,666) for different dts')
     
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_dts,transpose(vcorr_psi1[:,:]))
     plt.legend((w_lengths))
     ax.set_xlabel('dt window [bins=0.256 ns]')
     ax.set_ylabel('00 + 11 probability')
     ax.set_title('Psi 1 RO corrected correlation vs dts from (637,666) for different w_length')
     
     plt.figure()
     ax=plt.subplot(111)
     plt.plot(w_dts,transpose(vcorr_psi2[:,:]))
     plt.legend((w_lengths))
     ax.set_xlabel(' dt window  [bins=0.256 ns]')
     ax.set_ylabel('01 + 10 probability')
     ax.set_title('Psi 2 RO corrected correlation vs dts from (637,666) vor different w_length')
    
#ax.

imax=where(vcorr_psi1==amax(vcorr_psi1))[0][0]
jmax=where(vcorr_psi1==amax(vcorr_psi1))[1][0]
print 'vcorr_psi1 max:', amax(vcorr_psi1), \
        w_lengths[imax], w_dts[jmax], corr_00[imax,jmax].sum()+corr_11[imax,jmax].sum()
imax=where(vcorr_psi2==amax(vcorr_psi2))[0][0]
jmax=where(vcorr_psi2==amax(vcorr_psi2))[1][0]
print 'vcorr_psi2 max:', amax(vcorr_psi2), \
        w_lengths[imax], w_dts[jmax], corr_01[imax,jmax].sum()+corr_10[imax,jmax].sum()
imax=where(((vcorr_psi1+vcorr_psi2)/2.)==amax(((vcorr_psi1+vcorr_psi2)/2.)))[0][0]
jmax=where(((vcorr_psi1+vcorr_psi2)/2.)==amax(((vcorr_psi1+vcorr_psi2)/2.)))[1][0]
print 'vcorr_psi1+2/2 max:', amax(((vcorr_psi1+vcorr_psi2)/2.)), \
        w_lengths[imax], w_dts[jmax], corr[imax,jmax].sum()  
print 'vcorr max:',amax(vcorr), w_lengths[where(vcorr==amax(vcorr))[0][0]], w_dts[where(vcorr==amax(vcorr))[1][0]] 
