from analysis.lib.lde import tail_cts_per_shot_v4
from analysis.lib.fitting import rabi

path_20 = r'D:\measuring\data\LDE\analysis_data\opt_rabi_vs_CR\2012-10-24-CRvsRabi_LT1\183600_rabi_vs_cr\20nW'

path_60 = r'D:\measuring\data\LDE\analysis_data\opt_rabi_vs_CR\2012-10-24-CRvsRabi_LT1\183600_rabi_vs_cr\60nW'

path_250 = r'D:\measuring\data\LDE\analysis_data\opt_rabi_vs_CR\2012-10-24-CRvsRabi_LT1\183600_rabi_vs_cr\250nW'

output_folder = r'D:\measuring\data\LDE\analysis_output\20121104-opt-rabi-vs-CR'


paths = [path_20, path_60, path_250]

ave_fR = []
ave_tau = []

all_fR = np.zeros((0,14))
all_fR_err = np.zeros((0,14))
all_tau = np.zeros((0,14))
all_tau_err = np.zeros((0,14))

for path in paths:

    files = [ f for f in os.listdir(path) if f[:1] != '.' and f[-4:] == '.txt' ]

    tau_list=[]
    tau_err_list=[]
    fR_list=[]
    fR_err_list=[]
    CR=[]

    for i in files:
        print i
        fR, fR_err, tau, tau_err = tail_cts_per_shot_v4.optical_rabi_resonant_fit_hannes(
                os.path.join(path,i),34,61, use_channel=0, binsize=0.128)
        tau_list.append(tau)
        tau_err_list.append(tau_err)
        fR_list.append(fR*1e3)
        fR_err_list.append(fR_err*1e3)
        CR.append(float(i[:2]))


    fig=plt.figure()
    plt.errorbar(CR,tau_list,yerr=tau_err_list,fmt='ro')
    plt.xlabel('CR threshold')
    plt.ylabel('Decay time [ns]')
    plt.title(path)
    fig.savefig(os.path.join(path,'CRvsTau.png'))



    fig2=plt.figure()
    plt.errorbar(CR,fR_list,yerr=fR_err_list,fmt='ro')
    plt.xlabel('CR threshold')
    plt.ylabel('Rabi frequency [MHz]')
    plt.title(path)
    fig2.savefig(os.path.join(path,'CRvsRabiF.png'))

    average_f = np.mean(fR_list[7:15])
    average_tau = np.mean(tau_list[7:15])

    ave_fR.append(average_f)
    ave_tau.append(average_tau)

    all_fR=np.vstack((all_fR,fR_list[0:14]))
    all_tau=np.vstack((all_tau,tau_list[0:14]))
    all_fR_err=np.vstack((all_fR_err,fR_err_list[0:14]))
    all_tau_err=np.vstack((all_tau_err,tau_err_list[0:14]))

fig3 = plt.figure()
plt.plot(ave_fR,ave_tau,'bo')
plt.xlabel('average Rabi frequency (saturated part)')
plt.ylabel('average tau')

CRs = arange(14)

fig4 = plt.figure()

for i in CRs:
    plt.errorbar(all_fR[:,i],all_tau[:,i],
            xerr=all_fR_err[:,i],yerr=all_tau_err[:,i], fmt = 'o',
            label = 'CR threshold: %i' % CR[i])

plt.xlabel('fitted Rabi frequency [MHz]')
plt.ylabel('fitted Tau [ns]')
plt.legend()

np.savez(os.path.join(output_folder, 'LT1_fit_without_slope_range_34_61.npz'),
        powers_nW = [20,60,250], 
        rabi_f = all_fR, rabi_f_err = all_fR_err,
        damping = all_tau, damping_err = all_tau_err,
        CR = CR[0:14])
