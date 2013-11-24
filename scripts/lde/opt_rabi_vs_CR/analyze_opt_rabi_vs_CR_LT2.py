import cPickle as pickle
from analysis.lib.lde import tail_cts_per_shot_v4


data_folder = r'D:\measuring\data\LDE\analysis_data\opt_rabi_vs_CR\2012-10-18_CRvsRabi'
output_folder = r'D:\measuring\data\LDE\analysis_output\20121104-opt-rabi-vs-CR'

##high rabi f measurements with long pulse:
files_hp_long = [
#            r'1525H_0001.txt',
            r'1523H_0101.txt',
            r'1521H_0202.txt',
            r'1520H_0403.txt',
            r'1518H_0604.txt',
            r'1516H_0806.txt',
            r'1514H_1008.txt',
            r'1512H_1209.txt',
            r'1511H_1612.txt',
            r'1509H_1814.txt',
            r'1534H_1814.txt',
            r'1537H_2217.txt',
            r'1541H_2720.txt',
            r'1544H_3023.txt',
            #r'1526H_0001.txt',
            #r'1426H_0001.txt',
            ]


##high rabi f measurements with short pulse:
files_hp_short = [
            r'1424H_0101.txt',
            r'1421H_0202.txt',
            r'1419H_0806.txt',
            r'1417H_1814.txt',
            r'1413H_1411.txt'       
         ]

CR_hp_long = [1,2,4,6,8,10,12,16,18,18,22,27,30]
tau_long_hp=[]
tau_long_hp_err=[]
fR_long_hp=[]
fR_long_hp_err=[]

##low rabi f measurements:
files_lp = [
        #r'1615H_0001_lp.txt',
        r'1614H_0101_lp.txt',
        r'1612H_0202_lp.txt',
        r'1611H_0403_lp.txt',
        r'1609H_0806_lp.txt',
        r'1608H_1814_lp.txt',
        #r'1617H_1814_lp.txt',
        r'1617H_3023_lp.txt'
        ]
CR_lp = [1,2,4,8,18,30]
tau_lp=[]
tau_lp_err=[]
fR_lp=[]
fR_lp_err=[]


for i in files_hp_long:
    print i
    file_name = os.path.join(data_folder,i)
    fR, fR_err, tau, tau_err = tail_cts_per_shot_v4.optical_rabi_resonant_fit_hannes(file_name,
            5,33, use_channel=1, binsize=0.128, rebins = 4)
    tau_long_hp.append(tau)
    tau_long_hp_err.append(tau_err)
    fR_long_hp.append(fR)
    fR_long_hp_err.append(fR_err)

#plt.close('all')

fig=plt.figure()
plt.errorbar(CR_hp_long,tau_long_hp,yerr=tau_long_hp_err,fmt='ro')
plt.xlabel('CR threshold')
plt.ylabel('Decay time [ns]')
plt.title('decay time vs CR threshold high power (ca. 200 MHz rabi f)')

for f in files_lp:
    print f
    file_name = os.path.join(data_folder,f)
    fR, fR_err, tau, tau_err = tail_cts_per_shot_v4.optical_rabi_resonant_fit_hannes(file_name,
            8,35, use_channel=1, binsize=0.128, rebins = 4)
    tau_lp.append(tau)
    tau_lp_err.append(tau_err)
    fR_lp.append(fR)
    fR_lp_err.append(fR_err)

fig=plt.figure()
plt.errorbar(CR_lp,tau_lp,yerr=tau_lp_err,fmt='ro')
plt.xlabel('CR threshold')
plt.ylabel('Decay time [ns]')
plt.title('decay time vs CR threshold low power (ca. 75 MHz rabi f)')

np.savez(os.path.join(output_folder, 'LT2_fit_without_slope_8_35_lp_5_33_hp.npz'),
        powers_nW = [100,880], 
        rabi_f = [fR_lp,fR_long_hp], rabi_f_err = [fR_lp_err,fR_long_hp_err],
        damping = [tau_lp,tau_long_hp], damping_err = [tau_lp_err,tau_long_hp_err],
        CR = [CR_lp, CR_hp_long])



