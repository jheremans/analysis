from analysis.lib.fitting import fit, common
import os



LT2_datafolder = r'D:\measuring\data\LDE\analysis_data\homogeneous_scans\20120417'

LT2_homo = np.loadtxt(r'D:\measuring\data\LDE\analysis_data\homogeneous_scans\20120417\095929_Laserscan_sil9_without_gate\095929_Laserscan_sil9_without_gate_1.dat', skiprows = 15)

LT1_homo = np.loadtxt(r'D:\measuring\data\LDE\analysis_data\homogeneous_scans\20120518\115154_Laserscan_sil2_LT1_MW_SM_no_green\115154_Laserscan_sil2_LT1_MW_SM_no_green_1.dat',skiprows=15)


common_freq_axis = np.linspace(40,90,5001)

def find_nearest(array,value):
    idx=(abs(array-value)).argmin()
    return idx

def change_freq_axis(common_freq_axis, freq_axis, counts):
    new_counts=np.ones(len(common_freq_axis))*-1

    for i in np.arange(len(counts)):
        idx=find_nearest(common_freq_axis,freq_axis[i])
        try:
            if new_counts[idx]<0:
                new_counts[idx]=counts[i]
            else:
                new_counts[idx]+=counts[i]
        except:
            print 'something went wrong with index ', idx
    return new_counts

black_list = ['100621_Laserscan_sil9_300_1.dat','102627_Laserscan_sil9_305_1.dat',
        '105557_Laserscan_sil9_310_1.dat','111757_Laserscan_sil9_315_1.dat',
        '113108_Laserscan_sil9_320_1.dat','114625_Laserscan_sil9_325_1.dat',
        '121839_Laserscan_sil9_330_1.dat','124115_Laserscan_sil9_335_1.dat',
        '125107_Laserscan_sil9_340_1.dat']

d = np.loadtxt(r'D:\measuring\data\LDE\analysis_data\inhomogeneous_scans\20120424\180028_Laserscan_sil2_LT1_green\180028_Laserscan_sil2_LT1_green_1.dat',skiprows = 15)


LT2_folders = [i for i in os.listdir(LT2_datafolder) if 'Laserscan_sil9' in i and 'ionized' not in i]
file_names = []


all_counts = np.zeros(len(common_freq_axis))

fig1 = plt.figure()

for f in LT2_folders:
    file_name = [i for i in os.listdir(os.path.join(LT2_datafolder,f)) if i[-5:]=='1.dat' and i not in black_list]
    file_names.append(file_name)
    d = np.loadtxt(os.path.join(LT2_datafolder,f,file_name[0]),skiprows=15)
    freq = d[:,1]
    counts = d[:,2]
    
    plt.plot(freq,counts)

    all_counts += change_freq_axis(common_freq_axis,freq,counts)


##########LT2 homo
fit_range_LT2_homo = (56.2,56.8)  #GHz
center_LT2_homo = 56.5

freq_LT2_homo = LT2_homo[:,1]
counts_LT2_homo = LT2_homo[:,2]

idx_min_LT2_homo = find_nearest(freq_LT2_homo, fit_range[0])
idx_max_LT2_homo = find_nearest(freq_LT2_homo, fit_range[1])

freq_fit_LT2_homo = freq_LT2_homo[idx_min_LT2_homo:idx_max_LT2_homo]
counts_fit_LT2_homo = counts_LT2_homo[idx_min_LT2_homo:idx_max_LT2_homo]

fit_result_LT2_homo = fit.fit1d(freq_fit_LT2_homo,counts_fit_LT2_homo,
        common.fit_lorentz,0.,50.,center_LT2_homo,0.1,
        do_plot = True, newfig = True, ret = True)

LT2_homo_amp = fit_result_LT2_homo[0]['params_dict']['A']*2/(fit_result_LT2_homo[0]['params_dict']['gamma']*np.pi)
LT2_homo_center = fit_result_LT2_homo[0]['params_dict']['x0']


###########LT2 inhomo
freq_LT2_inhomo = common_freq_axis
counts_LT2_inhomo = all_counts

idx_min_LT2_inhomo = find_nearest(freq_LT2_inhomo, fit_range[0])
idx_max_LT2_inhomo = find_nearest(freq_LT2_inhomo, fit_range[1])

freq_fit_LT2_inhomo = freq_LT2_inhomo[idx_min_LT2_inhomo:idx_max_LT2_inhomo]
counts_fit_LT2_inhomo = counts_LT2_inhomo[idx_min_LT2_inhomo:idx_max_LT2_inhomo]


fit_result_LT2_inhomo = fit.fit1d(freq_fit_LT2_inhomo,counts_fit_LT2_inhomo,
        common.fit_gauss,0.,50.,56.45,0.1,
        do_plot = True, newfig = True, ret = True)

LT2_inhomo_center = fit_result_LT2_inhomo[0]['params_dict']['x0']


#########LT2 combined plot
x_offset = LT2_inhomo_center - LT2_homo_center

fig2 = plt.figure()
plt.plot(freq_fit_LT2_homo+x_offset, counts_fit_LT2_homo/LT2_homo_amp,'bo')
plt.plot(freq_fit_LT2_homo+x_offset, fit_result_LT2_homo[0]['fitdata']/LT2_homo_amp,'r-')

plt.plot(freq_fit_LT2_inhomo, counts_fit_LT2_inhomo/fit_result_LT2_inhomo[0]['params_dict']['A'],'ko')
plt.plot(freq_fit_LT2_inhomo, fit_result_LT2_inhomo[0]['fitdata']/fit_result_LT2_inhomo[0]['params_dict']['A'],'r-')
