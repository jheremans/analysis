from analysis.lib.fitting import fit, common
import cPickle as pickle
output_folder = r'D:\measuring\data\LDE\analysis_output\20121113-linewidths'
import os



LT2_datafolder = r'D:\measuring\data\LDE\analysis_data\homogeneous_scans\20120417'

LT2_homo = np.loadtxt(r'D:\measuring\data\LDE\analysis_data\homogeneous_scans\20120417\124115_Laserscan_sil9_335_1\124115_Laserscan_sil9_335_1_1.dat', skiprows = 15)



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
failed_fit = []


all_counts = np.zeros(len(common_freq_axis))
shifted_counts = np.zeros(len(common_freq_axis))
center_freq = 56.4
fit_range_LT2_homo = (56.2,56.8)  #GHz


fig1 = plt.figure()
ax = plt.subplot(111)

for f in LT2_folders:
    file_name = [i for i in os.listdir(os.path.join(LT2_datafolder,f)) \
            if i[-5:]=='1.dat' and i not in black_list]
    
    d = np.loadtxt(os.path.join(LT2_datafolder,f,file_name[0]),skiprows=15)
    freq = d[:,1]
    counts = d[:,2]
    idx_min = find_nearest(freq, fit_range_LT2_homo[0])
    idx_max = find_nearest(freq, fit_range_LT2_homo[1])

    if np.max(counts[idx_min:idx_max])>150:
        ax.plot(freq,counts)
                              
        freq_fit = freq[idx_min:idx_max]
        counts_fit = counts[idx_min:idx_max]
        center_guess = freq_fit[find_nearest(counts_fit,max(counts_fit))]


        fit_result = fit.fit1d(freq_fit,counts_fit,
                common.fit_lorentz,50.,140.,center_guess,0.02,
                do_plot = True, newfig = True, ret = True)
        if fit_result:
            offset = fit_result[0]['params_dict']['x0']-center_freq
            print offset, fit_result[0]['params_dict']['x0']
            new_counts = change_freq_axis(common_freq_axis, freq-offset, counts)
            shifted_counts+=new_counts
            all_counts += change_freq_axis(common_freq_axis,freq,counts)   
            file_names.append(file_name)
            idx_min = find_nearest(common_freq_axis, fit_range_LT2_homo[0])
            idx_max = find_nearest(common_freq_axis, fit_range_LT2_homo[1])
            print common_freq_axis[find_nearest(new_counts,max(new_counts[idx_min:idx_max]))]

            print file_name
            plt.title(file_name)
        else:
            failed_fit.append(file_name)



##########LT2 homo
fit_range_LT2_homo = (56.2,56.8)  #GHz
center_LT2_homo = 56.5

freq_LT2_homo = LT2_homo[:,1]#common_freq_axis
counts_LT2_homo = LT2_homo[:,2]#shifted_counts

idx_min_LT2_homo = find_nearest(freq_LT2_homo, fit_range_LT2_homo[0])
idx_max_LT2_homo = find_nearest(freq_LT2_homo, fit_range_LT2_homo[1])

freq_fit_LT2_homo = freq_LT2_homo[idx_min_LT2_homo:idx_max_LT2_homo]
counts_fit_LT2_homo = counts_LT2_homo[idx_min_LT2_homo:idx_max_LT2_homo]

fit_result_LT2_homo = fit.fit1d(freq_fit_LT2_homo,counts_fit_LT2_homo,
        common.fit_lorentz,0.,5000.,center_LT2_homo,0.1,
        do_plot = True, newfig = True, ret = True)

LT2_homo_amp = fit_result_LT2_homo[0]['params_dict']['A']*2\
        /(fit_result_LT2_homo[0]['params_dict']['gamma']*np.pi)+fit_result_LT2_homo[0]['params_dict']['a']
LT2_homo_center = fit_result_LT2_homo[0]['params_dict']['x0']


###########LT2 inhomo
freq_LT2_inhomo = common_freq_axis
counts_LT2_inhomo = all_counts

idx_min_LT2_inhomo = find_nearest(freq_LT2_inhomo, fit_range_LT2_homo[0])
idx_max_LT2_inhomo = find_nearest(freq_LT2_inhomo, fit_range_LT2_homo[1])

freq_fit_LT2_inhomo = freq_LT2_inhomo[idx_min_LT2_inhomo:idx_max_LT2_inhomo]
counts_fit_LT2_inhomo = counts_LT2_inhomo[idx_min_LT2_inhomo:idx_max_LT2_inhomo]


fit_result_LT2_inhomo = fit.fit1d(freq_fit_LT2_inhomo,counts_fit_LT2_inhomo,
        common.fit_gauss,1000.,50000.,56.45,0.1,
        do_plot = True, newfig = True, ret = True)

LT2_inhomo_center = fit_result_LT2_inhomo[0]['params_dict']['x0']
LT2_inhomo_amp = fit_result_LT2_inhomo[0]['params_dict']['A']+fit_result_LT2_inhomo[0]['params_dict']['a']

#########LT2 combined plot
x_offset = LT2_inhomo_center - LT2_homo_center

fig2 = plt.figure()
plt.plot(freq_fit_LT2_homo+x_offset, counts_fit_LT2_homo/LT2_homo_amp,'bo')
plt.plot(fit_result_LT2_homo[0]['fitx']+x_offset, fit_result_LT2_homo[0]['fity']/LT2_homo_amp,'r-')

plt.plot(freq_fit_LT2_inhomo, counts_fit_LT2_inhomo/LT2_inhomo_amp,'ko')
plt.plot(fit_result_LT2_inhomo[0]['fitx'], fit_result_LT2_inhomo[0]['fity']/LT2_inhomo_amp,'r-')

f = open(os.path.join(output_folder, r'LT2_linewidths.pkl'),'wb')
pickle.dump([{'counts_homo' : counts_fit_LT2_homo,
                'homo_axis' : freq_fit_LT2_homo+x_offset,
                'norm_homo' : LT2_homo_amp,
                'homo_offset' : x_offset,
                'fit_result_homo' : fit_result_LT2_homo[0],
                'counts_inhomo': counts_fit_LT2_inhomo,
                'inhomo_axis' : freq_fit_LT2_inhomo,
                'norm_inhomo' : LT2_inhomo_amp,
                'fit_result_inhomo':fit_result_LT2_inhomo[0],
                'noof_scans' : len(file_names)}],f)
f.close()
