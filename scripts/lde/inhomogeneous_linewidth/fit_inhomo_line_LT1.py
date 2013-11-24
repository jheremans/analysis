from analysis.lib.fitting import fit, common
import cPickle as pickle
import os
output_folder = r'D:\measuring\data\LDE\analysis_output\20121113-linewidths'



LT1_datafolder = r'D:\measuring\data\LDE\analysis_data\homogeneous_scans\LT1\20120412'

LT1_homo = np.loadtxt(r'D:\measuring\data\LDE\analysis_data\homogeneous_scans\LT1\20120412\144729_Laserscan_10deg\scan3.txt')

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


LT1_folders = [i for i in os.listdir(LT1_datafolder) if 'Laserscan' in i]
file_names = []

blacklist = ['152400_Laserscan_35degscan3.txt',
                '152400_Laserscan_35degscan1.txt',
                '150924_Laserscan_25degscan2.txt',
                '150924_Laserscan_25degscan3.txt',
                '150203_Laserscan_20degscan3.txt',
                '145446_Laserscan_15degscan4.txt',
                '143011_Laserscan_0degscan2.txt',
                #'143011_Laserscan_0degscan0.txt',
                '151640_Laserscan_30degscan2.txt'
                ]

all_counts = np.zeros(len(common_freq_axis))

center_freq = 60.91
fit_range_LT1_homo = (60.2,61.2)  #GHz



shifted_counts = np.zeros(len(common_freq_axis))
failed_fit = []


fig1 = plt.figure()
ax = plt.subplot(111)

for f in LT1_folders:
    file_name = [i for i in os.listdir(os.path.join(LT1_datafolder,f)) \
            if i[-4:]=='.txt' \
            and 'raw' not in i and '2D' not in i and 'scan' in i and 'parameters' not in i]
    for s in file_name:
        d = np.loadtxt(os.path.join(LT1_datafolder,f,s))
        freq = d[:,0]/1000.
        counts = d[:,1]
        idx_min = find_nearest(freq, fit_range_LT1_homo[0])
        idx_max = find_nearest(freq, fit_range_LT1_homo[1])
        center_guess = freq[find_nearest(counts,max(counts))]

        if np.max(counts[idx_min:idx_max])>5000 and f+s not in blacklist:
            ax.plot(freq,counts)                                          
            freq_fit = freq[idx_min:idx_max]
            counts_fit = counts[idx_min:idx_max]

            fit_result = fit.fit1d(freq_fit,counts_fit,
                    common.fit_2lorentz,50.,400.,center_freq,0.02,30.,center_freq-0.2,0.02,
                    do_plot = True, newfig = True, ret = True)
            if fit_result:
                offset = fit_result[0]['params_dict']['x01']-center_freq
                new_counts = change_freq_axis(common_freq_axis, freq-offset, counts)
                all_counts += change_freq_axis(common_freq_axis,freq,counts)
                shifted_counts+=new_counts

                idx_min = find_nearest(common_freq_axis, fit_range_LT1_homo[0])
                idx_max = find_nearest(common_freq_axis, fit_range_LT1_homo[1])
                print f,s
                print offset, fit_result[0]['params_dict']['x01']
                print common_freq_axis[find_nearest(new_counts,max(new_counts[idx_min:idx_max]))]
                file_names.append(s)
                title=os.path.join(f,s)
                plt.title(title)
            else:
                failed_fit.append(os.path.join(f,s))



##########LT1 homo


freq_LT1_homo = LT1_homo[:,0]/1000# common_freq_axis
counts_LT1_homo = LT1_homo[:,1]#shifted_counts

idx_min_LT1_homo = find_nearest(freq_LT1_homo, fit_range_LT1_homo[0])
idx_max_LT1_homo = find_nearest(freq_LT1_homo, fit_range_LT1_homo[1])

freq_fit_LT1_homo = freq_LT1_homo[idx_min_LT1_homo:idx_max_LT1_homo]
counts_fit_LT1_homo = counts_LT1_homo[idx_min_LT1_homo:idx_max_LT1_homo]

fit_result_LT1_homo = fit.fit1d(freq_fit_LT1_homo,counts_fit_LT1_homo,
        common.fit_2lorentz,0.,5000.,center_freq,0.03,300.,60.7,0.03,fit_curve_points = 2001,
        do_plot = True, newfig = True, ret = True)

LT1_homo_amp = fit_result_LT1_homo[0]['params_dict']['A1']*2/(fit_result_LT1_homo[0]['params_dict']['gamma1']*np.pi)+fit_result_LT1_homo[0]['params_dict']['a1']
LT1_homo_center = fit_result_LT1_homo[0]['params_dict']['x01']


###########LT1 inhomo
freq_LT1_inhomo = common_freq_axis
counts_LT1_inhomo = all_counts

idx_min_LT1_inhomo = find_nearest(freq_LT1_inhomo, fit_range_LT1_homo[0])
idx_max_LT1_inhomo = find_nearest(freq_LT1_inhomo, fit_range_LT1_homo[1])

freq_fit_LT1_inhomo = freq_LT1_inhomo[idx_min_LT1_inhomo:idx_max_LT1_inhomo]
counts_fit_LT1_inhomo = counts_LT1_inhomo[idx_min_LT1_inhomo:idx_max_LT1_inhomo]


fit_result_LT1_inhomo = fit.fit1d(freq_fit_LT1_inhomo,counts_fit_LT1_inhomo,
        common.fit_2gauss,0.,20000.,60.9,0.1,2000,60.7,0.1,fit_curve_points = 2001,
        do_plot = True, newfig = True, ret = True)

LT1_inhomo_center = fit_result_LT1_inhomo[0]['params_dict']['x01']
LT1_inhomo_amp = fit_result_LT1_inhomo[0]['params_dict']['A1']+fit_result_LT1_inhomo[0]['params_dict']['a1']

#########LT1 combined plot
x_offset = LT1_inhomo_center - LT1_homo_center

fig2 = plt.figure()
plt.plot(freq_fit_LT1_homo+x_offset, counts_fit_LT1_homo/LT1_homo_amp,'bo')
plt.plot(fit_result_LT1_homo[0]['fitx']+x_offset, fit_result_LT1_homo[0]['fity']/LT1_homo_amp,'r-')

plt.plot(freq_fit_LT1_inhomo, counts_fit_LT1_inhomo/LT1_inhomo_amp,'ko')
plt.plot(fit_result_LT1_inhomo[0]['fitx'], fit_result_LT1_inhomo[0]['fity']/LT1_inhomo_amp,'r-')

f = open(os.path.join(output_folder, r'LT1_linewidths.pkl'),'wb')
pickle.dump([{'counts_homo' : counts_fit_LT1_homo,
                'homo_axis' : freq_fit_LT1_homo+x_offset,
                'norm_homo' : LT1_homo_amp,
                'fit_result_homo' : fit_result_LT1_homo[0],
                'homo_offset' : x_offset,
                'counts_inhomo': counts_fit_LT1_inhomo,
                'inhomo_axis' : freq_fit_LT1_inhomo,
                'norm_inhomo' : LT1_inhomo_amp,
                'fit_result_inhomo':fit_result_LT1_inhomo[0],
                'noof_scans' : len(file_names)}],f)
f.close()
