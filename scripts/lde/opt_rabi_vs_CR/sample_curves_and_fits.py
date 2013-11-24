import cPickle as pickle


LT2_data_folder = r'D:\measuring\data\LDE\analysis_data\opt_rabi_vs_CR\2012-10-18_CRvsRabi'

LT1_data_folder = r'D:\measuring\data\LDE\analysis_data\opt_rabi_vs_CR\2012-10-24-CRvsRabi_LT1\183600_rabi_vs_cr\60nW'

output_folder = r'D:\measuring\data\LDE\analysis_output\20121104-opt-rabi-vs-CR'


def find_nearest(array,value):
    idx=(abs(array-value)).argmin()
    return idx

###############
#example curve with fit
df_LT1=os.path.join(LT1_data_folder,r'40_60nW_2006.txt')
df_LT2=os.path.join(LT2_data_folder,r'1608H_1814_lp.txt')

data_LT1 = np.loadtxt(df_LT1, skiprows=10)
data_LT2 = np.loadtxt(df_LT2, skiprows=10)

counts_LT1 = data_LT1[:,0]
counts_LT2 = data_LT2[:,1]


###fitting
fR_init_guess = 0.100 #[GHz]
tau_init_guess = 15 #[ns]
ampl_init_guess = 1000.
offset_init_guess = 0.0
slope_init_guess = 0.0#1.2


#########LT1
min_fit_bnd = 32
max_fit_bnd = 76

x = np.arange(0,len(counts_LT1))
time = x*0.128

fit_range = range(find_nearest(time,min_fit_bnd),
    find_nearest(time,max_fit_bnd))
time_fit = time[fit_range]
counts_fit = counts_LT1[fit_range]

figure = plt.figure(figsize=(8., 6.))
plt.hold(True)

offset_init_guess = 30.6

fit_result = fit.fit1d(time_fit, counts_fit, rabi.fit_rabi_damped_exp_with_offset_on_linslope, 
        fR_init_guess,ampl_init_guess, offset_init_guess, 
        tau_init_guess, time_fit[0], slope_init_guess, fixed=[],
        do_plot = True, newfig = False, ret = True,
        plot_fitonly = True, color = 'r', linewidth = 2.0)

plt.plot(time, counts_LT1, '.', color = 'k')

f = open(os.path.join(output_folder, r'LT1_sample_curve.pkl'),'wb')
pickle.dump([{'counts':counts_LT1,'time_fit':time_fit,'time_axis':time,'fit_result':fit_result}],f)
f.close()


######LT2
min_fit_bnd = 6
max_fit_bnd = 49

x = np.arange(0,len(counts_LT2))
time = x*0.128

fit_range = range(find_nearest(time,min_fit_bnd),
    find_nearest(time,max_fit_bnd))
time_fit = time[fit_range]
counts_fit = counts_LT2[fit_range]

figure = plt.figure(figsize=(8., 6.))
plt.hold(True)

offset_init_guess = 4.8

fit_result = fit.fit1d(time_fit, counts_fit, rabi.fit_rabi_damped_exp_with_offset_on_linslope, 
        fR_init_guess,ampl_init_guess, offset_init_guess, 
        tau_init_guess, time_fit[0], slope_init_guess, fixed=[],
        do_plot = True, newfig = False, ret = True,
        plot_fitonly = True, color = 'r', linewidth = 2.0)

plt.plot(time, counts_LT2, '.', color = 'k')

f = open(os.path.join(output_folder, r'LT2_sample_curve.pkl'),'wb')
pickle.dump([{'counts':counts_LT2,'time_fit':time_fit,'time_axis':time,'fit_result':fit_result}],f)
f.close()
