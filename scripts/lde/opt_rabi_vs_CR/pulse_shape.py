from analysis.lib.lde import tail_cts_per_shot_v4

def find_nearest(array,value):
    idx=(abs(array-value)).argmin()
    return idx


d=numpy.loadtxt(r'D:\measuring\data\LDE\analysis_data\opt_rabi_vs_CR\2012-10-24-CRvsRabi_LT1\183600_rabi_vs_cr\pulseshape_50ns_250nW.txt', skiprows = 10)

counts0=d[:,0]
counts1=d[:,1]


counts_rebin0=tail_cts_per_shot_v4.rebin(counts0,4)
counts_rebin1=tail_cts_per_shot_v4.rebin(counts1,4)

time_ax = arange(len(counts_rebin0))*4*0.128

fit_min = 58.0 #ns
fit_max = 100.0 #ns

idx_min = find_nearest(time_ax,fit_min)
idx_max = find_nearest(time_ax,fit_max)

fig=plt.figure()

plt.plot(time_ax,counts_rebin1,'bo')

plt.xlabel('time (ns)')
plt.ylabel('counts')
plt.title('pulse shape')

fit_result = fit.fit1d(time_ax[idx_min:idx_max],
        counts_rebin1[idx_min:idx_max],
        common.fit_line,2000,0,do_plot = True, 
        newfig = False, plot_fitonly = True, color = 'r', linewidth = 2.0)
