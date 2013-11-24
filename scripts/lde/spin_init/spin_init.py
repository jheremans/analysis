from analysis.lib.fitting import fit, common

ms1_file = np.load(r'Z:\data\20121012\063628_ADwin_SSRO_SIL2_LT1\ADwin_SSRO-001_SP_histogram.npz')

ms0_file = np.load(r'Z:\data\20121012\063628_ADwin_SSRO_SIL2_LT1\ADwin_SSRO-000_SP_histogram.npz')


ms1_time = ms1_file['time']/1000.
ms1_counts = ms1_file['counts']/20000.*1e6
ms1_file.close()

ms0_time = ms0_file['time']/1000.
ms0_counts = ms0_file['counts']/20000.*1e6
ms0_file.close()

fig1 = plt.figure()
plt.plot(ms1_time[:100],ms1_counts[:100],'ro')
fit_result = fit.fit1d(ms1_time[4:100],ms1_counts[4:100],
        common.fit_exp_decay_shifted_with_offset,
        1000,6000,2,4,fixed=[3],
        do_plot=True,newfig = False, plot_fitonly = True, 
        color = 'b', linewidth = 2.0, ret=True)
plt.xlabel('time [us]')
plt.ylabel('countrate [Hz]')
plt.title(r'20121012\063628_ADwin_SSRO_SIL2_LT1')

A=fit_result[0]['params_dict']['A']
B=fit_result[0]['params_dict']['a']
A_err=fit_result[0]['error_dict']['A']
B_err=fit_result[0]['error_dict']['a']
spin_init = B/(A+B)
spin_init_err = np.sqrt((-B/(A+B)**2*A_err)**2+(-B/(A+B)**2*B_err)**2)

plt.text(ms1_time.max()*0.2,ms1_counts.max()*0.8,
        'init error %.2f +/- %.3f %%' % (spin_init*100, spin_init_err*100),fontsize=16)

fig2 = plt.figure()
plt.plot(ms0_time[:100],ms0_counts[:100],'ro')
fit_result = fit.fit1d(ms0_time[3:100],ms0_counts[3:100],
        common.fit_exp_decay_shifted_with_offset,
        1000,6000,2,3,fixed=[3],
        do_plot=True,newfig = False, plot_fitonly = True, 
        color = 'b', linewidth = 2.0, ret=True)
plt.xlabel('time [us]')
plt.ylabel('countrate [Hz]')
plt.title(r'20121012\063628_ADwin_SSRO_SIL2_LT1')

A=fit_result[0]['params_dict']['A']
B=fit_result[0]['params_dict']['a']
A_err=fit_result[0]['error_dict']['A']
B_err=fit_result[0]['error_dict']['a']
spin_init = B/(A+B)
spin_init_err = np.sqrt((-B/(A+B)**2*A_err)**2+(-B/(A+B)**2*B_err)**2)

plt.text(ms0_time.max()*0.2,ms0_counts.max()*0.8,
        'init error %.2f +/- %.3f %%' % (spin_init*100, spin_init_err*100),fontsize=16)
