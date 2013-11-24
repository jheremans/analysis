from analysis.lib.m2.ssro import ssro
import numpy as np


folder=r'D:\measuring\data\20130619\141204_AdwinSSRO_SSROCalibration_sil9_A2'#tb.latest_data('AdwinSSRO')
a = ssro.SSROAnalysis(folder)
  
n='ms0'    
a.get_run(n)
a.charge_hist(a.cr_counts, name=n)
a.spinpumping(a.sp_time, a.sp_counts, a.reps, a.binsize, name=n)

v_reps=len(np.where(a.cr_counts>0)[0])
v_ro_counts=a.ro_counts[a.cr_counts>0]
a.cpsh_hist(v_ro_counts, v_reps, name=n)
a.readout_relaxation(a.ro_time, v_ro_counts, v_reps, a.binsize, name=n)
a.fidelity(v_ro_counts, v_reps, a.binsize, 0, name=n)

print float(len(np.where(np.sum(v_ro_counts,axis=1)==0)[0]))/v_reps