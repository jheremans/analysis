import os, sys, time
import numpy as np
from matplotlib import pyplot as plt

LT1_file=r'K:\ns\qt\Diamond\Data\LDE\analysis_data\laserscans\2012-10-24-132200_laserscan_SIL2_paper_scans_,LT15nW_mw_True_0nW_green_0gate_phase0,-1\132200_laserscan_SIL2_paper_scans_,LT15nW_mw_True_0nW_green_0gate_phase0,-1_1.dat'

LT2_file=r'K:\ns\qt\Diamond\Data\LDE\analysis_data\laserscans\2012-10-18-192833_LT25nW_mw_True_0nW_green_0gate_+0\192833_LT25nW_mw_True_0nW_green_0gate_+0_1.dat'

lt1_scan=np.loadtxt(LT1_file)
lt2_scan=np.loadtxt(LT2_file)

fig = figure(figsize=(3,1))

ax = plt.subplot(111)
#ax.locator_params(nbins=6)
plt.plot(lt1_scan[:,1],lt1_scan[:,2]+50,linewidth=2, label='NV a')
plt.plot(lt2_scan[:,1],lt2_scan[:,2],linewidth=2, label='NV b')
plt.legend()
plt.xlim(50,72)
plt.ylim(-5,400)

fig.savefig(r'D:\experiments\LDE\ldeplots\laserscan2.pdf')
