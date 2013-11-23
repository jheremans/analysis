
from numpy import *
import pylab as plt
import numpy as np
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
#from analysis.lib.spin import spin_control as sc

'''
#datafolders=['1154','1252','1258','1303','1310','1347','1351','1316','1326', '1453']
#RO_time=[0,1,2,3,4,5,6,7,9, 11]
SIL9=[61.8,66.6,r'D:\measuring\data\20130627\122204_LaserFrequencyScan_red_scan_coarse_gv_0.0\122204_LaserFrequencyScan_red_scan_coarse_gv_0.0.npz',9]
SIL10=[56.9,57.3,r'D:\measuring\data\20130701\140903_LaserFrequencyScan_SIL10_red_mw_gv_0.0\140903_LaserFrequencyScan_SIL10_red_mw_gv_0.0.npz',10]
SIL11=[53.62,55.61,r'D:\measuring\data\20130702\093927_LaserFrequencyScan_SIL11_green_gv_0.0\093927_LaserFrequencyScan_SIL11_green_gv_0.0.npz',11]
SIL18=[55.6,59.5,r'D:\measuring\data\20130702\102752_LaserFrequencyScan_SIL18_green_gv_0.0\102752_LaserFrequencyScan_SIL18_green_gv_0.0.npz',18]
SIL15=[59.6,63.2,r'D:\measuring\data\20130702\112428_LaserFrequencyScan_SIL15_Green_gv_0.0\112428_LaserFrequencyScan_SIL15_Green_gv_0.0.npz',15]
SIL12=[0,0,r'D:\measuring\data\20130702\175735_LaserFrequencyScan_SIL12_green_gv_0.0\175735_LaserFrequencyScan_SIL12_green_gv_0.0.npz',12]
SIL8=[55.5,60.4,r'D:\measuring\data\20130703\132950_LaserFrequencyScan_SIL8_red_gv_0.0\132950_LaserFrequencyScan_SIL8_red_gv_0.0.npz',8]
SIL7=[52.7,59.5,r'D:\measuring\data\20130702\122103_LaserFrequencyScan_SIL7_Green_gv_0.0\122103_LaserFrequencyScan_SIL7_Green_gv_0.0.npz',7]
SIL5=[40.8,50.3,r'D:\measuring\data\20130702\173219_LaserFrequencyScan_SIL5_green_gv_0.0\173219_LaserFrequencyScan_SIL5_green_gv_0.0.npz',5]
SIL1 = [71,78.4,r'D:\measuring\data\20130702\171520_LaserFrequencyScan_SIL1_red_gv_0.0\171520_LaserFrequencyScan_SIL1_red_gv_0.0.npz',1]
datafolders=[SIL1,SIL5,SIL7,SIL8,SIL9,SIL10,SIL11,SIL15,SIL18]

SIL1 = [32.6,49.3,r'D:\measuring\data\20130812\132500_LaserFrequencyScan_red_scan_coarse_Hans_SIL1_1nW_mw_gv_0.0\132500_LaserFrequencyScan_red_scan_coarse_Hans_SIL1_1nW_mw_gv_0.0.npz',1]
SIL3 = [39.4,56.7,r'D:\measuring\data\20130812\135005_LaserFrequencyScan_red_scan_coarse_Hans_SIL3_1nW_mw_gv_0.0\135005_LaserFrequencyScan_red_scan_coarse_Hans_SIL3_1nW_mw_gv_0.0.npz',3]
SIL4 = [61,63.17,r'D:\measuring\data\20130812\142112_LaserFrequencyScan_red_scan_coarse_Hans_SIL4_1nW_mw_gv_0.0\142112_LaserFrequencyScan_red_scan_coarse_Hans_SIL4_1nW_mw_gv_0.0.npz',4]
SIL5 = [35.2,50.5,r'D:\measuring\data\20130812\154026_LaserFrequencyScan_red_scan_coarse_Hans_SIL5_1nW_mw_gv_0.0\154026_LaserFrequencyScan_red_scan_coarse_Hans_SIL5_1nW_mw_gv_0.0.npz',5]
SIL7 = [68,74.5,r'D:\measuring\data\20130812\160400_LaserFrequencyScan_red_scan_coarse_Hans_SIL7_1nW_mw_gv_0.0\160400_LaserFrequencyScan_red_scan_coarse_Hans_SIL7_1nW_mw_gv_0.0.npz',7]
'''

SIL9 = [55.6,62.1,r'X:\data\20131015\121526_LaserFrequencyScan_red_scan_coarse_lt2_sil9_MWs_gv_0']
SIL10 = [55.6,62.1,r'X:\data\20131015\132133_LaserFrequencyScan_red_scan_coarse_lt2_sil10_Green_gv_0']
SIL18 = [50.12,55.0,r'X:\data\20131015\134109_LaserFrequencyScan_red_scan_coarse_lt2_sil18_Green_gv_0']
SIL5 = [39.6,48.97,r'X:\data\20131015\135226_LaserFrequencyScan_red_scan_coarse_lt2_sil5_Green_gv_0']



datafolders=[SIL9]#,SIL8,SIL11,SIL15,SIL18]
plt.figure(9)
j=0
n=25.5
for i in datafolders:
    
    result = np.load(i[2])
    r=result['data']
    plt.plot(r[:,1],r[:,2]/float(max(r[:,2]))+j)
    if (i[0]!=0):
        plt.text(i[0]-1,j+0.5,'Ey',fontsize=8)
    if (i[1]!=0):
        plt.text(i[1]+0.3,j+0.5,'Ex',fontsize=8)
    name='SIL'+str(i[3])
    plt.text(n,j+0.1,name,fontsize=10)
    j=j+1
    result.close()
plt.xlim([25,75])
plt.ylim([-0.1,j])
plt.xlabel ('Relative Frequency [GHz]')   
plt.ylabel ('Counts (Normalized)')   


