
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc

datafolders=['no_Ey','1u_Ey','2u_Ey','3u_Ey']
RO_time=[0,1,2,3]
for i in datafolders:
    result= sc.plot_rabi(sc.get_latest_data(i))
    amp.append(2*result[0]['params'][1])
    phase.append(result[0]['params'][3])
    

plt.figure(1)
plt.plot(RO_time,amp,'bo')
plt.xlabel ('RO time [us]', fontsize = 16)
plt.ylabel ('Contrast', fontsize = 16)   
plt.ylim ([0, 1])
plt.show()

plt.figure(2)
plt.plot(RO_time,phase,'bo')
plt.xlabel ('RO time [us]', fontsize = 16)
plt.ylabel ('Phase [degree]', fontsize = 16)   
plt.ylim ([0, 1])
plt.show()
