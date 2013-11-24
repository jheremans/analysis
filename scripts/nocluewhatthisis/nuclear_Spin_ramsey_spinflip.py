
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot
from analysis.lib.spin import spin_control as sc

#datafolders=['1154','1252','1258','1303','1310','1347','1351','1316','1326', '1453']
#RO_time=[0,1,2,3,4,5,6,7,9, 11]
datafolders=['144831','145714','150926','151701','152503','153402','154024']
RO_power=[0,12,5,20,30,25,8]
amp=[]
phase=[]
for i in datafolders:
    result= sc.plot_rabi(sc.get_latest_data(i))
    result= sc.plot_rabi(sc.get_latest_data(i))
    #qt.sleep(4)
    amp.append(abs(2*result[0]['params'][1]))
    phase.append(result[0]['params'][3])
'''
datafoldersCond=['1502','4us','6us','8us', '10us', '12us']
RO_timeCond=[2,4,6,8, 10, 12]
ampCond=[]
phaseCond=[]
for i in datafoldersCond:
    result= sc.plot_rabi(sc.get_latest_data(i,datapath=r'D:\measuring\data\20121219'))
    #qt.sleep(4)
    ampCond.append(abs(2*result[0]['params'][1]))
    phaseCond.append(result[0]['params'][3])
'''
result=fit.fit1d(np.array(RO_power), np.array(amp/amp[0]), common.fit_exp_decay_with_offset, 
            0, 1, 5,
            do_plot = True, do_print = True, newfig = False,ret=True)

#resultcond=fit.fit1d(np.array(RO_timeCond), np.array(ampCond/amp[0]), common.fit_exp_decay_with_offset, 
#            0, 1, 5,
#            do_plot = True, do_print = True, newfig = False,ret=True)
print 'normal RO'
print result[0]['params']
#print 'cond RO'
#print resultcond[0]['params']

plt.figure(9)
plt.plot(RO_power,amp/amp[0],'bo')
plt.plot(result[0]['fitx'],result[0]['fity'],'b-')
#plt.plot(resultcond[0]['fitx'],resultcond[0]['fity'],'r-')
#plt.plot (RO_timeCond, ampCond/amp[0], 'r<')
plt.xlabel ('RO time [us]', fontsize = 16)
plt.ylabel ('Contrast', fontsize = 16)   
plt.ylim ([0, 1])
plt.show()

plt.figure(11)
plt.plot(RO_power,mod(phase,180),'bo')
#plt.plot (RO_timeCond, mod(phaseCond, 360), 'r<')
plt.xlabel ('RO Power [nW]', fontsize = 16)
plt.ylabel ('Phase [degree]', fontsize = 16)   
plt.ylim ([0, 360])
plt.show()


