import os, sys, time
import numpy as np

analysis_lib = 'Z:\\'
os.chdir(analysis_lib)
sys.path.append(analysis_lib)
##print os.listdir('.')

# define parameter columns
fcol = 0
acol = 1
tcol = 2
ocol = 3

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return (idx,array[idx])

def fitfunc_powertovoltage(p,a,tau,o):
    return o+a*np.exp((p)/(tau))

def fitfunc_voltagetopower(v,a,tau,o):
    return tau*(np.log((v-o)/(a)))

def load_calibrationfile():
    base_dir = 'Z:\\calibrations'
    file_name = '203918_herotek_calibration_fitparams.dat'
    file_name = os.path.join(base_dir,file_name)
    ##    fp = open(file_name,'r')
    ##    colheader = fp.readline()
    ##    fp.close()
    data = np.loadtxt(file_name,skiprows = 1)
    return data

def convert_VtoP(freq,volt):
    data = load_calibrationfile()
    (idx,nearest_freq) = find_nearest(data[:,fcol], freq)
    power = fitfunc_voltagetopower(volt,data[idx,acol],data[idx,tcol],data[idx,ocol])
    return (power,nearest_freq)

def convert_PtoV(freq,power):
    data = load_calibrationfile()
    (idx,nearest_freq) = find_nearest(data[:,fcol], freq)
    voltage = fitfunc_powertovoltage(power,data[idx,acol],data[idx,tcol],data[idx,ocol])
    return (voltage,nearest_freq)

