import h5py
import os
import linecache
import numpy as np
from Tkinter import Tk

#A2
#X:\data\20130313\134621_AdwinSSRO_SIL2_SSROCalibration


#FB
#X:\data\20130313\141122_AdwinSSRO_SIL2_SSROCalibration
#X:\data\20130313\141931_AdwinSSRO_SIL2_SSROCalibration_forbidden2

def process_data(folder):
    
    for (path, dirs, files) in os.walk(folder):
        for file in files:
            if (str(file).find('.hdf5')!=-1) and (str(file).find('AdwinSSRO')!=-1):
                hdffn=os.path.join(path,file)
                timestamp=str(file)[0:6]
            if str(file)=='ssro.py':
                pyfn=os.path.join(path,file)
    
    #-------------sp----------
    
    a_sp_power=None
    sp_powerlns=find_lines(pyfn,"m.params['A_SP_amplitude']",164,20)
    for ln in sp_powerlns:
        sp_power=float(ln[ln.find('=')+1:ln.find('#')])
        if sp_power!=0.0:
            a_sp_power=sp_power
        
    f = h5py.File(hdffn, 'r')
    cal=f[f.keys()[0]]
    ms0=cal['ms0']
    
    a_sp_duration=ms0.attrs['SP_duration']
    sp_hist=ms0['SP_hist'].value
    sp_max_sat=float(max(sp_hist))/np.average(sp_hist[-4:])
    
    #-----------------cr-----------------
    #-----------------cr-----------------
    
    a_cr_power=cal.attrs['A_CR_amplitude']
    cr_threshold=cal.attrs['CR_preselect']
           
    cr_after=ms0['CR_after'].value
    p0_cr=float(len(np.where(cr_after==0)[0]))/len(cr_after)
    
    #----------------ro-------------------
    
    ro_counts = ms0['RO_data'].value.reshape(
                ms0.attrs['SSRO_repetitions'], ms0.attrs['SSRO_duration'])
    d = np.sum(ro_counts[:,:], axis=1)
    p0_ro = len(np.where(d==0)[0])/float(ms0.attrs['SSRO_repetitions'])
    
        
    ret= [timestamp,1-p0_ro,a_sp_power,a_sp_duration,sp_max_sat,a_cr_power,cr_threshold,p0_cr]
    copy_array_to_clipboard(ret)
    return ret
    
def find_line(fn,str,line,searchrange):
    for i in range(searchrange):
        for j in [-1,1]:
            ln=linecache.getline(fn,line+int(i*j))
            if ln.find(str) != -1:
                return ln
    return None
    
def find_lines(fn,str,line,searchrange):
    found=[]
    for i in range(searchrange):
        for j in [-1,1]:
            ln=linecache.getline(fn,line+int(i*j))
            if ln.find(str) != -1:
                found.append(ln)
    return found
    
def copy_array_to_clipboard(arr):
    to_str=''
    for a in arr:
        to_str=to_str+str(a)+'\t'
    r = Tk()
    r.withdraw()
    r.clipboard_clear()
    r.clipboard_append(to_str)
    r.destroy()