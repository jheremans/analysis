import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt
from analysis import config

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'
config.datadir = r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde'

srcfolder = config.datadir
savefolder = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-ssro-stats')

times = np.array([])
lt1durations = np.array([])
lt2durations = np.array([])

for (path,dirs,files) in os.walk(srcfolder):
    base,dn = os.path.split(path)
    if 'adwin_lt2-000' in dn:   
        f = np.load(os.path.join(path, 'lde_params.npz'))
        duration_lt2 = f['SSRO_duration']
        duration_lt1 = f['SSRO_duration_lt1']
        f.close()

        _b,timefolder = os.path.split(base)
        _b,datefolder = os.path.split(_b)
        timestr = timefolder[:6]
        datestr = datefolder
        datetimestr = datestr+timestr
        t = int(time.mktime(time.strptime(datetimestr, '%Y%m%d%H%M%S')))
        
        times = np.append(times, t)
        lt1durations = np.append(lt1durations, duration_lt1)
        lt2durations = np.append(lt2durations, duration_lt2)

if not os.path.exists(savefolder):
    os.makedirs(savefolder)

np.savez(os.path.join(savefolder, 'ssro_durations'), 
        times=times, lt1durations=lt1durations, lt2durations=lt2durations)







