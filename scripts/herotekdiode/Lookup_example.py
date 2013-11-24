import os,sys,time

analysis_lib = 'Z:\\'
os.chdir(analysis_lib)
sys.path.append(analysis_lib)

from analysis.scripts.herotekdiode import HerotekDiodeLookup as herotek

# Converting for V to P
freq = 4 # GHz
volt = -2.5e-3 # mV
(power,nfreq) = herotek.convert_VtoP(freq,volt)

print ('VtoP converter: at %s GHz, %s V => %s dBm' % (nfreq,volt,power))


# Converting for P to V
freq = 4 # GHz
power = -15 # dBm
(volt,nfreq) = herotek.convert_PtoV(freq,power)

print ('PtoV convert: at %s GHz, %s dBm => %s mV' % (nfreq,power,volt*1e3))



