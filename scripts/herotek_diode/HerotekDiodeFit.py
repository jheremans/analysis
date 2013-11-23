import os
import numpy
import matplotlib
import matplotlib.pyplot as plt

basefile_dir = 'C:\Users\Felix\Dropbox\Awschalom Lab Share\Calibration Files\HerotekCalibrations'
os.chdir(basefile_dir)
file_name = '0 _100MHz_Calibration.dat'
file_name = os.path.join(basefile_dir,file_name)
fp = open(file_name,'r')

data = fp.readlines()

SigGenPower = numpy.zeros(len(data))
HerotekVoltage = numpy.zeros(len(data))

ii = 0
for line in data:
    fields = line.split('\t')
    SigGenPower[ii] = fields[0]
    print(fields[0]) # prints the first fields value

    HerotekVoltage[ii] = fields[1]
    print(fields[1]) # prints the second fields value

    ii = ii + 1
fp.close()


# Plot one waveform

Xdata = SigGenPower
Ydata = HerotekVoltage

plt.plot(Xdata,Ydata,linestyle = 'None',marker='o')
plt.title('Herotek Calibration')
plt.ylabel('Herotek Voltage (V)')
plt.xlabel('Signal Generator Power (dBm)')
plt.grid(True)
plt.show()