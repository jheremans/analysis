import os,sys
analysis_lib = 'C:\Users\Felix\Documents\GitHub'
os.chdir(analysis_lib)
sys.path.append(analysis_lib)
print os.listdir('.')

from analysis.lib import fitting
from analysis.lib.fitting import fit
from analysis.lib.tools import plot

from matplotlib import pyplot as plt
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


# Colors
SteelBlue = '#4682b4'
Maroon = '#800000'

# Plot one waveform
Xdata = SigGenPower
Ydata = HerotekVoltage

x = Xdata
y = Ydata

# fit data
a = fit.Parameter(-0.2, 'a')
##o = fit.Parameter(0., 'o')
c = fit.Parameter(7, 'c')
##x0 = fit.Parameter(16,'x0')
fitfunc_str = 'y(x) = a*exp(x/c)'

def fitfunc(x):
    return a() * numpy.exp((x) / (c()))


fig = plt.figure()
ax = fig.add_subplot(111)

fit_result = fit.fit1d(x,y, None, p0=[c,a], fixed=[], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
ax = plot.plot_fit1d(fit_result, numpy.linspace(x[0],x[-1],201), ax=ax, ret = 'ax',
        plot_data=True)

def fixedfunc(a,c,x):
    return a * numpy.exp((x) / (c))


x = Xdata
y = fit_result['y']

##lines = plt.plot(Xdata,Ydata,x,y)
l1,l2 = ax.lines

plt.setp(l1,'marker','o','markersize',4,'linestyle','none','markeredgecolor','none','color',SteelBlue)
plt.setp(l2,'marker','','linestyle','-','linewidth',2,'color',Maroon)

ax.set_title('Herotek Calibration')
ax.set_ylabel('Herotek Voltage (V)')
ax.set_xlabel('Signal Generator Power (dBm)')
ax.grid(True)
plt.show()