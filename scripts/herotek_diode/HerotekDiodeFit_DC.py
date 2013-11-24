import os,sys,time

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

# set-up file load
basefile_dir = 'C:\Users\Felix\Dropbox\Awschalom Lab Share\Python Codes'
os.chdir(basefile_dir)
calibfilename = '203918_herotek_calibration'
file_name =('%s.dat' % calibfilename)
filepath = os.path.join(basefile_dir,file_name)
fcol = 0
xcol=1
ycol=2

# load file
d = numpy.loadtxt(filepath)
rows, cols = d.shape
xlen = 1
while d[xlen,xcol] != d[0,xcol]:
    xlen+=1
flen = d.shape[0]/xlen

# declare variables
fvals = numpy.zeros(flen)
xvals = numpy.zeros((flen,xlen))
yvals = numpy.zeros((flen,xlen))

# separate values
for ii in range (0,flen):
    fvals[ii] = d[ii*xlen,fcol]
    for jj in range(0,xlen):
        xvals[ii,jj] = d[ii*xlen+jj][xcol]
        yvals[ii,jj] = d[ii*xlen+jj][ycol]

# Colors
SteelBlue = '#4682b4'
Maroon = '#800000'

# fit function

# fit data
def fitDiodeData(x,y):
    o = fit.Parameter(0,'offset')
    a = fit.Parameter(-0.2, 'amplitude')
    tau = fit.Parameter(7, 'tau')
    fitfunc_str = 'y(x) = o + a*exp(x/tau)'

    def fitfunc(x):
        return o() + a() * numpy.exp((x) / (tau()))

    fit_result = fit.fit1d(x,y, None, p0=[o,tau,a], fixed=[], fitfunc=fitfunc,
            fitfunc_str=fitfunc_str, do_print=True, ret=True)
    return fit_result


# select parameters to fit
freq = fvals
power = xvals


plt.ion()
plt.show()
fig = plt.figure(figsize=[12,6])
ax = fig.add_subplot(111)

ampl = numpy.zeros(len(freq))
offset = numpy.zeros(len(freq))
tau = numpy.zeros(len(freq))

##for ff in range (0,len(fvals)):
for ff in range (0,len(freq)):
    x = xvals[ff,:]
    y = yvals[ff,:]

    # set up plot

    # fit diode data
    fit_result = fitDiodeData(x,y)

    plot.plot_fit1d(fit_result, numpy.linspace(x[0],x[-1],201), ax=ax, ret = 'ax',
        plot_data=True)

    l1,l2 = ax.lines

    ampl[ff] = fit_result['params'][2]
    offset[ff] = fit_result['params'][0]
    tau[ff] = fit_result['params'][1]

    plt.setp(l1,'marker','o','markersize',4,'linestyle','none','markeredgecolor','none','color',SteelBlue)
    plt.setp(l2,'marker','','linestyle','-','linewidth',2,'color',Maroon)

    ax.set_title(('Herotek Calibration, Freq: %s GHz' % freq[ff]))
    ax.set_ylabel('Herotek Voltage (V)')
    ax.set_xlabel('Signal Generator Power (dBm)')
    ax.grid(True)
    plt.draw()
##    plt.show()
    time.sleep(1)
    ax.cla()


base_dir = 'C:\Users\Felix\Dropbox\Awschalom Lab Share\Python Codes'
file_name = ('%s_fitparams.dat'% calibfilename)
file_name = os.path.join(base_dir,file_name)
fp = open(file_name,'w')
fp.write('Freq. \t Ampl. \t Tau \t Offs\n')
for ff in range(0,len(freq)):
    fp.write('%s\t%s\t%s\t%s\n' % (freq[ff],ampl[ff],tau[ff],offset[ff]))
fp.close()

