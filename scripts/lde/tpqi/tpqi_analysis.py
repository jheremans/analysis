import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams

from analysis.lib.fitting import fit
from analysis.lib.lde import tpqi
from analysis.lib.pq import hht3

from analysis import config

config.outputdir = r'D:\bjhensen\data\20120814_HWP_ZZ\tpqi'
config.datadir = r'D:\bjhensen\data\20120814_HWP_ZZ'

rcParams['mathtext.default'] = 'regular'
rcParams['legend.numpoints'] = 1


class Coincidences:
    '''
    Class to extract TPQI coincidences from the LDE data.
    for all deltas (i.e., the peaks in the pulsed g2 function),
    we hold coincidences in the form [[ch0 time], [ch1 time]].
    this allows for later filtering on the times (starts, windowlengths).
    '''

    def __init__(self):
        self.savedir = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-tpqi')
        # self.ldesrcdir = '/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde/X-X/2012-10-06_X-X' 
        #self.ldesrcdir = '/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde'
        self.ldesrcdir = r'D:\bjhensen\data\20120814_HWP_ZZ'
        self.windowstart = (640, 670)
        self.windowlength = 150
        self.fnpattern = 'hhp_data'

        self.deltas = [-1,0,1]
        self.coincidences = [ np.zeros((0,2)) for d in self.deltas ]

    def load_coincidences_from_events(self):
        for (path,dirs,files) in os.walk(self.ldesrcdir):
            for fn in files:
                if self.fnpattern in fn and fn[0] != '.' and fn[-4:] == '.npz':
                    
                    print os.path.join(path,fn)
                    
                    d = np.load(os.path.join(path, fn))
                    data = d['hhdata']
                    d.close()

                    data = hht3.filter_timewindow(data, 0, 
                            mintime=0,
                            maxtime=2400)
                    data = hht3.filter_timewindow(data, 1,
                            mintime=0,
                            maxtime=2400)
                    
                    d,c = tpqi.coincidence_times(data, self.deltas)
                    self.coincidences = tpqi.add_coincidence_times(
                            self.coincidences, c)

    def save_coincidences(self):
        if not os.path.isdir(self.savedir):
            os.makedirs(self.savedir)

        np.savez(os.path.join(self.savedir, 'coincidences'),
                deltas=self.deltas,
                coincidences=self.coincidences)

    def load_coincidences(self, folder):
        d = np.load(os.path.join(config.outputdir, folder, 'coincidences.npz'))
        self.coincidences = d['coincidences']
        self.deltas = d['deltas']

    def filtered_coincidences(self):
        '''
        returns the time-filtered coincidences.
        for subsequent processing, only dt is needed (instead of ch0 time, and
        ch1 time). therefore, coincidences have the format [dt] for each delta.
        '''

        dts = []
        for i,d in enumerate(self.deltas):
            _c = self.coincidences[i]
            good0 = np.logical_and(_c[:,0]>=self.windowstart[0],
                    _c[:,0]<self.windowstart[0]+self.windowlength)
            good1 = np.logical_and(_c[:,1]>=self.windowstart[1],
                    _c[:,1]<self.windowstart[1]+self.windowlength)
            good = np.logical_and(good0,good1)
            dts.append(_c[good,0]-_c[good,1])
        return dts


class TPQIAnalysis:

    def __init__(self):
        self.savedir = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-tpqi')
        self.binwidth = 0.256
        self.pulsesep = 600.
        self.offset = 30.
        self.central = 1
        self.tau = 12.1 / self.binwidth
        self.binedges = np.append(0, np.arange(2,501,2)) # np.array([0,20,50,100,150])
        self.meandts = self.binedges[:-1] + (self.binedges[1:]-self.binedges[:-1])/2.
        self.cumulative = True
        
        self.deltas = [-1,0,1]
        self.coincidences = [ np.array([]) for d in self.deltas ]
    
    def g2_hist(self, range=(-520,520), bins=2*52):
        # NOTE bins should be an even number!
        
        self.peaks = {}
        self.amplitudes = []
        self.normpeaks = {}

        for i,d in enumerate(self.deltas):
            hy,hx = np.histogram(self.coincidences[i]+self.offset, 
                    bins=bins, range=range)
            self.peaks[d] = (hy,hx)
            self.fitdt = hx[1]-hx[0]

            if d != 0:
                A = fit.Parameter(max(hy))

                def fitfunc(x):
                    return A()*np.exp(-abs(x)/self.tau)

                fit.fit1d(hx[:-1]+(hx[1]-hx[0])/2., hy, None, fitfunc=fitfunc, p0=[A], do_print=True, 
                        fixed=[])

                self.amplitudes.append(A())
                    
        # normalize peaks
        self.meanamp = np.mean(self.amplitudes)
        for d in self.peaks:
            self.normpeaks[d] = (self.peaks[d][0]/self.meanamp, self.peaks[d][1])
                       
        return True   

    
    # NOTE here the deltas are actually hardcoded; fine for the LDE experiment
    def plot_g2(self, fits=True):
        
        self.g2fig,(ax1,ax2,ax3) = plt.subplots(1,3,sharey=True,figsize=(12,8))
              
        for d,ax in zip([-1,0,1], (ax1,ax2,ax3)):
            ax.plot(self.normpeaks[d][1][1:],self.normpeaks[d][0],
                    'k', drawstyle='steps')
            if fits:
                x = np.linspace(-300,300,101)
                if d == 0:                    
                    # model for the shape of the central peak (from Bas,
                    # according to Legero et al.)
                    A = fit.Parameter(2.)
                    dw = fit.Parameter(0.04)
                    def ff(x):
                        return A()*np.exp(-abs(x)/self.tau)*\
                                (1.-np.exp(-0.25*dw()**2*x**2))

                    dx = self.normpeaks[0][1][1]-self.normpeaks[0][1][0]
                    
                    res=fit.fit1d(self.normpeaks[d][1][:-1]+dx/2., 
                            self.normpeaks[d][0], None, fitfunc=ff, p0=[A,dw],
                            do_print=True, fixed=[0],ret=True)
                    
                    ax.plot(x, 2.*np.exp(-abs(x)/self.tau), 'r--', lw=2)
                    ax.plot(x, ff(x), 'r-', lw=2)                    
                    ax.text(-290,1.8, ('$A = %.1f$ (fixed)'+'\n'+\
                            r'$\delta\omega = 2\pi\times(%.0f\pm%.0f)$ MHz') % \
                            (A(), dw()*0.256*1e3, res['error'][0]*0.256*1e3),
                            color='r')
                    
                else:
                    ax.plot(x, np.exp(-abs(x)/self.tau), 'r-', lw=2)

        ax1.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax1.yaxis.tick_left()
        ax3.yaxis.tick_right()
        ax2.yaxis.set_ticks_position('none')
        ax1.set_ylabel(r'$g^{(2)}(\tau)$')
        ax2.set_xlabel(r'$\tau$ (bins)')
        plt.subplots_adjust(wspace=0.15)

    def analysis(self, crratio=1.1):
        
        self.visibility = np.array([])
        self.u_visibility = np.array([])
        self.fidelity = np.array([])
        self.u_fidelity = np.array([])
        self.counts = np.array([])

        zp = self.peaks[0]
        amp = 2*self.meanamp 
        dip = 2*2*crratio/(1+crratio)**2
        
        # rebin the counts for more appropriate resolution of the visibility
        # use these also for the fidelity measurements
        self.zprebinned, _be = np.histogram(
                np.array(self.coincidences[self.central]), 
                bins=np.append(-self.binedges[:0:-1], self.binedges))

        length = len(self.meandts)

        # NOTE formulas see labbook
        for i in range(0,length):

            if self.cumulative:
                tpqicts = np.sum(self.zprebinned[length-i-1:length+i+1])
            else:
                tpqicts = np.sum(self.zprebinned[length+i:length+i+1]) \
                        + np.sum(self.zprebinned[length-i-1:length-i])
            
            def notpqicts(x0,x1):
                return -2*dip*amp/self.fitdt*self.tau*(np.exp(-abs(x1)/self.tau)\
                        - np.exp(-abs(x0)/self.tau)) # factor two: absolute!
            
            if self.cumulative:
                n = notpqicts(0, self.binedges[i+1])
            else:
                n = notpqicts(self.binedges[i], self.binedges[i+1])
            
            v = (n-tpqicts)/n
            u_v = np.sqrt(tpqicts)/n
            F = 0.5*(1+v)
            u_F = 0.5*u_v

            self.visibility = np.append(self.visibility, v)
            self.u_visibility = np.append(self.u_visibility, u_v)
            self.fidelity = np.append(self.fidelity, F)
            self.u_fidelity = np.append(self.u_fidelity, u_F)
            self.counts = np.append(self.counts, tpqicts)        
        
        ### plot visibility and projected fidelity 
        self.fidfig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(self.binedges, np.append(self.visibility[0],self.visibility), 
                'k', drawstyle='steps')
        ax.errorbar(self.meandts, self.visibility, fmt='o', mec='k', mfc='w',
                yerr=self.u_visibility, ecolor='k')

        ax.set_xlabel(r'$| \tau |$ (bins)')
        ax.set_ylabel('visibility')
        ax.set_xlim(0,self.binedges.max())
        ax.set_ylim(0,1.1)

        ax2 = ax.twinx()
        ax2.plot(self.binedges, np.append(self.fidelity[0], self.fidelity),
                 'r', drawstyle='steps')
        ax2.errorbar(self.meandts, self.fidelity, fmt='o', mec='r', mfc='w',
                 yerr=self.u_fidelity, ecolor='r')
        ax2.set_ylabel('fidelity', color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')

        ax2.set_xlim(0,self.binedges.max())
        ax2.set_ylim(0,1.1)
        ax2.grid(color='r')

    def save(self):
        if not os.path.isdir(self.savedir):
            os.makedirs(self.savedir)

        suffix = '_ch0%d_win%d' % (self.ch0start, self.windowlength)
        prefix = 'cumulative_' if self.cumulative else ''

        try:       
            self.g2fig.savefig(os.path.join(self.savedir,'g2'+suffix+'.png'))
        except AttributeError:
            print 'no g2 plot, skip...'

        try:       
            self.fidfig.savefig(os.path.join(self.savedir,prefix+'visibility'+suffix+'.png'))
        except AttributeError:
            print 'no fidelity plot, skip...'

        
        f = open(os.path.join(self.savedir, 'peak_histograms'+suffix+'.pkl'), 'wb')
        pickle.dump(self.peaks, f)
        f.close()

        f = open(os.path.join(self.savedir, 'peak_histograms_normalized'+suffix+'.pkl'), 'wb')
        pickle.dump(self.normpeaks, f)
        f.close()

        np.savez(os.path.join(self.savedir, prefix+'tpqi'+suffix),
                visibility=self.visibility, 
                u_visibility=self.u_visibility, 
                binedges=self.binedges,
                meandts=self.meandts,
                fidelity=self.fidelity,
                u_fidelity=self.u_fidelity,
                zprebinned=self.zprebinned,
                )


def sweep():
    c = Coincidences()
    c.load_coincidences('20121122-tpqi')

    for ch0 in range(636, 651, 2):
        for l in range(100, 301, 100):
            c.windowstart = (ch0, ch0+30)
            c.windowlength = l
            
            a = TPQIAnalysis()
            a.coincidences = c.filtered_coincidences()
            a.ch0start, a.ch1start = c.windowstart
            a.windowlength = c.windowlength
            a.cumulative = False

            a.g2_hist()
            a.plot_g2()
            a.analysis()
            a.save()

            plt.close('all')

if __name__ == '__main__':
    c = Coincidences()
    c.load_coincidences('20121122-tpqi')
    c.windowstart = (640,670)
    c.windowlength = 1000

    a = TPQIAnalysis()
    a.coincidences = c.filtered_coincidences()
    a.ch0start, a.ch1start = c.windowstart
    a.windowlength = c.windowlength
    a.cumulative = True
    a.g2_hist()
    a.plot_g2()
    a.analysis()
    a.save()
