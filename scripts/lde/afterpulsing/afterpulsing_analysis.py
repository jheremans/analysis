import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt

from analysis.lde import hht3, spcorr

class AfterPulsingAnalysis:
    
    # rawfolder = r"D:\Analysis\2012-09-21_afterpulsing\20120923\083802_LDE_test-after-pulsing\rawdata-000"
    # datfolder,_tmp = os.path.split(rawfolder)
    # prepared_fname = os.path.join(datfolder, 'hhdata.npz')

    def __init__(self):
        self.srcfolder = '/Volumes/MURDERHORN/TUD/LDE/data/afterpulsing/20120927'
        self.fnpattern = 'hh_data_3windows-'
        self.rawpattern = 'rawdata-'
        self.suffixlen = 3
        self.savefolder = '.'

    def prepare_hhdata(self):

        print 'start processing rawdata'

        for (path,dirs,files) in os.walk(self.srcfolder):
            for dn in dirs:
                if self.rawpattern in dn:
                    
                    print 'processing', os.path.join(path,dn)
                    
                    suffix = dn[-self.suffixlen:]
                    d,w1,w2,w3 = spcorr.load_hhdata(
                            os.path.join(path,dn), windows=[1,2,3],
                            ch0maxtime=2400, ch1maxtime=2400, 
                            do_filter_crap=False, print_status=True)

                    np.savez(os.path.join(path,self.fnpattern+suffix),
                            hhdata=d, w1=w1, w2=w2, w3=w3)

    
    def process_data(self):
        
        self.ap = [ np.zeros((0,4), dtype=np.uint) for i in range(2) ]
        self.noap = np.zeros((0,4), dtype=np.uint)

        for (path,dirs,files) in os.walk(self.srcfolder):
            for fn in files:
                if self.fnpattern in fn and fn[:1] != '.':
                    fp = os.path.join(path,fn)
                    print 'process file', fp

                    d = np.load(fp)
                    hhdata = d['hhdata']
                    w1 = d['w1']
                    w2 = d['w2']
                    w3 = d['w3']
                    d.close()

                    ap, noap = self.get_ap([w1,w2,w3], [850, 850+2344])
                    self.noap = np.vstack((self.noap, noap))
                    for i,a in enumerate(ap):
                        self.ap[i] = np.vstack((self.ap[i], a))

        np.savez(os.path.join(self.savefolder, 'clicks.npz'), noap=self.noap,
                ap1 = self.ap[0], ap2 = self.ap[1])
        np.savez(os.path.join(self.savefolder, 'click_times_ch1.npz'), 
                noap=hht3.get_click_times(self.noap)[1],
                ap1=hht3.get_click_times(self.ap[0])[1],
                ap2=hht3.get_click_times(self.ap[1])[1])

    def load_data(self, clicksfn='clicks.npz'):
        f = np.load(clicksfn)
        self.ap = []
        self.ap.append(f['ap1'])
        self.ap.append(f['ap2'])
        self.noap = f['noap']
        f.close()


    def get_stats(self):
        stats = {
                'seqstarts' : 0,
                }
         
        for (path,dirs,files) in os.walk(self.srcfolder):
            for fn in files:
                if 'LDE-' in fn and fn[:1] != '.' and fn[-4:] == '.npz':
                    fp = os.path.join(path,fn)
                    print 'process file', fp

                    d = np.load(fp)
                    d2 = d['adwin_lt2_pars']
                    startsidx = d2[:,0] == 'get_noof_seq_starts'
                    stats['seqstarts'] += int(d2[startsidx,1][0])
                    d.close()

        pprint.pprint(stats)
        f = open(os.path.join(self.savefolder, 'stats.pkl'), 'wb')
        pickle.dump(stats, f)
        f.close()

    def count_histogram(self, ch=1, rng=(0,3*2400), bins=(3*40+1)):
        t_noap = hht3.get_click_times(self.noap)[ch]
        t_ap = [ hht3.get_click_times(ap)[ch] for ap in self.ap ]

        apcolors = ['r', 'b', 'g']
        aphatches = ['', '', '']
        aplabels = ['AP pulse 1', 'AP pulse 2', 'AP pulse 3']
        apalphas = [ 0.5, 0.5, 0.5 ]

        fig = plt.figure()
        ax = plt.subplot(111)
        
        ax.hist(t_noap, bins=bins, range=rng, log=True, label='no AP',
                color='w', hatch='/')
        
        for t,c,h,l,a in zip(t_ap, apcolors, aphatches, aplabels, apalphas):
            ax.hist(t, bins=bins, label=l, range=rng, color=c,
                    hatch=h, alpha=a)

        ax.legend()

        fig.savefig(os.path.join(self.savefolder, 'click_histogram.pdf'),
            format='pdf')

        return
    
    def get_ap(self, windows, thresholdtimes, dt=2343.75):
        photons = np.zeros((0,4), dtype=np.uint)
        for i,w in enumerate(windows):
            wphotons = w[w[:,3]==0]
            wphotons[:,1] += int(i*dt+0.5)
            wphotons[:,0] -= i
            photons = np.vstack((photons,wphotons))

        noap = photons.copy()
        ap = []
        for i,t in enumerate(thresholdtimes):
            apguess = noap[noap[:,1]>t]
            laserguess = noap[noap[:,1]<t]
            realap = np.logical_and(np.logical_and(np.in1d(noap[:,0], laserguess[:,0]),
                np.in1d(noap[:,0], apguess[:,0])), noap[:,1]>t)
            ap.append(noap[realap])
            noap = noap[logical_not(realap)]

        return ap, noap


if __name__ == '__main__':
    a = AfterPulsingAnalysis()
    #a.prepare_hhdata()
    a.load_data()
    a.count_histogram()
    # a.get_stats()


