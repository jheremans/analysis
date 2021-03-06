import os, sys, time
import pickle
import pprint
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import ImageGrid

from analysis.lib.fitting import fit
from analysis.lib.lde import sscorr, fidelities
from analysis import config

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'
config.datadir = r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde'

rcParams['mathtext.default'] = 'regular'
rcParams['legend.numpoints'] = 1


hhpfilebase = 'hhp_data'
chmaxtime = 2300

# columns of the data field
T1 = 0
T2 = 1
CH1 = 2
CH2 = 3
RO1 = 4
RO2 = 5
CR1 = 6
CR2 = 7
GP = 8

# indices used for bases
ZZidx = 0
XXidx = 1
XmXidx = 2

class EntanglementEventAnalysis:
    """
    use to get possible entanglement events from hhp-data files as 
    generated by the sscorr-module. can then filter those events
    based on allowed detection windows.
    """

    def __init__(self, name=''):
        self.events = np.zeros((0,9))
        self.savedir = os.path.join(config.outputdir, 
                time.strftime('%Y%m%d')+'-lde')
        self.name = name

    def add_events(self, fp, start=630, stop=880):
        idx = os.path.splitext(fp)[0][-3:]
        hhpf = np.load(fp)
        hhp = hhpf['hhp']
        hhpf.close()

        folder,fn = os.path.split(fp)
        ssro1, ssro2, cr1, cr2, gate = sscorr.load_adwin_data(
                folder, DATIDX=int(idx))

        photons = np.logical_and(np.logical_and(hhp[:,3] == 0, 
            hhp[:,1] >= start), hhp[:,1] <= stop)
        mrkrs = hhp[:,3] == 1
        for _i, nsync in np.ndenumerate(hhp[hhp[:,3]==1,0]):
            i = _i[0]
            w1 = hhp[np.logical_and(photons, hhp[:,0] == nsync-1)]
            w2 = hhp[np.logical_and(photons, hhp[:,0] == nsync)]           
            
            if len(w1) == 1 and len(w2) == 1:
                self.events = np.vstack((self.events, 
                    np.array([w1[0,1], w2[0,1], w1[0,2], w2[0,2], 
                        ssro1[i], ssro2[i], cr1[i], cr2[i], gate[i]])))

    def save_events(self, fn):
        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        np.savez(os.path.join(self.savedir,fn), events=self.events)

    def load_events(self, folder, fn):
        f = np.load(os.path.join(config.outputdir, folder, fn))
        self.events = f['events']
        f.close()

    def get_entanglement_events(self, folder):
        for (path, dirs, files) in os.walk(
                os.path.join(config.datadir,folder)):
            for fn in files:
                if hhpfilebase in fn and fn[0] != '.':
                    self.add_events(os.path.join(path, fn))

        print 'found %d events in total' % len(self.events)

    def correlations(self, plot=True, save=False, failsafe=True, **kw):
        """
        gives correlation counts in the form [ms=00, 01, 10, 11]
        """
        self.psi1= np.zeros(4, dtype=int)
        self.psi2= np.zeros(4, dtype=int)
        for e in self.events:
            idx = int(e[RO1] == 0)*2**1 + int(e[RO2] == 0)*2**0
            if e[CH1] == e[CH2]:
                self.psi1[idx] += 1
            else:
                self.psi2[idx] += 1

        if failsafe and (self.psi1.sum() == 0 or self.psi2.sum() == 0):
            return False

        self.readout_correction(**kw)

        if save:
            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)
            np.savez(os.path.join(self.savedir, self.name+'_correlations'), 
                    psi1=self.psi1, psi2=self.psi2,
                    psi1_corrected=self.psi1_corrected, 
                    u_psi1_corrected=self.u_psi1_corrected,
                    psi2_corrected=self.psi2_corrected, 
                    u_psi2_corrected=self.u_psi2_corrected)
        
        if plot:
            self.plot_correlations()

        return True
        
    def plot_correlations(self, save=False):
        fig,[[ax1,ax2],[ax3,ax4]] = plt.subplots(2,2,figsize=(8,8))

        for ax,state,name in zip([ax1,ax2], [self.psi1, self.psi2], 
                ['Psi1', 'Psi2']):

            ax.bar(np.arange(4)-0.4, state, color='w', ec='k', width=0.8)
            ax.set_xticks(np.arange(4))
            ax.set_xticklabels(['00', '01', '10', '11'])
            ax.set_xlabel('State')
            ax.set_ylabel('Occurences')
            ax.set_xlim(-0.5,3.5)
            ax.set_title(name)

        for ax,state,err,name in zip([ax3,ax4], 
                [self.psi1_corrected, self.psi2_corrected],
                [self.u_psi1_corrected, self.u_psi2_corrected], 
                ['Psi1 corrected', 'Psi2 corrected']):
           
            ax.bar(np.arange(4)-0.4, state, color='w', ec='k', width=0.8,
                    yerr=err, ecolor='k')
            ax.set_xticks(np.arange(4))
            ax.set_xticklabels(['00', '01', '10', '11'])
            ax.set_xlabel('State')
            ax.set_ylabel('Probability')
            ax.set_xlim(-0.5,3.5)
            ax.set_title(name)
        
        plt.suptitle(self.name)
        plt.tight_layout()
        
        if save:
            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)
            fig.savefig(os.path.join(self.savedir, 
                self.name+'_correlations.png'))

    def filter_times(self, mintimes=(640,670), window1=200, window2=200, dtmin=0, dtmax=50):
        offset = (0,mintimes[1]-mintimes[0])
        
        badidxs = []
        for i,e in enumerate(self.events):
            ch1 = int(e[CH1])
            ch2 = int(e[CH2])
            if e[T1] < mintimes[ch1] or e[T1] > mintimes[ch1]+window1:
                badidxs.append(i)
                continue

            if e[T2] < mintimes[ch2] or e[T2] > mintimes[ch2]+window2:
                badidxs.append(i)
                continue
            
            dt = abs(e[T2] - offset[ch2] - e[T1] + offset[ch1])
            if dt < dtmin or dt > dtmax:
                badidxs.append(i)
                continue

        self.events = np.delete(self.events, badidxs, 0)

    def filter_gate_phase(self):
        self.events = self.events[self.events[:,8]>0]

    def readout_correction(self, F0_1=0.917, uF0_1=0.003, F1_1=0.993, 
            uF1_1=0.001, F0_2=0.809, uF0_2=0.007, F1_2=0.985, uF1_2=0.011):
        
        # print ''
        # print 'psi1:', self.psi1[::-1]
        p1, up1 = sscorr.ssro_correct_twoqubit_state_photon_numbers(
                self.psi1[::-1], F0_1, F0_2, F1_1, F1_2, dF0a=uF0_1, dF1a=uF1_1,
                dF0b=uF0_2, dF1b=uF1_2, return_error_bars=True, verbose=False)
        self.psi1_corrected = p1[::-1].reshape(-1)
        self.u_psi1_corrected = up1[::-1]

        # print 'psi2:', self.psi2[::-1]
        p2,up2 = sscorr.ssro_correct_twoqubit_state_photon_numbers(
                self.psi2[::-1], F0_1, F0_2, F1_1, F1_2, dF0a=uF0_1, dF1a=uF1_1,
                dF0b=uF0_2, dF1b=uF1_2, return_error_bars=True, verbose=False)
        self.psi2_corrected = p2[::-1].reshape(-1)
        self.u_psi2_corrected = up2[::-1]


class FidelityAnalysis:
    """
    Walks through parameter space spanned by the channel starts, window 
    length and maximum dt and calculates the fidelity for each point.
    can save correlations and fidelity.
    can plot color maps of the fidelity, vs parameters, incl information on
    statistical significance.
    """
 
    def __init__(self, name=''):
        self.savedir = os.path.join(config.outputdir, 
                time.strftime('%Y%m%d')+'-ldefidelity')
        self.name = name

        self.mode = None
        self.dtslices = False
       
        self.psi1_ZZ = None
        self.psi1_XX = None
        self.psi1_XmX = None
        self.psi2_ZZ = None
        self.psi2_XX = None
        self.psi2_XmX = None
        
        #self.psi1_XmX_pt2 = None
        #self.psi2_XmX_pt2 = None

        ### these are the values corrected for imperfect initialization
        self.F0_1 = 0.921; self.u_F0_1 = 0.003
        self.F1_1 = 0.997; self.u_F1_1 = 0.001
        self.F0_2 = 0.822; self.u_F0_2 = 0.008
        self.F1_2 = 0.989; self.u_F1_2 = 0.010

        self.F1_2_upper = 0.998; self.F1_2_lower = 0.962
        
        # self.dtvals = np.arange(11).astype(int)*20 + 20
        self.dtvals = np.arange(10,151,2) # np.array([20,50,100,150], dtype=int)
        self.win1vals = np.array([150]) # np.arange(5).astype(int)*50 + 100
        self.win2vals = np.array([75,150]) # np.arange(11).astype(int)*25 + 50
        self.ch0start = 640
        # self.ch0starts = np.arange(3).astype(int)*2 + 636 
        # self.ch1starts = np.arange(3).astype(int)*2 + 668
        self.ch1offset = 30
        self.dtmin = 0
        
    def _get_corrected_correlations(self, correlations, mod, *rofidelities):
        _c,_u_c = sscorr.ssro_correct_twoqubit_state_photon_numbers(
                correlations[::-1], *rofidelities, verbose=False)
        """
        perform readout correction on a histogram, with given SSRO fidelities
        and uncertainties. depending on mode, enforces positive probabilities,
        keeping the sum at 1 at the expense of others (depends on mode).
        """
        c = _c[::-1].reshape(-1)
        u_c = _u_c[::-1].reshape(-1)

        def mod_none(c,u_c):
            return c, u_c

        def mod_nonegatives_compensate_proportionally(c,u_c):
            neg = c < 0.
            total = c[neg].sum()
            c[np.logical_not(neg)] /= (1.+float(total))
            c[neg] = 0.
            return c, u_c
        
        # NOTE not universal; would not be too hard to rewrite, but since not
        # necessary, i'm too lazy now; -wolfgang
        def mod_nonegatives_compensate_opposite_parity(c,u_c):
            if c[3] < 0:
                c[1] -= abs(c[3]/2.)
                c[2] -= abs(c[3]/2.)
                c[3] = 0
            return c, u_c

        modifiers = {
                None : mod_none,
                'lowerbound' : mod_nonegatives_compensate_opposite_parity,
                'bestguess' : mod_none,
                }

        return modifiers[mod](c,u_c)
    
    def _fidelity_from_hists(self, zzhist, xxhist, xmxhist, state, mod, *rofids):
        
        if mod != 'raw' and mod != 'rawlowerbound' and mod != 'rawbestguess':
            c_zz, u_c_zz = self._get_corrected_correlations(zzhist, mod, *rofids)
            c_xx, u_c_xx = self._get_corrected_correlations(xxhist, mod, *rofids)
            c_xmx, u_c_xmx = self._get_corrected_correlations(xmxhist, mod, *rofids)
        else:
            c_zz = zzhist.astype(float)/zzhist.sum()
            u_c_zz = c_zz*(c_zz.sum()-c_zz)/c_zz.sum()
            c_xx = xxhist.astype(float)/xxhist.sum()
            u_c_xx = c_xx*(c_xx.sum()-c_xx)/c_xx.sum()
            c_xmx = xmxhist.astype(float)/xmxhist.sum()
            u_c_xmx = c_xmx*(c_xmx.sum()-c_xmx)/c_xmx.sum()

        if state == 'psi1':
            self.psi1_zz_corrected, self.u_psi1_zz_corrected = c_zz, u_c_zz
            self.psi1_xx_corrected, self.u_psi1_xx_corrected = c_xx, u_c_xx
            self.psi1_xmx_corrected, self.u_psi1_xmx_corrected = c_xmx, u_c_xmx
        elif state == 'psi2':    
            self.psi2_zz_corrected, self.u_psi2_zz_corrected = c_zz, u_c_zz
            self.psi2_xx_corrected, self.u_psi2_xx_corrected = c_xx, u_c_xx
            self.psi2_xmx_corrected, self.u_psi2_xmx_corrected = c_xmx, u_c_xmx

        zz_term = c_zz[[1,2]].sum()
        u_zz_term = fidelities.ro_c_dF(zzhist[::-1], *rofids)
        
        _s1 = c_zz[0] if c_zz[0] > 0. else 0.
        _s2 = c_zz[3] if c_zz[3] > 0. else 0.
        sqrt_term = np.sqrt(_s1*_s2)
        u_sqrt_term = fidelities.ro_c_dF_S(zzhist[::-1], *rofids)
        
        _xx = c_xx[0]+c_xx[3] if state=='psi1' else c_xx[1]+c_xx[2]
        _u_xx = fidelities.ro_c_dF(xxhist[::-1], *rofids)
        _xmx = c_xmx[1]+c_xmx[2] if state=='psi1' else c_xmx[0]+c_xmx[3]
        _u_xmx = fidelities.ro_c_dF(xmxhist[::-1], *rofids)
        nxx = float(sum(xxhist))/sum(xxhist+xmxhist)
        nxmx = float(sum(xmxhist))/sum(xxhist+xmxhist)
        xx_term = nxx*_xx + nxmx*_xmx
        u_xx_term = np.sqrt((nxx*_u_xx)**2 + (nxmx*_u_xmx)**2)

        if mod=='bestguess' or mod=='rawbestguess':
            sqrt_term = 0.
            u_sqrt_term = 0.
        
        F = zz_term/2. - sqrt_term + (xx_term - 0.5)
        u_F = np.sqrt(.25*u_zz_term**2 + u_sqrt_term**2 + u_xx_term**2)

        return F, u_F
    
    def _get_fidelity(self, zzhist, xxhist, xmxhist, state='psi1', mod=None):
        rofids_regular = [self.F0_1, self.F0_2, self.F1_1, self.F1_2, 
                self.u_F0_1, self.u_F0_2, self.u_F1_1, self.u_F1_2]
        
        rofids_extremal1 = [self.F0_1, self.F0_2, self.F1_1, self.F1_2_upper, 
                self.u_F0_1, self.u_F0_2, self.u_F1_1, 0.]

        rofids_extremal2 = [self.F0_1, self.F0_2, self.F1_1, self.F1_2_lower, 
                self.u_F0_1, self.u_F0_2, self.u_F1_1, 0.]

        if mod == None or mod == 'bestguess' or mod == 'raw' or mod == 'rawbestguess':
            return self._fidelity_from_hists(zzhist, xxhist, xmxhist, 
                    state, mod, *rofids_regular)
        
        if mod == 'lowerbound' or mod == 'rawlowerbound':
            F1,u_F1 = self._fidelity_from_hists(zzhist, xxhist, xmxhist,
                    state, mod, *rofids_extremal1)
            F2,u_F2 = self._fidelity_from_hists(zzhist, xxhist, xmxhist,
                    state, mod, *rofids_extremal2)

            if F1 < F2:
                # this is to set the correlations back to the ones
                # corresponding to the lower fidelity
                F1,u_F1 = self._fidelity_from_hists(zzhist, xxhist, xmxhist,
                    state, mod, *rofids_extremal1)
                return F1,u_F1 
            else: 
                return F2,u_F2

    def fidelity(self, mod=None):       
        self.F_psi1, self.u_F_psi1 = self._get_fidelity(
                self.psi1_ZZ, self.psi1_XX, self.psi1_XmX, 
                state='psi1', mod=mod)
        
        self.F_psi2, self.u_F_psi2 = self._get_fidelity(
                self.psi2_ZZ, self.psi2_XX, self.psi2_XmX, 
                state='psi2', mod=mod)

    def get_fidelities(self, folder='20121030-lde', 
            ZZfn='ZZ.npz', XXfn='XX.npz', XmXfn='X-X.npz'):

        eZZ = EntanglementEventAnalysis('ZZ')
        eZZ.load_events(folder, ZZfn)
        eZZ.filter_gate_phase()
        ZZ0 = eZZ.events.copy()

        eXX = EntanglementEventAnalysis('XX')
        eXX.load_events(folder, XXfn)
        eXX.filter_gate_phase()
        XX0 = eXX.events.copy()
        
        eXmX = EntanglementEventAnalysis('X-X')
        eXmX.load_events(folder, XmXfn)
        eXmX.filter_gate_phase()
        XmX0 = eXmX.events.copy()
        
        self.psi1fids = np.zeros((len(self.dtvals), len(self.win1vals), 
            len(self.win2vals)))
        self.u_psi1fids = np.zeros((len(self.dtvals), len(self.win1vals), 
            len(self.win2vals)))
        self.psi2fids = np.zeros((len(self.dtvals), len(self.win1vals), 
            len(self.win2vals)))
        self.u_psi2fids = np.zeros((len(self.dtvals), len(self.win1vals), 
            len(self.win2vals)))

        self.rawpsi1correlations = np.zeros((len(self.dtvals), 
            len(self.win1vals), len(self.win2vals), 3, 4))
        self.rawpsi2correlations = np.zeros((len(self.dtvals), 
            len(self.win1vals), len(self.win2vals), 3, 4))
        self.correctedpsi1correlations = np.zeros((len(self.dtvals), 
            len(self.win1vals), len(self.win2vals), 3, 4))
        self.correctedpsi2correlations = np.zeros((len(self.dtvals), 
            len(self.win1vals), len(self.win2vals), 3, 4))
        self.u_correctedpsi1correlations = np.zeros((len(self.dtvals), 
            len(self.win1vals), len(self.win2vals), 3, 4))
        self.u_correctedpsi2correlations = np.zeros((len(self.dtvals), 
            len(self.win1vals), len(self.win2vals), 3, 4))

        ch0start = self.ch0start
        ch1start = self.ch0start + self.ch1offset
        
        cnt = 0
        for i,dt in enumerate(self.dtvals):
            for j,window1 in enumerate(self.win1vals):
                
                cnt += 1
                print 'chunk', cnt, '/', len(self.win1vals)*len(self.dtvals)

                for k,window2 in enumerate(self.win2vals):

#                 for k,ch0start in enumerate(self.ch0starts):
#                     cnt += 1
#                     print 'chunk', cnt, '/', \
#                             len(self.dtvals)*len(self.winvals)*len(self.ch0starts)
# 
#                     for l,ch1start in enumerate(self.ch1starts):

                    if self.dtslices and i>0:
                        dtmin = self.dtvals[i-1]
                    else:
                        dtmin = self.dtmin

                    eZZ.events = ZZ0.copy()
                    eZZ.filter_times(dtmax=dt, dtmin=dtmin, window1=window1, 
                            window2=window2, mintimes=(ch0start,ch1start))
                    eXX.events = XX0.copy()
                    eXX.filter_times(dtmax=dt, dtmin=dtmin, window1=window1, 
                            window2=window2, mintimes=(ch0start,ch1start))
                    eXmX.events = XmX0.copy()
                    eXmX.filter_times(dtmax=dt, dtmin=dtmin, window1=window1, 
                            window2=window2, mintimes=(ch0start,ch1start))

                    # readout correction goes wrong when the sum of
                    # coincidences is zero. we check this here.
                    _zzcorr = eZZ.correlations(plot=False, save=False,
                            failsafe=True)
                    _xxcorr = eXX.correlations(plot=False, save=False,
                            failsafe=True)
                    _xmxcorr = eXmX.correlations(plot=False, save=False,
                            failsafe=True)
                    if not _zzcorr or not _xxcorr or not _xmxcorr:
                        continue
                    
                    self.psi1_ZZ = eZZ.psi1
                    self.psi1_XX = eXX.psi1
                    self.psi1_XmX = eXmX.psi1

                    self.psi2_ZZ = eZZ.psi2
                    self.psi2_XX = eXX.psi2
                    self.psi2_XmX = eXmX.psi2

                    self.rawpsi1correlations[i,j,k,ZZidx] = self.psi1_ZZ
                    self.rawpsi1correlations[i,j,k,XXidx] = self.psi1_XX
                    self.rawpsi1correlations[i,j,k,XmXidx] = self.psi1_XmX
                    self.rawpsi2correlations[i,j,k,ZZidx] = self.psi2_ZZ
                    self.rawpsi2correlations[i,j,k,XXidx] = self.psi2_XX
                    self.rawpsi2correlations[i,j,k,XmXidx] = self.psi2_XmX
                    
                    self.fidelity(mod=self.mode)                
                    self.psi1fids[i,j,k] = self.F_psi1
                    self.u_psi1fids[i,j,k] = self.u_F_psi1
                    self.psi2fids[i,j,k] = self.F_psi2
                    self.u_psi2fids[i,j,k] = self.u_F_psi2
                    
                    self.correctedpsi1correlations[i,j,k,ZZidx] = self.psi1_zz_corrected
                    self.correctedpsi1correlations[i,j,k,XXidx] = self.psi1_xx_corrected
                    self.correctedpsi1correlations[i,j,k,XmXidx] = self.psi1_xmx_corrected
                    self.correctedpsi2correlations[i,j,k,ZZidx] = self.psi2_zz_corrected
                    self.correctedpsi2correlations[i,j,k,XXidx] = self.psi2_xx_corrected
                    self.correctedpsi2correlations[i,j,k,XmXidx] = self.psi2_xmx_corrected

                    self.u_correctedpsi1correlations[i,j,k,ZZidx] = self.u_psi1_zz_corrected
                    self.u_correctedpsi1correlations[i,j,k,XXidx] = self.u_psi1_xx_corrected
                    self.u_correctedpsi1correlations[i,j,k,XmXidx] = self.u_psi1_xmx_corrected
                    self.u_correctedpsi2correlations[i,j,k,ZZidx] = self.u_psi2_zz_corrected
                    self.u_correctedpsi2correlations[i,j,k,XXidx] = self.u_psi2_xx_corrected
                    self.u_correctedpsi2correlations[i,j,k,XmXidx] = self.u_psi2_xmx_corrected

                                           
    def save_fidelities(self):
        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        
        suffix = '' if self.mode == None else '_'+self.mode
        suffix += '' if self.dtslices == False else '_dtslices'
        suffix += ('_dtmin%d' % (self.dtmin)) if not self.dtslices else ''

        np.savez(os.path.join(self.savedir, 'fidelities'+suffix), 
                psi1fids=self.psi1fids,
                u_psi1fids=self.u_psi1fids, 
                psi2fids=self.psi2fids,
                u_psi2fids=self.u_psi2fids, 
                dtvals=self.dtvals, 
                win1vals=self.win1vals,
                win2vals=self.win2vals,
                ch0start=self.ch0start, 
                ch1offset=self.ch1offset)

        np.savez(os.path.join(self.savedir, 'correlations'+suffix), 
                rawpsi1correlations = self.rawpsi1correlations,
                rawpsi2correlations = self.rawpsi2correlations,
                correctedpsi1correlations = self.correctedpsi1correlations,
                correctedpsi2correlations = self.correctedpsi2correlations,
                u_correctedpsi1correlations = self.u_correctedpsi1correlations,
                u_correctedpsi2correlations = self.u_correctedpsi2correlations, 
                dtvals=self.dtvals, 
                win1vals=self.win1vals,
                win2vals=self.win2vals,
                ch0start=self.ch0start, 
                ch1offset=self.ch1offset)

    def load_fidelities(self, folder, fns=['fidelities', 'correlations']):
        self.savedir = os.path.join(config.outputdir, folder)

        suffix = '' if self.mode == None else '_'+self.mode
        suffix += '' if self.dtslices == False else '_dtslices'
        suffix += ('_dtmin%d' % (self.dtmin)) if not self.dtslices else ''

        for fn in fns:
            f = np.load(os.path.join(config.outputdir, folder, fn+suffix+'.npz'))
            for k in f.keys():
                setattr(self, k, f[k])
            f.close()

    
    def plot_fidelity_map(self, psi1slice, psi2slice, 
            xticks, yticks, xlabel='x', ylabel='y'):
        
        fig = plt.figure(figsize=(20,8))
        grid = ImageGrid(fig, 111, # similar to subplot(111)
                nrows_ncols = (2,5),
                axes_pad = 0.75,
                # add_all=True,
                label_mode = "L",
                cbar_mode = 'each',
                cbar_size="5%",
                cbar_pad="1%",
                )
        
        vmins = [0.5, 0. , None, None, None, 0.5, 0., None, None, None ]
        vmaxs = [None, None, None, None, None, None, None, None, None, None ]
        titles = ['F (psi1)', 'sigmas (psi1)', 'N', '$N p_{ZZ}(00)$', '$N p_{ZZ}(11)$',
                'F (psi2)', 'sigmas (psi2)', 'N', '$N p_{ZZ}(00)$', '$N p_{ZZ}(11)$' ]

        im0 = self.psi1fids[psi1slice]
        im1 = (self.psi1fids[psi1slice]-0.5)/self.u_psi1fids[psi1slice]
        im2 = self.rawpsi1correlations[psi1slice].sum(-1).sum(-1)
        im3 = self.rawpsi1correlations[psi1slice][...,ZZidx,0] 
        im4 = self.rawpsi1correlations[psi1slice][...,ZZidx,-1]
        
        im5 = self.psi2fids[psi2slice]
        im6 = (self.psi2fids[psi2slice]-0.5)/self.u_psi2fids[psi2slice]
        im7 = self.rawpsi2correlations[psi2slice].sum(-1).sum(-1)
        im8 = self.rawpsi2correlations[psi2slice][...,ZZidx,0] 
        im9 = self.rawpsi2correlations[psi2slice][...,ZZidx,-1]
                
        for i,im in enumerate([im0,im1,im2,im3,im4,im5,im6,im7,im8,im9]):
            img = grid[i].imshow(im, cmap=cm.gist_earth, origin='lower',
                    interpolation='nearest', vmin=vmins[i], vmax=vmaxs[i])
            fig.colorbar(img, cax=grid.cbar_axes[i])
            grid[i].set_title(titles[i])

        grid[5].set_xlabel(xlabel)
        grid[5].set_ylabel(ylabel)
        xt = grid[5].get_xticks().astype(int)[:-1]
        grid[5].set_xticklabels(xticks[xt])
        yt = grid[5].get_yticks().astype(int)[:-1]
        grid[5].set_yticklabels(yticks[yt])

        grid[0].set_yticklabels(yticks[yt])
        for i in range(5,10):
            grid[i].set_xticklabels(xticks[xt])

        return fig

    def plot_correlations(self, state, dt, win1, win2, save=True):
        bases = ['ZZ', 'XX', 'XmX']

        def idx(arr, val):
            return np.argmin(abs(arr-val))

        dtidx = idx(self.dtvals, dt)
        win1idx = idx(self.win1vals, win1)
        win2idx = idx(self.win2vals, win2)

        if state == 'psi1':
            raw = self.rawpsi1correlations[dtidx,win1idx,win2idx,:,:]
            corr = self.correctedpsi1correlations[dtidx,win1idx,win2idx,:,:]
            u_corr = self.u_correctedpsi1correlations[dtidx,win1idx,win2idx,:,:]
            fid = self.psi1fids[dtidx,win1idx,win2idx]
            u_fid = self.u_psi1fids[dtidx,win1idx,win2idx]

        else:
            raw = self.rawpsi2correlations[dtidx,win1idx,win2idx,:,:]
            corr = self.correctedpsi2correlations[dtidx,win1idx,win2idx,:,:]
            u_corr = self.u_correctedpsi2correlations[dtidx,win1idx,win2idx,:,:]
            fid = self.psi2fids[dtidx,win1idx,win2idx]
            u_fid = self.u_psi2fids[dtidx,win1idx,win2idx]

        
        ind = np.arange(4)
        w = 0.8
        
        fig,axs = plt.subplots(3,2,figsize=(8,12))
        for i,base in enumerate(bases):
            
            maxheight = 0
            rects = axs[i,0].bar(ind,raw[i,:], w, color='w', hatch='///')
            axs[i,0].set_xticks(ind+w/2)
            axs[i,0].set_xticklabels(['00', '01', '10', '11'])
            axs[i,0].set_xlabel('State')
            axs[i,0].set_xlim(-0.1, 3.9)
            axs[i,0].set_ylabel('Occurrences')
            axs[i,0].set_title(state+', '+base+', raw')
            for j,r in enumerate(rects):
                h = r.get_height()
                if h > maxheight: 
                    maxheight = h
                axs[i,0].text(ind[j]+w/2, h+1, str(int(h)), ha='center', va='bottom')
            
            axs[i,0].set_ylim(0,maxheight+10)

            rects = axs[i,1].bar(ind, corr[i,:], w, color='w', hatch='/',
                    yerr=u_corr[i,:], ecolor='k')
            axs[i,1].set_xticks(ind+w/2)
            axs[i,1].set_xticklabels(['00', '01', '10', '11'])
            axs[i,1].set_xlabel('State')
            axs[i,1].set_xlim(-0.1, 3.9)
            axs[i,1].set_ylabel('Fraction')
            axs[i,1].set_ylim(-0.1,0.75) 
            axs[i,1].hlines([0.], -0.1, 3.9, color='k')
            axs[i,1].set_title(state+', '+base+', corrected')
            for j,r in enumerate(rects):
                h = r.get_height()
                axs[i,1].text(ind[j]+w/2, h+0.05, '%.2f' % h, ha='center', va='bottom')

            if i == 0:
                axs[i,1].text(3.8, 0.7, 'F = %.3f $\pm$ %.3f' % (fid, u_fid),
                        ha='right', va='center')

        plt.tight_layout()

        if save:
            suffix = '' if self.mode == None else '_'+self.mode
            suffix += '' if self.dtslices == False else '_dtslices'

            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)
            fig.savefig(os.path.join(self.savedir, 'correlations%s_%s_win1%d_win2%d_dt%d.png' % \
                    (suffix, state, win1, win2, dt)))

if __name__ == '__main__':
    
    #### get all fidelities
    fid = FidelityAnalysis('Fidelity')
    
    fid.mode = None
    fid.dtslices = False
    
    fid.get_fidelities()
    fid.save_fidelities()
    
    # fid.load_fidelities('20121128-ldefidelity')
    # fid.plot_map_starts()
    # fid.plot_map_window(ch0start=641,ch1start=670)
    #fid.plot_correlations('psi1', 50, 150, 75)
    #fid.plot_correlations('psi2', 100, 150, 150)


    #### use this way to extract (and filter) entanglement events from the hhp-data
    # e = EntanglementEventAnalysis('X-X')
    # e.get_entanglement_events('X-X')
    # e.save_events('X-X')

    # e = EntanglementEventAnalysis('ZZ')
    # e.get_entanglement_events('ZZ')
    # e.save_events('ZZ')

    # e = EntanglementEventAnalysis('XX')
    # e.get_entanglement_events('XX')
    # e.save_events('XX')

    #### filtering and correlations
    #e.load_events('20121029-lde', 'X-X_highdc.npz')
    #e.filter_gate_phase()
    #e.filter_times()
    #e.correlations(save=False, F1_2=0.959, uF1_2=0.004)
