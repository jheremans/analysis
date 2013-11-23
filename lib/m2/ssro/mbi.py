import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.spin import Nspin_twolevel_correction
from analysis.lib.spin import Nspin_correction
from analysis.lib.spin import N_and_e_spin_correction
reload(N_and_e_spin_correction)
from analysis.lib.m2 import m2
from analysis.lib.m2.ssro import ssro
from analysis.lib.math import error
from analysis.lib.tools import toolbox

class MBIAnalysis(m2.M2Analysis):
    def get_readout_results(self, name=''):
        """
        Get the readout results.
        self.ssro_results contains the readout results (sum of the photons for 
            a given sweep point and readout in a sequence)
        self.normalized_ssro contains the normalized result (i.e., probability
            for getting a photon)
        """
        self.result_corrected = False

        adwingrp = self.adwingrp(name)
        self.adgrp = adwingrp

        self.reps = adwingrp.attrs['reps_per_ROsequence']
        self.pts = adwingrp.attrs['sweep_length']
        self.readouts = adwingrp.attrs['nr_of_ROsequences']
        discards = len(adwingrp['ssro_results'].value) % (self.pts*self.readouts)
        self.reps = int(len(adwingrp['ssro_results'].value) / (self.pts*self.readouts))      
        
        if discards > 0:
            results = adwingrp['ssro_results'].value.reshape(-1)[:-discards]
        else:
            results = adwingrp['ssro_results'].value

        self.ssro_results = results.reshape((-1,self.pts,self.readouts)).sum(axis=0)
        self.normalized_ssro = self.ssro_results/float(self.reps)
        self.u_normalized_ssro = \
            (self.normalized_ssro*(1.-self.normalized_ssro)/self.reps)**0.5

    def get_correlations(self, name=''):
        """
        computes the readout correlations between all readouts in a sequence.
        we encode correlations in binary: 
            C = \sum_i{ 2^i \times r_i }
        where i is the i-th readout (i=0 is the first), and r is the SSRO result
        of that readout, with r = 1 for ms=1, and r = 0 for ms=0.
        (Yes, that's opposite to the normalized_ssro result, where 0 means
        always 0 photons and 1 means always > 0 photons!) 
        """
        adwingrp = self.adwingrp(name)
        self.result_correlation_corrected = False
        
        self.reps = self.g.attrs['reps_per_ROsequence']
        self.pts = adwingrp.attrs['sweep_length']
        self.readouts = adwingrp.attrs['nr_of_ROsequences']
        
        discards = len(adwingrp['ssro_results'].value) % (self.pts*self.readouts)
        self.reps = int(len(adwingrp['ssro_results'].value) / (self.pts*self.readouts))
        
        if discards > 0:
            results = adwingrp['ssro_results'].value.reshape(-1)[:-discards]
        else:
            results = adwingrp['ssro_results'].value

        _res = results.reshape((-1,self.pts,self.readouts))
        c = np.zeros((self.reps, self.pts))
        
        #Here we make sure that outcomes are ordered such that in the example of 2 RO's 
        #(RO1,RO2)=('e','N') -> c -> repr  is 
        #(1,1) -> 3 -> 11; (1,0) -> 1 -> 01 ; (0,1) -> 2 -> 10 , (0,0) ->0 -> 00 
        for i in range(self.readouts):
            c += 2**i * (1-_res[:,:,i])

        corrvals = range(2**self.readouts)
        self.correlations = np.zeros((self.pts, len(corrvals)))
        for v in corrvals:
            for p in range(self.pts):
                self.correlations[p,v] = len(np.where(c[:,p]==v)[0])
        
        #The correlations are ordered in the order Ne, in binary rep: 00,01,10,11 respectively.
        _cnames = [np.binary_repr(
            v, width=self.readouts)[::-1] for v in corrvals]
        self.correlation_names = []
        for cn in _cnames:
            n = ''
            for j in range(self.readouts):
                if j > 0:
                    n += ' & '
                n += r'{}:$|{}\rangle$'.format(j+1, cn[j])

            self.correlation_names.append(n)
        
        self.normalized_correlations = self.correlations / float(self.reps)
        self.u_normalized_correlations = (
            self.normalized_correlations*(1.-self.normalized_correlations)/self.reps)**0.5

    def get_sweep_pts(self):
        self.sweep_name = self.g.attrs['sweep_name']
        self.sweep_pts = self.g.attrs['sweep_pts']

    def get_electron_ROC(self, ssro_calib_folder=''):
        if ssro_calib_folder == '':
            ssro_calib_folder = toolbox.latest_data('SSRO')

        self.p0 = np.zeros(self.normalized_ssro.shape)
        self.u_p0 = np.zeros(self.normalized_ssro.shape)
        
        ro_durations = self.g.attrs['E_RO_durations']
        roc = error.SingleQubitROC()

        for i in range(len(self.normalized_ssro[0])):
            roc.F0, roc.u_F0, roc.F1, roc.u_F1 = \
                ssro.get_SSRO_calibration(ssro_calib_folder, 
                        ro_durations[i])
            p0, u_p0 = roc.num_eval(self.normalized_ssro[:,i],
                    self.u_normalized_ssro[:,i])
            
            self.p0[:,i] = p0
            self.u_p0[:,i] = u_p0
        
        self.result_corrected = True

    def get_N_ROC(self, P_min1=1, u_P_min1=0, P_0=1, u_P_0=0,
            F0_RO_pulse=1, u_F0_RO_pulse=0,
            F1_RO_pulse=1, u_F1_RO_pulse=0,
            ssro_calib_folder=None):

        if ssro_calib_folder == None:
            ssro_calib_folder = toolbox.latest_data('SSRO')
        
        self.p0 = np.zeros(self.normalized_ssro.shape)
        self.u_p0 = np.zeros(self.normalized_ssro.shape)
        
        ro_durations = self.g.attrs['E_RO_durations']
        roc = Nspin_correction.NuclearSpinROC()
        roc.P_min1 = P_min1
        roc.u_P_min1 = u_P_min1
        roc.P_0 = P_0
        roc.u_P_0 = u_P_0
       
        roc.F0_RO_pulse = F0_RO_pulse
        roc.u_F0_RO_pulse = u_F0_RO_pulse
        roc.F1_RO_pulse = F1_RO_pulse
        roc.u_F1_RO_pulse = u_F1_RO_pulse
        
        for i in range(len(self.normalized_ssro[0])):
            roc.F0_ssro, roc.u_F0_ssro, roc.F1_ssro, roc.u_F1_ssro = \
                ssro.get_SSRO_calibration(ssro_calib_folder,
                        ro_durations[i])
        
            p0, u_p0 = roc.num_evaluation(self.normalized_ssro[:,i],
                    self.u_normalized_ssro[:,i])
            
            self.p0[:,i] = p0
            self.u_p0[:,i] = u_p0
        
        self.result_corrected = True

    def get_correlation_ROC(self, P_min1=1, u_P_min1=0, P_0=0, u_P_0=0,
            F0_RO_pulse=1, u_F0_RO_pulse=0,
            F1_RO_pulse=1, u_F1_RO_pulse=0,
            ssro_calib_folder=None):

        if ssro_calib_folder == None:
            ssro_calib_folder = toolbox.latest_data('SSRO')
        
        #The shape below is [sweep_pts,4].
        self.p_correlations= np.zeros(self.normalized_correlations.shape)
        self.u_p_correlations = np.zeros(self.normalized_correlations.shape)
        
        ro_durations = self.g.attrs['E_RO_durations']
        roc = N_and_e_spin_correction.CorrelationsROC()
        roc.P_min1 = P_min1
        roc.u_P_min1 = u_P_min1
        roc.P_0 = P_0
        roc.u_P_0 = u_P_0
       
        roc.F0_RO_pulse = F0_RO_pulse
        roc.u_F0_RO_pulse = u_F0_RO_pulse
        roc.F1_RO_pulse = F1_RO_pulse
        roc.u_F1_RO_pulse = u_F1_RO_pulse
        
        #The electron RO correction for the first readout (with ro_durations[0])
        roc.F0_e_ssro, roc.u_F0_e_ssro, roc.F1_e_ssro, roc.u_F1_e_ssro = \
            ssro.get_SSRO_calibration(ssro_calib_folder,
                       ro_durations[0])
        
        #The electron RO correction for the second readout (with ro_durations[1])  
        roc.F0_N_ssro, roc.u_F0_N_ssro, roc.F1_N_ssro, roc.u_F1_N_ssro = \
            ssro.get_SSRO_calibration(ssro_calib_folder,
                   ro_durations[1])        
        
        #Here the correction for the double readout is calculation. 
        p_correlations, u_p_correlations = roc.num_evaluation(self.normalized_correlations[:,:],
                self.u_normalized_correlations[:,:])
        
        #The output p_correlations above is tuple, we transform it to array 
        #and Transpose it to retain the right shape.
        self.p_correlations = np.array(p_correlations).T
        self.u_p_correlations = np.array(u_p_correlations).T
        
        self.result_correlation_corrected = True

    def save(self,correction, fname='analysis.hdf5'):
        f = h5py.File(os.path.join(self.folder, fname), 'w')
        g = f.create_group('readout_result')
        
        g.attrs['sweep_name'] = self.sweep_name

        f['/readout_result/sweep_pts'] = self.sweep_pts
        f['/readout_result/normalized_ssro'] = self.normalized_ssro
        f['/readout_result/u_normalized_ssro'] = self.u_normalized_ssro
        
        if correction != 'correlation':
            if self.result_corrected:
                f['readout_result/p0'] = self.p0
                f['readout_result/u_p0'] = self.u_p0
        else:
            f['/readout_result/normalized_correlations'] = self.normalized_correlations
            f['/readout_result/u_normalized_correlations'] = self.u_normalized_correlations
            if self.result_correlation_corrected:
                f['readout_result/p_correlations'] = self.p_correlations
                f['readout_result/u_p_correlations'] = self.u_p_correlations
                
        f.close()
        

    def plot_results_vs_sweepparam(self, name='', save=True, **kw):
        mode = kw.get('mode', 'ssro_results')
        labels = kw.get('labels', None)
        ret = kw.get('ret', None)
        ylim = kw.get('ylim', (-0.05, 1.05))
        ax = kw.get('ax', None)

        if not hasattr(self, 'sweep_pts'):
            self.sweep_pts = np.arange(self.pts) + 1
            self.sweep_name = 'sweep parameter'
        
        if ax == None:
            fig = self.default_fig(figsize=(6,4))
            ax = self.default_ax(fig)
        else:
            ax.set_title(self.timestamp+'\n'+self.measurementstring)
        
        if labels == None:
            if mode != 'correlations':
                labels = ['RO #%d' % i for i in range(self.readouts)]
            else:
                labels = self.correlation_names
        
        if mode != 'correlations':
            for i in range(self.readouts):
                if not self.result_corrected:
                    ax.errorbar(self.sweep_pts, self.normalized_ssro[:,i], fmt='o-',
                        yerr=self.u_normalized_ssro[:,i], label=labels[i])
                else:
                    ax.errorbar(self.sweep_pts, self.p0[:,i], fmt='o',
                        yerr=self.u_p0[:,i], label=labels[i])
                    ax.axhspan(0,1,fill=False,ls='dotted')
        
            if not self.result_corrected:
                ax.set_ylabel('avg (uncorrected) outcome')
            else:
                ax.set_ylabel(r'$F(|0\rangle)$')

        else:
            for i in range(len(self.correlation_names)):
                if not self.result_correlation_corrected:
                    ax.errorbar(self.sweep_pts, self.normalized_correlations[:,i], 
                        fmt='o', yerr=self.u_normalized_correlations[:,i], 
                        label=labels[i])
                else:
                    ax.errorbar(self.sweep_pts, self.p_correlations[:,i], fmt='o',
                        yerr=self.u_p_correlations[:,i], label=labels[i])
                    ax.axhspan(0,1,fill=False,ls='dotted')
            
            if not self.result_correlation_corrected:
                ax.set_ylabel('avg (uncorrected) outcome')
            else:
                ax.set_ylabel('Probability')
            
        
        ax.set_xlabel(self.sweep_name)       
        ax.set_ylim(ylim)
        
        if self.readouts > 1:
            ax.legend(loc='best')

        if save and ax != None:
            try:
                fig.savefig(
                    os.path.join(self.folder, '{}_vs_sweepparam.'.format(mode) + self.plot_format),
                    format=self.plot_format)
            except:
                print 'Figure has not been saved.'

        if ret == 'ax':
            return ax
        if ret == 'fig':
            return fig

    
    def finish(self):
        self.f.close()

def analyze_single_sweep(folder, name='', correction='electron', **kw):
    P_min1 = kw.pop('P_min1',0.95)
    u_P_min1 = kw.pop('u_P_min1',0.03)
    P_0 = kw.pop('P_0',0.04)
    u_P_0 = kw.pop('u_P_0',0.03)
    
    F0_RO_pulse = kw.pop('F0_RO_pulse', 0.98)
    u_F0_RO_pulse = kw.pop('u_F0_RO_pulse', 0.0)
    F1_RO_pulse = kw.pop('F1_RO_pulse', 0.99)
    u_F1_RO_pulse = kw.pop('u_F1_RO_pulse', 0.0)
    
    ret = kw.pop('ret', None)
    mode = kw.pop('mode', 'ssro_results')

    a = MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name)

    
    if mode == 'correlations':
        a.get_correlations(name)

    if correction == 'electron':
        a.get_electron_ROC()
    elif correction == 'N':
        a.get_N_ROC(P_min1, u_P_min1, P_0, u_P_0, F0_RO_pulse, u_F0_RO_pulse, F1_RO_pulse,
                u_F1_RO_pulse)
    elif correction == 'correlation':
        a.get_correlation_ROC(P_min1, u_P_min1, P_0, u_P_0, F0_RO_pulse, u_F0_RO_pulse, F1_RO_pulse,
                u_F1_RO_pulse)

    a.save(correction)    
    a.plot_results_vs_sweepparam(mode = mode, **kw)
    a.finish()
    
    if ret == 'obj':
        return a
        
    
