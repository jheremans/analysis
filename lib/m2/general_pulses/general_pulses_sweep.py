import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
from analysis.lib.m2 import m2
#from measurement.scripts.lt2_scripts.m2.general_pulses import AdwinPulses as ap

class GeneralPulsesSweepAnlysis(m2.M2Analysis):

    def plot(self, save_plots=True):
        
        self.grp=self.g['pulses']
        self.max_element=self.grp.attrs['max_element']
        
        self.max_sweep=self.grp.attrs['max_sweep']
        self.cycle_duration=self.grp.attrs['cycle_duration']/300*1e-6 #s
        
        if (self.grp.attrs['do_sweep_duration'] > 0):
            self.sweep_axis=np.array(self.grp.attrs['sweep_durations'])*self.cycle_duration*1e6
            self.sweep_str='sweep duration [us]'
        else:
            self.raw_sweep_axis=np.array(self.grp.attrs['sweep_voltages'])
            self.sweep_axis=np.array(self.g.attrs['sweep_axis'])*1e9
            self.raw_sweep_str='dac_voltage [V]'
            self.sweep_str='yellow power [nW]'
        
        self.durations=np.array(self.grp.attrs['element_durations'])*self.cycle_duration
        self.repetitions=self.grp['repetition_counter'].value
        raw_result=self.grp['results'].value
        raw_histogram=self.grp['histogram'].value
        raw_counts=self.grp['counter'].value/float(self.repetitions)
        
        
        self.result=np.split(raw_result,self.max_element)
        
        self.histograms=np.split(raw_histogram,self.max_sweep)# (or len(raw_results)/max_sweep)
        self.max_counts=len(self.histograms[0])
        #plot result
        plt.figure()
        plt.title(self.default_plot_title)
        for i,x in enumerate(np.arange(self.max_element)):
            sp=plt.subplot(self.max_element,1,i+1)
            y=self.result[i]/(self.durations[i]*float(self.repetitions))
            plt.plot(self.sweep_axis,y)
            plt.xlabel(self.sweep_str)
            plt.ylabel('counts [Hz]')
            plt.text(self.max_sweep,0,'element'+str(x+1))
        
        if save_plots:
            plt.savefig(os.path.join(self.folder,'result.png'))
            plt.close('all')
            
        #plot histograms
        plt.figure()
        off=0
        offstep=1.2*np.max(self.histograms)
        for i,x in enumerate(self.sweep_axis):
            plt.bar(np.arange(self.max_counts),self.histograms[i], bottom=off)
            plt.text(0,off,self.sweep_str+':'+str(x))
            off=off+offstep
        if save_plots:
            plt.savefig(os.path.join(self.folder,'histogram.png'))
            plt.close('all')
        

def analyse_pulses(fp):
    a=GeneralPulsesSweepAnlysis(fp)
    a.plot()
    return a