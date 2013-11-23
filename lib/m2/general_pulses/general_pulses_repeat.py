import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
from analysis.lib.m2 import m2
#from measurement.scripts.lt2_scripts.m2.general_pulses import AdwinPulses as ap

class GeneralPulsesRepeatAnlysis(m2.M2Analysis):

    def plot(self, save_plots=True):
        
        self.grp=self.g['pulses']
        self.max_element=self.grp.attrs['max_element']
        
        self.cycle_duration=self.grp.attrs['cycle_duration']/300*1e-6 #s
        
                
        self.durations=np.array(self.grp.attrs['element_durations'])*self.cycle_duration
        self.repetitions=self.grp['repetition_counter'].value
        self.raw_result=self.grp['results'].value
        self.raw_histogram=self.grp['histogram'].value
        self.timed_element_length=len(self.raw_histogram)/self.repetitions
        self.histogram=np.reshape(self.raw_histogram,(self.repetitions,self.timed_element_length))
        self.histogram=self.histogram/float(self.repetitions)/self.cycle_duration
        
        self.result=np.reshape(self.raw_result,(self.repetitions,self.max_element))
        self.first_count=self.grp['first_count'].value
        
        plot result
        plt.figure()
        plt.title(self.default_plot_title)
        #for i,x in enumerate(np.arange(self.max_element)):
        #    sp=plt.subplot(self.max_element,1,i+1)
        #    y=self.result[i]/(self.durations[i]*float(self.repetitions))
        #   plt.plot(self.sweep_axis,y)
        #    plt.xlabel(self.sweep_str)
        #    plt.ylabel('counts [Hz]')
        #    plt.text(self.max_sweep,0,'element'+str(x+1))
        
        #if save_plots:
        #    plt.savefig(os.path.join(self.folder,'result.png'))
        #    plt.close('all')
            
        #plot histograms
        #plt.figure()
        
        plt.plot(np.arange(len(self.histogram[0]))*self.cycle_duration,np.sum(self.histogram,axis=0), 'o')
        plt.xlabel('time [s]')
        plt.ylabel('counts [Hz]')
        if save_plots:
            plt.savefig(os.path.join(self.folder,'histogram.png'))
            plt.close('all')
        

def analyse_pulses(fp):
    a=GeneralPulsesRepeatAnlysis(fp)
    a.plot()
    return a