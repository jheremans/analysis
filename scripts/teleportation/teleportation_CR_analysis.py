import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import ssro
from analysis.lib.math import error
from analysis.lib.m2 import m2
from measurement.lib.tools import toolbox

class TeleportationAnalysis(m2.M2Analysis):
    def load_hdf5data(self):
        self.h5filepath = toolbox.measurement_filename(self.folder)
        self.f = h5py.File(self.h5filepath, 'r')
        self.name = self.f.keys()[5] #nr 5 has adwin data
        self.g = self.f[self.name]

        self.measurementstring = os.path.split(self.folder)[1]
        self.timestamp = os.path.split(os.path.split(self.folder)[0])[1] \
                + '/' + self.measurementstring[:6]
        self.measurementstring = self.measurementstring[7:]        
        self.default_plot_title = self.timestamp+'\n'+self.measurementstring

    def get_cr_results(self, name='', plot=True):
        adwingrp = self.adwingrp(name)   
        self.cr_hist_all = adwingrp['CR_hist_all'].value
        self.cr_hist_all = [float(i) for i in self.cr_hist_all]
        self.cr_prob_hist_all = self.cr_hist_all/ np.sum( self.cr_hist_all)

        if plot:
            self.plot_cr_results('cr hist all', self.cr_prob_hist_all)


    def plot_cr_results(self,name,cr_counts):   
        title_suffix = ': '+name if name != '' else ''
        fn_suffix = '_'+name if name != '' else ''

        fig = self.default_fig(figsize=(6,4))
        ax = self.default_ax(fig)

        ax.bar(cr_counts, 'o')
        ax.set_xlabel('counts during CR check sequence')
        ax.set_ylabel('number of counts')
        ax.set_title(self.default_plot_title + title_suffix)

def analyze_cr_teleportation(folder, name = ''):
        a = TeleportationAnalysis(folder)
        a.get_cr_results(name, plot = True)
        a.finish()

        

