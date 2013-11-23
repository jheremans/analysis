### imports
import sys, os, time
import logging
import numpy as np
import h5py
from matplotlib import pyplot as plt
from analysis.lib.tools import toolbox


class M2Analysis:

    plot_format = 'png'

    def __init__(self, folder):
        self.folder = folder
        self.load_hdf5data()
        

    def load_hdf5data(self):
        self.h5filepath = toolbox.measurement_filename(self.folder)
        self.f = h5py.File(self.h5filepath, 'r')
        self.name = self.f.keys()[0]
        self.g = self.f[self.name]

        self.measurementstring = os.path.split(self.folder)[1]
        self.timestamp = os.path.split(os.path.split(self.folder)[0])[1] \
                + '/' + self.measurementstring[:6]
        self.measurementstring = self.measurementstring[7:]        
        self.default_plot_title = self.timestamp+'\n'+self.measurementstring

    def finish(self):
        self.f.close()

    def adwingrp(self, name=''):
        if name != '':
            adwingrpname = name
        else:
            if len(self.g.keys()) == 1:
                adwingrpname = self.g.keys()[0]
            else:
                logging.error("More than one measurement. Please give a name")
                return False

        return self.g[adwingrpname]

    def analysis_h5data(self, name='analysis'):
        if not os.path.exists(os.path.join(self.folder, name+'.hdf5')):
            mode = 'w'    
        else:
            mode = 'r+'            
            
        return h5py.File(os.path.join(self.folder, name+'.hdf5'), mode)

    def default_fig(self, **kw):
        figsize = kw.pop('figsize', (4,4))

        return plt.figure(figsize=figsize, **kw)

    def default_ax(self, fig=None, *arg, **kw):
        if fig == None:
            fig = self.default_fig(*arg, **kw)

        ax = fig.add_subplot(111)
        ax.set_title(self.timestamp+'\n'+self.measurementstring)

        return ax




