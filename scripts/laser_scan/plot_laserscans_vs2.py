import os
import pylab as plt
import numpy as np
import matplotlib
from analysis.lib.fitting import fit, common
from analysis.lib.tools import toolbox, plot
from analysis.lib.spin import spin_control as sc

from analysis.lib.nv import nvlevels

class LaserScanAnalysis:

    def init(self, datadir):
        self.gate_voltage = np.zeros(len(self.folders(datadir)))
        self.strain = np.zeros(len(self.folders(datadir)))
        self.offset = np.zeros(len(self.folders(datadir)))
        #self.EyExGate = np.zeros([len(self.folders(datadir)),3])
        #self.Ey_line = np.zeros(len(self.folders(datadir)))
        #self.Ex_line = np.zeros(len(self.folders(datadir)))
        
    def folders(self, datadir):    
        allcontent = os.listdir(datadir) 

        folders = []
        for i in allcontent:
            if 'LaserFrequencyScan' in i:
                folders.append(i)  

        if len(folders) == 0:
            logging.warning('No data found in datadir')

        return folders

    def result(self, folder):
        data = toolbox.measurement_filename(directory = folder, ext = 'npz')
        result = np.load(data)

        return result


    def plot_figure(self, folder, j):
        result = self.result(folder)
        r=result['data']

        plt.figure(j+1)
        plt.plot(r[:,1],r[:,2]/float(max(r[:,2])))    

        #plt.draw()
        plt.suptitle (str(os.path.split(folder)[1]), fontsize = 6)
        plt.xlabel ('Relative Frequency [GHz]')   
        plt.ylabel ('Counts (Normalized)')  
        
        plt.show()

        result.close()

    def fill_in_EyEx(self, folder, j):
        r = self.result(folder)['data']
        #self.gate_voltage[j] = r[0,3] #this is hte gate voltage as saved
        try:
            self.Ey_line=float(raw_input('Ey line?')) #GHz
            self.Ex_line=float(raw_input('Ex line?')) #GHz
            if self.plot_sils == False:
                self.gate_voltage[j] = float(raw_input('Gate voltage?'))
        except ValueError:
            print 'Could not understand input for lines'
            self.Ey_line = 0
            self.Ex_line = 0
            pass

        #self.EyExGate[j,:] = [self.Ey_line, self.Ex_line, self.gate_voltage]
        self.strain[j]=abs(self.Ex_line-self.Ey_line)/2.0
        self.offset[j]=np.min([self.Ey_line,self.Ex_line])+self.strain[j]
        print self.strain
        print self.offset

    def plot_all_scans(self, datadir, folder, j):
        r = self.result(folder)['data']

        plt.figure(len(self.folders(datadir))+1)    

        plt.plot(r[:,1],r[:,2]/float(max(r[:,2]))+j)
        plt.text(self.Ey_line-1,j+0.5,'Ey',fontsize=8)
        plt.text(self.Ex_line+0.3,j+0.5,'Ex',fontsize=8)


        if self.plot_sils:
            SIL_nr = float(raw_input('SIL number?'))  
            name='SIL'+str(SIL_nr)
        else:
            name = 'Gate voltage '+str(self.gate_voltage[j])
        plt.text(38.5,j+0.1,name,fontsize=5)

        plt.xlim([38,80])
        plt.ylim([-0.1,j+1])
        plt.xlabel ('Relative Frequency [GHz]')   
        plt.ylabel ('Counts (Normalized)')  
        


        self.result(folder).close()


    def plot_all_lines(self, datadir, folder):
        r = self.result(folder)['data']

        plt.savefig(os.path.join(os.path.split(folder)[0], 'all_scans.png'),
                format='png')

        plt.figure(len(self.folders(datadir))+2)    
        plt.xlabel('gate voltage (V)')
        plt.ylabel('(GHz)')

        strain=plt.plot(np.array(self.gate_voltage),np.array(self.strain), label='strain')
        offset=plt.plot(np.array(self.gate_voltage),np.array(self.offset), label='offset')
        plt.legend((strain, offset),('strain','offset'))

        plt.savefig(os.path.join(os.path.split(folder)[0], 'strain_offset_vs_gate_voltage.png'),
                format='png')

        self.result(folder).close()

  

def plot():
    """
    plots the laser scans in the given directory in different figures
    """
    a = LaserScanAnalysis()

    a.plot_sils = False
    datadir = r'\\Tudelft.net\staff-groups\tnw\ns\qt\Diamond\Documents\Documents Suzanne\LaserScans2'#r'd:\measuring\data\201307\20130729\LaserFrequencyScans' #folder that has all the data

    plt.close('all')
    folders = a.folders(datadir)
   
    for i,j in enumerate(folders):
        folder = os.path.join(datadir,j)
        #result = a.result(folder)
        print 'plot {} of {}'.format(i+1, len(folders))
        a.plot_figure(folder, i)

def find_EyEx():
    """
    asks for the Ey and Ex values of the through 'plot' plotted laserscans.
    With this is plots all laserscans in one figure.
    """
    a = LaserScanAnalysis()

    a.plot_sils = False
    datadir = r'\\Tudelft.net\staff-groups\tnw\ns\qt\Diamond\Documents\Documents Suzanne\LaserScans2'#r'd:\measuring\data\201307\20130729\LaserFrequencyScans' #folder that has all the data

    a.init(datadir)
    folders = a.folders(datadir)
   
    for i,j in enumerate(folders):
        folder = os.path.join(datadir,j)
        print 'Fill in the Ey and Ex lines of figure number {}'.format(i+1)
        a.fill_in_EyEx(folder, i)
        a.plot_all_scans(datadir, folder, i)
        a.plot_all_lines(datadir,folder)






