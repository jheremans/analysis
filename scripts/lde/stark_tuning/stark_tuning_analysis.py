import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm



class StarkAnalysis:

    def __init__(self):

        self.srcfolder = r'D:\measuring\data\LDE\analysis_data\stark_tuning'
        self.scanidf = r'LT25nW_mw_True_0nW_green'
        self.save_folder = r'D:\measuring\data\LDE\analysis_output\20121019-stark-tuning'
        self.common_freq_axis = np.linspace(40,90,5001)

    def find_nearest(self, array,value):
        idx=(abs(array-value)).argmin()
        return idx

    def get_data(self,file_name):
        d=np.loadtxt(os.path.join(self.srcfolder, file_name))
        freq=d[:,1]
        counts=d[:,2]
        return freq, counts

    def change_freq_axis(self, freq_axis, counts):
        self.new_counts=np.ones(len(self.common_freq_axis))*-1

        for i in np.arange(len(counts)):
            idx=self.find_nearest(self.common_freq_axis,freq_axis[i])
            try:
                if self.new_counts[idx]<0:
                    self.new_counts[idx]=counts[i]
                else:
                    self.new_counts[idx]+=counts[i]
            except:
                print 'something went wrong with index ', idx


        return self.new_counts

    def get_files(self):
        self.files=[]
        for d,sd,f in os.walk(self.srcfolder):
            self.files=np.hstack((self.files,
                [os.path.abspath(os.path.join(d,i)) \
		    for i in f if self.scanidf in i \
			and i[:1] != '.' and i[-5:]== '1.dat']))

        return self.files

    def mk_scan_array(self):
	self.scan_array=np.zeros([0,len(self.common_freq_axis)])
	for file_name in self.files:
	    freq, counts = self.get_data(file_name)
	    self.scan_array = np.vstack((self.scan_array, self.change_freq_axis(freq,counts)))
	return self.scan_array

    def plot_all_scans(self,scan_array,freq_min=51.,freq_max=72.):
	idx_min = self.find_nearest(self.common_freq_axis,freq_min)
	idx_max = self.find_nearest(self.common_freq_axis,freq_max)
	self.scans_positive = [scan_array[int(i*2),idx_min:idx_max] for i in range(int(len(scan_array)/2.))]
	print len(self.scans_positive)
	self.scans_negative = [scan_array[int(i*2+1),idx_min:idx_max] for i in range(int((len(scan_array))/2.))]
	print len(self.scans_negative)
	fig = figure()
	subplot(211)
	imgplot = plt.imshow(self.scans_positive, cmap='bone_r', 
	    aspect=(idx_max-idx_min)/100./len(self.scans_positive), extent=[freq_min,freq_max,18,1])#, norm=LogNorm())
	plt.xlabel('frequency [GHz]')
	plt.ylabel('positive voltages [V]')
	imgplot.set_clim(1.,50.)
	subplot(212)
	imgplot = plt.imshow(self.scans_negative, cmap='bone_r', 
	    aspect=(idx_max-idx_min)/100./len(self.scans_negative), extent=[freq_min,freq_max,-18,-1])
	plt.xlabel('frequency [GHz]')
	plt.ylabel('negative voltages [V]')
#	plt.colorbar()
	imgplot.set_clim(1.,50.)

	    

if __name__ == '__main__':
    first_time = False
    a = StarkAnalysis()
    if first_time:
        a.get_files()
	scan_array = a.mk_scan_array()
    a.plot_all_scans(scan_array=scan_array)

    np.savez(os.path.join(a.save_folder,'stark_tuning.npz'),
            all_scans=scan_array,
            positiveV=a.scans_positive, 
            negativeV=a.scans_negative,frequency=a.common_freq_axis)


