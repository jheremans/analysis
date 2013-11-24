# we're only concerned with LT2 here
#

import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt
from analysis import config

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'
config.datadir = r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/ssro-stats_withbasis/ssro_all'

class SSROStatAnalysis:

    def __init__(self):
        self.srcfolder = config.datadir
        self.savedir = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-ssro-stats_bybasis')
        self.fname0 = 'ADwin_SSRO-000_fid_vs_ROtime.npz'
        self.fname1 = 'ADwin_SSRO-001_fid_vs_ROtime.npz'
        self.lt1suffix = 'LT1'
        self.lt2suffix = 'LT2'

        self.fid = {
                'ZZ' : {
                    0 : [],
                    1 : [],
                    },
                'XX' : {
                    0 : [],
                    1 : [],
                    },
                'X-X' : {
                    0 : [],
                    1 : [],
                    },
                }

        self.lt2idx = 25

    def load_data(self):
        
        print 'start reading data'

        for (path,dirs,files) in os.walk(self.srcfolder):
            for fn in files:
                if fn == self.fname0 or fn == self.fname1:
                    d = np.load(os.path.join(path, fn))
                    fid = d['fid']
                    fiderr = d['fid_error']
                    d.close()

                    b,folder = os.path.split(path)
                    timestamp = folder[:4]
                    b,date = os.path.split(b)[-4:]
                    
                    if path[-3:] == self.lt1suffix:
                        setup = 'LT1'
                        continue
                    else:
                        setup = 'LT2'

                    if 'ZZ' in path:
                        basis = 'ZZ'
                    elif 'XX' in path:
                        basis = 'XX'
                    elif 'X-X' in path:
                        basis = 'X-X'
                    else:
                        print 'missed a set:', path

                    state = 0 if fn == self.fname0 else 1
                    idx = self.lt1idx if path[-3:] == self.lt1suffix else self.lt2idx

                    self.fid[basis][state].append(fid[idx])

    def stats(self):
        
        self.meanfid = {
                'ZZ' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'XX' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'X-X' : {
                    0 : 0.,
                    1 : 0.,
                    },
                }

        self.stdev = {
                'ZZ' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'XX' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'X-X' : {
                    0 : 0.,
                    1 : 0.,
                    },
                }

        for basis in self.fid:
            for state in self.fid[basis]:
                print basis, state
                print 'mean = %.4f' % np.mean(np.array(self.fid[basis][state]))
                print 'statistical uncertainty = %.4f' % \
                        np.std(np.array(self.fid[basis][state]))
                print ''

                self.meanfid[basis][state] = np.mean(np.array(self.fid[basis][state]))
                self.stdev[basis][state] = np.std(np.array(self.fid[basis][state]))
                

        fig,ax = plt.subplots(1,2, figsize=(8,4))
        
        for i in [0,1]:
            state = i
            ax[i].set_title('LT2, ms = ' + str(i))
            
            for basis, color in zip(['ZZ', 'XX', 'X-X'], ['k', 'r', 'b']):
                ax[i].hist(self.fid[basis][state], color=color, histtype='step', hatch='', label=basis,
                        lw=2)
            ax[i].legend()            
            ax[i].locator_params(nbins=4)
            ax[i].set_ylim(0,6)
            xpos = ax[i].get_xlim()[0] + (ax[i].get_xlim()[1]-ax[i].get_xlim()[0])/10.
            ypos = 5. # ax[i].get_ylim()[1] - (ax[i].get_ylim()[1]-ax[i].get_ylim()[0])/10.

            for j,basis,color in zip([0,1,2], ['ZZ', 'XX', 'X-X'], ['k', 'r', 'b']):
                ax[i].text(xpos,ypos,2*j*'\n'+'%.1f $\pm$ %.1f %%' % (self.meanfid[basis][state]*100.,
                    self.stdev[basis][state]*100.), color=color, va='top', ha='left')

            ax[i].set_xlabel('Fidelity')
            ax[i].set_ylabel('Occurrences')
            
        plt.tight_layout()

if __name__ == '__main__':
    a = SSROStatAnalysis()
    a.load_data()
    a.stats()

#     if not os.path.exists(a.savedir):
#         os.makedirs(a.savedir)
# 
#     f = open(os.path.join(a.savedir, 'mean_fidelities.pkl'), 'wb')
#     pickle.dump(a.meanfid, f)
#     f.close()
# 
#     f = open(os.path.join(a.savedir, 'mean_stdevs.pkl'), 'wb')
#     pickle.dump(a.fidstdev, f)
#     f.close()
# 
#     f = open(os.path.join(a.savedir, 'fidelities.pkl'), 'wb')
#     pickle.dump(a.fid, f)
#     f.close()
# 
#     f = open(os.path.join(a.savedir, 'stdevs.pkl'), 'wb')
#     pickle.dump(a.fiderr, f)
#     f.close()




