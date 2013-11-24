import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt
from analysis import config

from analysis.lib.ssro import readout_correction as roc

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'
config.datadir = r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/ssro-stats_withbasis/ssro_all'

rcParams['mathtext.default'] = 'regular'
rcParams['legend.numpoints'] = 1
rcParams['legend.frameon'] = False

class SSROStatAnalysis:

    def __init__(self):
        self.srcfolder = config.datadir
        self.lde_ssro_durations_fp = os.path.join(config.outputdir, 
                '20121105-ssro-stats', 'ssro_durations.npz')
        self.savedir = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-ssro-stats')
        self.fname0 = 'ADwin_SSRO-000_fid_vs_ROtime.npz'
        self.fname1 = 'ADwin_SSRO-001_fid_vs_ROtime.npz'
        self.lt1suffix = 'LT1'
        self.lt2suffix = 'LT2'

        # initialization (in)fidelities made in the calibration measurements
        self.initialization = {
                'LT1' : {
                    'p0' : 0.9951,
                    'u_p0' : 0.0001,
                    'p1' : 0.9970,
                    'u_p1' : 0.0002,
                    },
                'LT2' : {
                    'p0' : 0.9828,
                    'u_p0' : 0.0006,
                    'p1' : 0.9964,
                    'u_p1' : 0.0003,
                    },
                }

        self.fid = {
                'LT1' : {
                    0 : [],
                    1 : [],
                    },
                'LT2' : {
                    0 : [],
                    1 : [],
                    },
                }

        self.fiderr = {
                'LT1' : {
                    0 : [],
                    1 : [],
                    },
                'LT2' : {
                    0 : [],
                    1 : [],
                    },
                }

        self.pathtime = {
                'LT1' : [],
                'LT2' : [],
                }

        self.lt1idx = 14
        self.lt2idx = 25

    def load_data(self):
        
        print 'start reading data'

        for (path,dirs,files) in os.walk(self.srcfolder):           
            if self.fname0 in files and self.fname1 in files:
                
                # get information from the folder (time!), to be able to use
                # the correct readout time
                b,folder = os.path.split(path)
                timestamp = folder[:6]
                b,_date = os.path.split(b)
                date = _date[:8]

                if path[-3:] == self.lt1suffix:
                    setup = 'LT1'
                else:
                    setup = 'LT2'
               
                datetimestr = date+timestamp
                t = int(time.mktime(time.strptime(datetimestr, '%Y%m%d%H%M%S')))
                self.pathtime[setup].append(t)

                lt1duration, lt2duration = self.get_readoutduration(t)
                idx = lt1duration-1 if path[-3:] == self.lt1suffix else lt2duration-1
                
                # get the probabilities for 0/gt0 photons
                d0 = np.load(os.path.join(path,self.fname0))
                fid0 = d0['fid']
                fiderr0 = d0['fid_error']
                d0.close()

                _f0 = fid0[idx]
                _u_f0 = fiderr0[idx]

                d1 = np.load(os.path.join(path,self.fname1))
                fid1 = d1['fid']
                fiderr1 = d1['fid_error']
                d1.close()

                _f1 = fid1[idx]
                _u_f1 = fiderr1[idx]
                
                
                # calculate the real fidelities
                i = self.initialization[setup]
                F0,F1,u_F0,u_F1 = roc.readout_fidelity(i['p0'], i['p1'],
                        _f0, _f1, i['u_p0'], i['u_p1'], _u_f0, _u_f1)

                self.fid[setup][0].append(F0)
                self.fid[setup][1].append(F1)
                self.fiderr[setup][0].append(u_F0)
                self.fiderr[setup][1].append(u_F1)
    
    def get_readoutduration(self, ssrotime):
        f = np.load(self.lde_ssro_durations_fp)
        t = f['times']
        lt1 = f['lt1durations']
        lt2 = f['lt2durations']
        f.close()
        
        _t = (t-ssrotime) > 0
        closestidx = argmin(abs(t[_t]-ssrotime))

        print time.strftime('%m%d / %H%M', time.localtime(ssrotime)), '--->', 
        print time.strftime('%m%d / %H%M', time.localtime(t[closestidx]))
        
        return int(lt1[closestidx]), int(lt2[closestidx])
    
    def stats(self):
        
        self.meanfid = {
                'LT1' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'LT2' : {
                    0 : 0.,
                    1 : 0.,
                    },
                }

        self.fidstdev = {
                'LT1' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'LT2' : {
                    0 : 0.,
                    1 : 0.,
                    },
                }

        for setup in self.fid:
            for state in self.fid[setup]:
                print setup, state
                print 'mean = %.4f' % np.mean(np.array(self.fid[setup][state]))
                print 'statistical uncertainty = %.4f' % \
                        np.std(np.array(self.fid[setup][state]))
                print ''

                self.meanfid[setup][state] = np.mean(np.array(self.fid[setup][state]))
                self.fidstdev[setup][state] = np.std(np.array(self.fid[setup][state]))
                

        fig = plt.figure(figsize=(6,6))
        
        for i, setup in zip(range(1,3), ['LT1', 'LT2']):

            ax = plt.subplot(2,1,i)
            ax.set_title(setup)
            ax.hist(self.fid[setup][0], color='k', histtype='step', hatch='/',
                label='ms = 0', bins=5)
            plt.text(self.meanfid[setup][0],4.,'%.4f +/- %.4f' % \
                (self.meanfid[setup][0], self.fidstdev[setup][0]),
                color='k', backgroundcolor='w',
                ha='center')

            ax.hist(self.fid[setup][1], color='r', histtype='step', hatch='/',
                label='ms = 1', bins=5)
            plt.text(self.meanfid[setup][1],6.,'%.4f +/- %.4f' % \
                (self.meanfid[setup][1], self.fidstdev[setup][1]),
                color='r', backgroundcolor='w',
                ha='right')
            
            ax.legend(loc=2)
            ax.locator_params(nbins=4)
            ax.set_xlabel('Fidelity')
            ax.set_ylabel('Occurrences')
            ax.set_xlim(0.78,1)

        plt.tight_layout()

        # fig = plt.figure(figsize=(10,10))

        # for i, setup, state in zip(range(1,5), ['LT1', 'LT1', 'LT2', 'LT2'], [0,1,0,1]):

        #     ax = plt.subplot(2,2,i)
        #     ax.errorbar(np.arange(len(self.fid[setup][state])), np.array(self.fid[setup][state]),
        #             yerr=np.array(self.fiderr[setup][state]), fmt='ko')
        #     ax.axhline(y=self.meanfid[setup][state], ls=':')
        #     ax.set_title(setup + ', ms = ' + str(state))
        #     ax.set_xticks(np.arange(len(self.fid[setup][state])))
        #     ax.set_xticklabels(self.pathinfo[setup])
        #     labels = ax.get_xticklabels()
        #     plt.setp(labels, rotation=90)

        #     ax.set_ylabel('Fidelity')
        #     ax.set_xlabel('Run')

        # plt.tight_layout()


if __name__ == '__main__':
    a = SSROStatAnalysis()
    a.load_data()
    a.stats()

    if not os.path.exists(a.savedir):
        os.makedirs(a.savedir)

    f = open(os.path.join(a.savedir, 'mean_fidelities.pkl'), 'wb')
    pickle.dump(a.meanfid, f)
    f.close()

    f = open(os.path.join(a.savedir, 'mean_stdevs.pkl'), 'wb')
    pickle.dump(a.fidstdev, f)
    f.close()

    f = open(os.path.join(a.savedir, 'fidelities.pkl'), 'wb')
    pickle.dump(a.fid, f)
    f.close()

    f = open(os.path.join(a.savedir, 'stdevs.pkl'), 'wb')
    pickle.dump(a.fiderr, f)
    f.close()




