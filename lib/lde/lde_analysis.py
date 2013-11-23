<<<<<<< HEAD
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cbook
import pickle

from fitting import fit,common
from lde import spcorr, sscorr, tpqi
from pq import hht3

HHPFILEBASE='hhp_data'
CHMAXTIME=2300


class SingleLDEAnalysis:

    datfolder = ''
    rawfolder = ''
 
    def __init__(self):
        pass        

    def prepare_data(self, ch0maxtime=CHMAXTIME, ch1maxtime=CHMAXTIME, DATIDX = 0):
        self.hhdata, self.w1, self.w2 = sscorr.load_hhdata(self.rawfolder, ch0maxtime=ch0maxtime,
                ch1maxtime=ch1maxtime, do_filter_crap=False)
        self.hhp, self.hhp1, self.hhp2 = sscorr.get_hhpludata(self.hhdata)

       
        self.ssro1, self.ssro2, self.cr1, self.cr2, self.gate_phase = sscorr.load_adwin_data(
                self.datfolder, DATIDX=DATIDX)
        if (len(self.hhp)>0):
            self.noof_hh_ssros = len(self.hhp[np.logical_and(self.hhp[:,3]==1, 
                self.hhp[:,2]==2)])
        else:
            self.noof_hh_ssros=0
        self.noof_adwin_ssros = len(self.ssro2)

        if self.noof_hh_ssros > self.noof_adwin_ssros:
            print 'HH SSROs do not match Adwin!'
            raise(Exception('hharp vs ssro mismatch'))

        np.savez(os.path.join(self.datfolder,HHPFILEBASE+('-%.'+str(sscorr.DATIDXDIGITS)+'d') % DATIDX),\
                hhdata=self.hhdata, hhp=self.hhp, w1=self.w1, w2=self.w2)
        
    def get_adwin_data(self, DATIDX = 0):
        self.ssro1, self.ssro2, self.cr1, self.cr2, self.gate_phase = sscorr.load_adwin_data(
                self.datfolder, DATIDX=DATIDX)  
        return True
    
    def get_dict_from_npz(self,d):
        e={}
        for key in d.keys():
            e[key]=d[key]
        return e
        
    def get_adwin_statistics(self,DATIDX):
        fn = os.path.join(self.datfolder, 
                sscorr.ADWINLT1FOLDER+('-%.'+str(sscorr.DATIDXDIGITS)+'d') % DATIDX, 
                sscorr.ADWINDATAFN)
        d = np.load(fn)
        self.adwin_lt1_stats=self.get_dict_from_npz(d)
        d.close()
        
        sn = os.path.join(self.datfolder, 
                sscorr.ADWINLT2FOLDER+('-%.'+str(sscorr.DATIDXDIGITS)+'d') % DATIDX, 
                sscorr.ADWINDATAFN)
        s = np.load(sn)
        self.adwin_lt2_stats=self.get_dict_from_npz(s)
        s.close()


    def correlations(self, *arg, **kw):
        if (len(self.hhp)>0):
            self.w = sscorr.windowdata(self.hhp, w1_start=self.w1_start, 
                    w1_stop=self.w1_stop, w2_start=self.w2_start, 
                    w2_stop=self.w2_stop, dt_max=self.dt_max,dt_min=self.dt_min,ch1_offset=self.ch1_offset,
                    single_photon_range_start=self.single_photon_range_start,
                    single_photon_range_stop=self.single_photon_range_stop,
                    *arg, **kw)
            self.uncond_corr, self.corr_00, self.corr_01, self.corr_10, self.corr_11 =\
                    sscorr.correlations(self.ssro1, self.ssro2, self.w)
        else:
            self.uncond_corr=np.zeros(4)

            self.corr_00=np.zeros(4)
            self.corr_01=np.zeros(4)
            self.corr_10=np.zeros(4)
            self.corr_11=np.zeros(4)
            self.w=[]

        return True
    
    def tail_counts(self):
        """calculates the tail counts - needs get_adwin_statistics to be run first"""
        w1tailcts_ch0=len(np.where(
            np.logical_and(self.w1[:,2]==0,np.logical_and(self.w1[:,1]>self.w1_start[0],self.w1[:,1]<self.w1_stop[0])))[0])
        w1tailcts_ch1=len(np.where(
            np.logical_and(self.w1[:,2]==1,np.logical_and(self.w1[:,1]>self.w1_start[1],self.w1[:,1]<self.w1_stop[1])))[0])
        w2tailcts_ch0=len(np.where(
            np.logical_and(self.w2[:,2]==0,np.logical_and(self.w2[:,1]>self.w2_start[0],self.w2[:,1]<self.w2_stop[0])))[0])
        w2tailcts_ch1=len(np.where(
            np.logical_and(self.w2[:,2]==1,np.logical_and(self.w2[:,1]>self.w2_start[1],self.w2[:,1]<self.w2_stop[1])))[0])
        seq_starts=self.adwin_lt2_stats['get_noof_seq_starts']
        if float(seq_starts)>0:
            self.total_tail=1.*(w1tailcts_ch0+w1tailcts_ch1+w2tailcts_ch0+w2tailcts_ch1)/(float(seq_starts)*300.)
        else:
            self.total_tail=0

    def laser_counts(self, laser_length=16):
        """calculates the laser counts - needs get_adwin_statistics to be run first"""
        w1tailcts_ch0=len(np.where(
            np.logical_and(self.w1[:,2]==0,np.logical_and(self.w1[:,1]>self.w1_start[0]-laser_length,self.w1[:,1]<self.w1_start[0])))[0])
        w1tailcts_ch1=len(np.where(
            np.logical_and(self.w1[:,2]==1,np.logical_and(self.w1[:,1]>self.w1_start[1]-laser_length,self.w1[:,1]<self.w1_start[1])))[0])
        w2tailcts_ch0=len(np.where(
            np.logical_and(self.w2[:,2]==0,np.logical_and(self.w2[:,1]>self.w2_start[0]-laser_length,self.w2[:,1]<self.w2_start[0])))[0])
        w2tailcts_ch1=len(np.where(
            np.logical_and(self.w2[:,2]==1,np.logical_and(self.w2[:,1]>self.w2_start[1]-laser_length,self.w2[:,1]<self.w2_start[1])))[0])
        seq_starts=self.adwin_lt2_stats['get_noof_seq_starts']
        if float(seq_starts)>0:
            self.total_laser=1.*(w1tailcts_ch0+w1tailcts_ch1+w2tailcts_ch0+w2tailcts_ch1)/(float(seq_starts)*300.)
        else:
            self.total_laser=0
        

###NOTE does not consider gate phase yet
    def get_event_photons(self, windowdata, hhp):
        plu_idx=np.where(hhp[:,3]==1)[0][np.where(windowdata>0)]
        event_photons= np.zeros((0,4),dtype=np.uint32)
        
        for i in plu_idx:
            event_photons=np.vstack((event_photons,hhp[i-2]))
            event_photons=np.vstack((event_photons,hhp[i-1]))
            event_photons=np.vstack((event_photons,hhp[i]))
        return event_photons

         
    def get_g2_coincidences(self, deltas, **kw):
        # plotsavename = 'g2_ch0-%d-%d_ch1-%d-%d.png'
        data = kw.get('hhdata',self.hhdata)
        ch0window=kw.get('ch0window',(self.w1_start[0], self.w1_stop[0]))
        ch1window=kw.get('ch1window',(self.w1_start[1], self.w1_stop[1]))

        data = hht3.filter_timewindow(data, 0, 
                mintime=ch0window[0], maxtime=ch0window[1])
        data = hht3.filter_timewindow(data, 1, 
                mintime=ch1window[0], maxtime=ch1window[1])

        d,c = tpqi.coincidences(data, deltas)
        # tpqi.plot_g2(d,c,savepath=os.path.join(self.savedir,plotsavename))
        return c


    def print_SSRO(self, DATIDX=0):
        self.ssro1, self.ssro2, self.cr1, self.cr2 = \
                sscorr.load_adwin_data(self.datfolder, DATIDX=DATIDX)
        print self.ssro2

    def print_plu_vs_hh(self):
        if (len(self.hhp)>0):
            self.noof_hh_ssros= len(self.hhp[np.logical_and(self.hhp[:,3]==1, 
                    self.hhp[:,2]==2)])
        else:
            self.noof_hh_ssros=0
        if len(self.w)>0:    
            self.noof_valid_ssros=len(np.where(self.w>0)[0])
        else:
            self.noof_valid_ssros=0
        print 'Plu gave', self.noof_hh_ssros, 'entganglement events'
        print 'of which',self.noof_valid_ssros ,\
                'are double click events according to the HH'
        if self.noof_valid_ssros >  self.noof_hh_ssros:
            print 'w',self.w
            print '-------------------------------'
            if (len(self.hhp)>0):
                print 'hhp', self.hhp[np.logical_and(self.hhp[:,3]==1, 
                            self.hhp[:,2]==2)]

            raise Exception('more valid than total plu events!')

    def save(self,savefile):
        np.savez(savefile, corr=self.uncond_corr, 
               corr_psi1=self.corr_psi1,
               corr_psi2=self.corr_psi2,
               windows=self.w)
      

class LDEAnalysis:
    """This class defines an analysis of an LDE experiment consisting 
    of various runs, subruns etc. The main methoud lde_analysis() can be called 
    with a number of runs, subruns, and the resulting analysis these will be
    added to previous runs/subruns. The constructor can be given an existing
    LDEAnalysis datafile to add more runs to"""
    
    
    def __init__(self):
        
        self.all_events=np.zeros((0,4),dtype=np.uint32)
        self.all_double_clicks=np.zeros((0,4),dtype=np.uint32)
        self.all_photon_hist_x=None
        self.all_photon_hist_w1_h0y=None
        self.all_photon_hist_w1_h1y=None
        self.all_photon_hist_w2_h0y=None
        self.all_photon_hist_w2_h1y=None
               
        self.all_window_data=[]
        self.all_gate_phases=[]
        self.all_ssro1=[]
        self.all_ssro2=[]
        self.all_CR1=[]
        self.all_CR2=[]
        
        self.all_tails=[]
        self.all_laser=[]
        self.all_statistics_lt1=[]
        self.all_statistics_lt2=[]
        self.all_statistics_plu=[]
        self.all_double_click_statistics=[]
        
        self.save_statistics_lt1=['get_noof_CR_checks','get_probe_counts']
        self.save_statistics_lt2=['get_noof_seq_starts','get_noof_PLU_markers',
                                        'get_noof_CR_checks','get_probe_counts']
    
        self.total_corr=np.zeros(4,dtype=np.uint32)
        self.total_corr_00=np.zeros(4,dtype=np.uint32)
        self.total_corr_01=np.zeros(4,dtype=np.uint32)
        self.total_corr_10=np.zeros(4,dtype=np.uint32)
        self.total_corr_11=np.zeros(4,dtype=np.uint32)
        self.savename=''
        self.savedir=''
        self.all_runs=[]
        self.all_subruns=[]

     
    def analyse_lde_from_dir(self, startpath, **kw): 
        """This convenience method calls most of the interesting analysis 
        for an LDE experiment defined 'hhp-data-' files in a directory tree"""

        identifier = 'rawdata'
        runs=[]
        subruns=[]
        savedir=kw.pop('savedir',os.path.join(startpath,'analysis'))
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        for (path, dirs, files) in os.walk(startpath):
            for drn in dirs:
                if drn.find(identifier)==0: 
                    print drn
                    print path
                    idx=int(drn[-sscorr.DATIDXDIGITS:])
                    try:
                        run_nr = runs.index(path)
                        subruns[run_nr].append(idx)
                    except ValueError:
                        runs.append(path)
                        #run_nr=len(runs)-1
                        subruns.append([idx])
        print runs
        print subruns
        return self.analyse_lde(runs,subruns,savedir=savedir, **kw)
            
    def analyse_lde_from_runfile(self, runfile, **kw):
        """This convenience method calls most of the interesting analysis 
        for an LDE experiment defined by a runfile generated by the analysis deamon"""
        
        d=np.load(runfile)
        runs=d['runs']
        subruns=d['subruns']
        d.close()
        savedir=kw.pop('savedir',
                os.path.join(os.path.dirname(runfile),'analysis'))
        if not os.path.isdir(savedir):
            os.mkdir(savedir)
        return self.analyse_lde(runs,subruns,savedir=savedir, **kw)
    
    def reanalyse_lde(self, w_start = (637,666), 
            w_length=150, w_dt=-1, ch1_offset=-29, w_dt_min=-1, **kw):
        
        anal = SingleLDEAnalysis()

        anal.w1_start = w_start
        anal.w2_start = kw.get('w2_start', w_start)
        w1_length=w_length
        w2_length= kw.get('w2_length', w_length)
        anal.w1_stop = (anal.w1_start[0]+w1_length,anal.w1_start[1]+w1_length)
        anal.w2_stop = (anal.w2_start[0]+w2_length,anal.w2_start[1]+w2_length)
        anal.dt_max=w_dt #NOTE with this number you can adjust the time difference between the photons
        anal.dt_min=w_dt_min
        anal.ch1_offset=ch1_offset#-4
        
        filter_ap_length=kw.get('filter_ap', 0)
        anal.single_photon_range_start=(anal.w1_start[0]-filter_ap_length,anal.w1_start[1]-filter_ap_length)
        anal.single_photon_range_stop=anal.w1_stop 
        anal.hhp=self.all_events
        anal.ssro1=self.all_ssro1
        anal.ssro2=self.all_ssro2
        anal.correlations()
        apply_to_self=kw.pop('apply_to_self',True)
        if apply_to_self:
            self.total_corr,self.total_corr_00,self.total_corr_01,self.total_corr_10, self.total_corr_11 = \
                    anal.uncond_corr, anal.corr_00, anal.corr_01, anal.corr_10, anal.corr_11
        return anal.uncond_corr, anal.corr_00, anal.corr_01, anal.corr_10, anal.corr_11
        

    def analyse_lde(self, runs, subruns, savedir='', w_start = (637,666), 
            w_length=150, w_dt=-1, ch1_offset=-29, w_dt_min=-1, **kw):
        """
        This convenience method calls most of the interesting analysis 
        for an LDE experiment defined by an array of directories of LDE 'runs', 
        and a 2D-array with for each run the integers defining which subruns 
        to use. 
        The results are saved in the savedir, which defaults to the current 
        working directory.
        
        params: 
        - w_start(ch0,ch1) in bins
        - w_length in bins
        - w_dt max distance between clicks in bins
        - kw: * w2_start, w2_length = w_start, w_length
              * filter on afterpulsing, given as the length in bins before w_start
              to look for laser photons: filter_ap=0
              * 

        returns:
        SingleLDEAnalysis object of last subrun
        """
        do_analyse_g2=kw.get('analyse_g2',True)
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        self.savedir=savedir
        anal = SingleLDEAnalysis()

        anal.w1_start = w_start
        anal.w2_start = kw.get('w2_start', w_start)
        w1_length=w_length
        w2_length= kw.get('w2_length', w_length)
        anal.w1_stop = (anal.w1_start[0]+w1_length,anal.w1_start[1]+w1_length)
        anal.w2_stop = (anal.w2_start[0]+w2_length,anal.w2_start[1]+w2_length)
        anal.dt_max=w_dt #NOTE with this number you can adjust the time difference between the photons
        anal.dt_min=w_dt_min
        anal.ch1_offset=ch1_offset#-4
        
        filter_ap_length=kw.get('filter_ap', 0)
        anal.single_photon_range_start=(anal.w1_start[0]-filter_ap_length,anal.w1_start[1]-filter_ap_length)
        anal.single_photon_range_stop=anal.w1_stop
        

        if self.all_photon_hist_x == None:
            self.all_photon_hist_x=np.zeros(CHMAXTIME,dtype=np.uint32)
            self.all_photon_hist_w1_h0y=np.zeros(CHMAXTIME,dtype=np.uint32)
            self.all_photon_hist_w1_h1y=np.zeros(CHMAXTIME,dtype=np.uint32)
            self.all_photon_hist_w2_h0y=np.zeros(CHMAXTIME,dtype=np.uint32)
            self.all_photon_hist_w2_h1y=np.zeros(CHMAXTIME,dtype=np.uint32)


        self.savename = 'window1_ch0_%d-%d_ch1_%d-%d__window2_ch0_%d-%d_ch1_%d-%d-%d-dt_max-%d-dt_min-run1-4' % \
                    (anal.w1_start[0], anal.w1_stop[0], anal.w1_start[1], 
                            anal.w1_stop[1],anal.w2_start[0], anal.w2_stop[0], 
                            anal.w2_start[1], anal.w2_stop[1], anal.dt_max, 
                            anal.dt_min)
        
        
        
        self.g2_deltas = range(-2,3)
        self.g2_coincidences_all = [np.array([],dtype=int) \
                for d in self.g2_deltas]
        self.g2_coincidences_corr = [np.array([],dtype=int) \
                for d in self.g2_deltas]
        self.g2_coincidences_fullwindow = [np.array([],dtype=int) \
                for d in self.g2_deltas]
        
        for _i, datfolder in np.ndenumerate(runs):
            self.all_runs.append(datfolder)
            i=_i[0]
            self.all_subruns.append(subruns[i])
            for DATIDX in subruns[i]:
                DATIDXSTR=('-%.'+str(sscorr.DATIDXDIGITS)+'d') % DATIDX
                anal.datfolder = datfolder
                anal.rawfolder = os.path.join(datfolder, 'rawdata'+DATIDXSTR)
                hhpdata=os.path.join(datfolder, HHPFILEBASE+DATIDXSTR+'.npz')
                corrsavefile = os.path.join(datfolder, 
                        'corr_data_'+self.savename+DATIDXSTR)
                            
                ### data preparation or loading
                if not os.path.exists(hhpdata):
                    #print hhpdata
                    print 'run', datfolder, ', subrun', DATIDXSTR,\
                            'not yet prepared, starting preparation'
                    anal.prepare_data(DATIDX=DATIDX)  
                else:
                    print 'loading from hhp-data'
                    d=np.load(hhpdata)
                    anal.hhp=d['hhp']
                    anal.hhdata=d['hhdata']
                    anal.w1=d['w1']
                    anal.w2=d['w2']
                    d.close()
                if len(anal.hhdata)==0:
                    print 'run', datfolder, ', subrun', DATIDXSTR, 'EMPTY'
                    anal.hhdata=np.array([np.zeros(4,dtype=np.uint32)])
                    anal.hhp=np.array([np.zeros(4,dtype=np.uint32)])
                    anal.w1=np.array([np.zeros(4,dtype=np.uint32)])
                    anal.w2=np.array([np.zeros(4,dtype=np.uint32)])
                 
                ### adwin data 
                anal.get_adwin_data(DATIDX=DATIDX)
                
                ###statistics
              
                anal.get_adwin_statistics(DATIDX)            
                statistics_lt1 = [anal.adwin_lt1_stats[key] for key in self.save_statistics_lt1]
                statistics_lt2 = [anal.adwin_lt2_stats[key] for key in self.save_statistics_lt2]
                self.all_statistics_lt1.append(statistics_lt1)
                self.all_statistics_lt2.append(statistics_lt2)
                #print statistics_lt1
                
                ### g2
                if do_analyse_g2:
                    c_all = anal.get_g2_coincidences(self.g2_deltas, 
                            ch0window=(0,500), ch1window=(0,500))
                    c_corr = anal.get_g2_coincidences(self.g2_deltas)
                    c_fullwindow = anal.get_g2_coincidences(self.g2_deltas,
                            ch0window=(anal.w1_start[0],500), ch1window=(anal.w1_start[1],500))
                
                    self.g2_coincidences_all = tpqi.add_coincidences(
                            self.g2_coincidences_all, c_all)
                    self.g2_coincidences_corr = tpqi.add_coincidences(
                            self.g2_coincidences_corr, c_corr)
                    self.g2_coincidences_fullwindow = tpqi.add_coincidences(
                            self.g2_coincidences_fullwindow, c_fullwindow)
                
                ### the correlations            
                anal.correlations()
                self.total_corr+=anal.uncond_corr ##NOTE has to be brought into right order  
#                self.total_corr_psi1+=anal.corr_psi1
#                self.total_corr_psi2+=anal.corr_psi2
                self.total_corr_00+=anal.corr_00
                self.total_corr_01+=anal.corr_01
                self.total_corr_10+=anal.corr_10
                self.total_corr_11+=anal.corr_11

                ###save single analysis
                #anal.save(corrsavefile)

                if len(anal.w)>0:
                    if (len(anal.w[np.where(anal.w>0)])>0):
                        event_photons=anal.get_event_photons(anal.w,anal.hhp)
                       # print 'event photons', event_photons
                        self.all_events=np.vstack((self.all_events,event_photons))
                    x, w1_h0y, w1_h1y, w2_h0y, w2_h1y =\
                            self._get_time_histogram(anal.w1,anal.w2, CHMAXTIME,
                                    (0,CHMAXTIME-1))
                    self.all_photon_hist_x=x
                    self.all_photon_hist_w1_h0y+=w1_h0y
                    self.all_photon_hist_w1_h1y+=w1_h1y
                    self.all_photon_hist_w2_h0y+=w2_h0y
                    self.all_photon_hist_w2_h1y+=w2_h1y

                ###double_clicks check
                raw_double_clicks=sscorr.get_double_clicks(anal.hhdata,ch0_start=anal.w1_start[0], 
                        ch0_stop=anal.w1_start[0]+w1_length, \
                        ch1_start=anal.w1_start[1], 
                        ch1_stop=anal.w1_start[1]+w1_length)
                self.all_double_clicks=np.vstack((self.all_double_clicks,raw_double_clicks))
                double_click_statistics=sscorr.get_channel_statistics(raw_double_clicks)
                self.all_double_click_statistics.append(double_click_statistics)
                anal.print_plu_vs_hh()
                print 'The HH detected in total', len(raw_double_clicks)/2,\
                        'valid double click events'
                self.all_statistics_plu.append([anal.noof_hh_ssros, 
                    anal.noof_valid_ssros,len(raw_double_clicks)/2]) 
                ###tail_cts per run
                anal.tail_counts()
                self.all_tails = np.append(self.all_tails, anal.total_tail)
                anal.laser_counts()
                self.all_laser = np.append(self.all_laser, anal.total_laser)

                ### adding to all_arrays
                #NOTE when using indexing arrays, make sure that
                # a) the indexing array is not empty (or only consists only of False elements)
                # b) the array to be indexed is not empty either
                self.all_window_data =  np.append(self.all_window_data,anal.w[np.where(anal.w>0)] 
                      if len(anal.w) > 0 else [])
                self.all_gate_phases = np.append(self.all_gate_phases, 
                        anal.gate_phase[np.where(anal.w>0)] if len(anal.gate_phase) > 0 else [])
                self.all_ssro1 = np.append(self.all_ssro1,
                        anal.ssro1[np.where(anal.w>0)] if len(anal.ssro1) > 0 else [] )
                self.all_ssro2 = np.append(self.all_ssro2,
                        anal.ssro2[np.where(anal.w>0)] if len(anal.ssro2) > 0 else [] )
                self.all_CR1 = np.append(self.all_CR1,
                        anal.cr1[np.where(anal.w>0)] if len(anal.cr1) > 0 else [] )
                self.all_CR2 = np.append(self.all_CR2,
                        anal.cr2[np.where(anal.w>0)] if len(anal.cr2) > 0 else [] )
                
                
        print self.total_corr
        self.last_anal=anal
        return anal
        
    def save(self, filename=''):        
        if filename=='': 
            filename=os.path.join(self.savedir,'analysis_'+self.savename+'.pkl')
        #print filename
        f=open(filename,'wb')
        pickle.dump(self,f)
        f.close()


    def filter_on_gatephase(self, **kw):
        """returns the total correlations filtered on good gatephase. 
        kw: apply_to_self : boolean apply the filter to the total 
        correlations total_corr_ij"""
        if len(np.where(self.all_gate_phases>0)[0])==0:
            return [],[],[]
        gate_phases=self.all_gate_phases[self.all_gate_phases!=0]
        fltr = gate_phases>0
        
        return self.filter_correlations(fltr, **kw)

    def filter_on_CR(self, lt1_min, lt2_min, **kw):
        """returns the total correlations filtered on CR checks after RO. 
        Minimum counts during CR check for lt1 and lt2 should be given as arguments.
        kw: apply_to_self : boolean apply the filter to the total 
        correlations total_corr_ij"""

        if len(np.where(np.logical_and(self.all_CR1>lt1_min, self.all_CR2>lt2_min))[0])==0:
            return [],[],[]
        fltr = np.logical_and(self.all_CR1>lt1_min, self.all_CR2>lt2_min)
        return self.filter_correlations(fltr, **kw)
        
    def filter_correlations(self, filter, **kw):
        """returns the total correlations filtered on any Boolean array of length len(self.ssro1)
        kw: apply_to_self : boolean apply the filter to the total 
        correlations total_corr_ij"""
        apply_to_self=kw.pop('apply_to_self',True)
        if apply_to_self:
            self.total_corr,self.total_corr_00,self.total_corr_01,self.total_corr_10, self.total_corr_11 = \
            sscorr.correlations(self.all_ssro1[filter],
                self.all_ssro2[filter], 
                self.all_window_data[filter])
        return sscorr.correlations(self.all_ssro1[filter],
                self.all_ssro2[filter], 
                self.all_window_data[filter])

                
    def plot_correlations(self, F0LT2 = 0.805, F0LT1 = 0.905, F1LT2 = 0.998, F1LT1 = 0.9937,
            save_plots=True, **kw):

        """ Plots the correlations"""
        total_corr=kw.get('total_corr', self.total_corr) 
        total_corr_psi1=kw.get('total_corr_psi1', self.total_corr_00 + self.total_corr_11) 
        total_corr_psi2=kw.get('total_corr_psi2', self.total_corr_01 + self.total_corr_10)
        
        total_corr_err=sscorr.get_correlation_errors(total_corr)
        total_corr_psi1_err=sscorr.get_correlation_errors(total_corr_psi1)
        total_corr_psi2_err=sscorr.get_correlation_errors(total_corr_psi2)
      
        corrected_corr, corr_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(total_corr, 
                F0LT1, F0LT2, F1LT1, F1LT2, verbose = True, return_error_bars=True)
        corrected_corr_psi1, corr_psi1_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(total_corr_psi1, 
                F0LT1, F0LT2, F1LT1, F1LT2, verbose = True, return_error_bars=True)
        corrected_corr_psi2, corr_psi2_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(total_corr_psi2, 
                F0LT1, F0LT2, F1LT1, F1LT2, verbose = True, return_error_bars=True)

        sscorr.plot_uncorrected(total_corr, total_corr_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr'))
        sscorr.plot_uncorrected(total_corr_psi1, total_corr_psi1_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr_psi1'))
        sscorr.plot_uncorrected(total_corr_psi2, total_corr_psi2_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr_psi2'))
                 
        sscorr.plot_corrected(corrected_corr, corr_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr'))
        sscorr.plot_corrected(corrected_corr_psi1, corr_psi1_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr_psi1'))
        sscorr.plot_corrected(corrected_corr_psi2, corr_psi2_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr_psi2'))

        ### g2 plots
        if save_plots:
            tpqi.plot_g2(self.g2_deltas, self.g2_coincidences_all,
                delta_separation_time=200, 
                savepath=os.path.join(self.savedir, 'g2_all.png'))
            tpqi.plot_g2(self.g2_deltas, self.g2_coincidences_corr,
                delta_separation_time=200, 
                savepath=os.path.join(self.savedir, 'g2_correlationwindow.png'))
            tpqi.plot_g2(self.g2_deltas, self.g2_coincidences_fullwindow,
                delta_separation_time=200, 
                savepath=os.path.join(self.savedir, 'g2_fulltail.png'))
                
    def get_run_subrun_from_datapoint_idx(self,idx):
        """returns the measurement run (directory string) and subrun index 
        corresponding to a total subrun number in the analysis"""
        number=int(idx)
        totsubs=0
        run=0
        for subruns in self.all_subruns:
            if totsubs+len(subruns)>=number:
                break
            run+=1
            totsubs+=len(subruns)
        #print 'run', run, 'number', number, 'totsubs', totsubs
        return list(cbook.flatten(self.all_runs))[run], list(cbook.flatten(self.all_subruns))[idx]
    
    def plot_tailcts(self,save_plots=True):
        """plots total subrun number versus tailcounts and laser counts"""
        fig=plt.figure()
        plt.title('Combined counts per seq-start*300 from both \
                channels+windows')
        ax=plt.subplot(211)
        ax.plot(np.arange(len(list(cbook.flatten(self.all_subruns)))),self.all_tails)
        ax.set_xlabel('Subrun total #')
        ax.set_ylabel('Tailcounts')
        ax=plt.subplot(212)
        ax.plot(np.arange(len(list(cbook.flatten(self.all_subruns)))),self.all_laser)
        ax.set_ylabel('Lasercounts')
        
        if save_plots:
            fig.savefig(os.path.join(self.savedir,self.savename)+'tailcounts.png')
        
    def plot_statistics(self,setup='lt2',save_plots=True,**kw):
        """plots total subrun number versus adwin statistics saved 
        from in the lde_analysis
        pars:
        - setup: h=give setup to plot statistics for
        - known kw: plot_statistics: list of statistics to plot, 
                    defaults to all saved satistics for given setup.
                    note: statistics in plot_statistics must be 
                    contained in self.save_statistics for the given setup.
        """
        if setup=='lt1':
            stats=self.all_statistics_lt1
            saved_stats=self.save_statistics_lt1
            plot_statistics=kw.get('plot_statistics',self.save_statistics_lt1)
        elif setup=='lt2':
            stats=self.all_statistics_lt2
            saved_stats=self.save_statistics_lt2
            plot_statistics=kw.get('plot_statistics',self.save_statistics_lt2)
        else: 
            print 'unknown setup'
            return
        num_stat=len(plot_statistics)
        i=0
        fig=plt.figure()
        for stat in plot_statistics:
            try:
                stat_column=saved_stats.index(stat)
            except ValueError:
                print 'statistic', stat, 'not found in analysed statistcs'
                continue
            i+=1
            ax=plt.subplot(num_stat,1,i)
            ax.plot(np.arange(len(list(cbook.flatten(self.all_subruns)))),
                    [subrun[stat_column] for subrun in stats])
            ax.set_xlabel('Subrun total #')
            ax.set_ylabel(stat)
        if save_plots:
            fig.savefig(os.path.join(self.savedir,self.savename)+'statistics_'+setup+'.png')
    
    def plot_plu_vs_tail_all(self,range=None,bins=None,bins_plu=30,**kw):
        """Plots a histogram of the arrival times of the gated plu events vs the
        raw tail data of the last subrun
        keyword args:   range=last_anal.w1_(start,stop), bins=full_range,
                        bins_plu=30, log_plots=False, save_plots=False, save_path=''   
        """
        w1_photons_plu, w2_photons_plu =\
                        self._get_double_click_windows(self.all_events,2)
        last_anal=self.last_anal
        if range == None:
            range=(last_anal.w1_start[0], last_anal.w1_stop[0])
        if bins == None:
            bins = range[1]-range[0]
        #print 'length x before:', len(self.all_photon_hist_x)
        #print 'length x:', len(self.all_photon_hist_x[range[0]:range[1]])
        #print 'range:', range
        #print 'length y:', len(self.all_photon_hist_w1_h1y[range[0]:range[1]])
        self._plot_windowed_photons(w1_photons_plu, w2_photons_plu,
                    bins=bins_plu,
                    range=range,
                    x2=self.all_photon_hist_x[range[0]:range[1]+1],
                    w1_h0y_2=self.all_photon_hist_w1_h0y[range[0]:range[1]], 
                    w1_h1y_2=self.all_photon_hist_w1_h1y[range[0]:range[1]], 
                    w2_h0y_2=self.all_photon_hist_w2_h0y[range[0]:range[1]], 
                    w2_h1y_2=self.all_photon_hist_w2_h1y[range[0]:range[1]],
                    c1='r', c2='g',
                     **kw)

    def plot_plu_vs_tail_last_subrun(self,**kw):
        """Plots a histogram of the arrival times of the gated plu events vs the
        raw tail data of the last subrun
        keyword args:   range=last_anal.w1_(start,stop),bins=full_range,
                        bins_plu=30, log_plots=False, save_plots=False, save_path=''   
        """
        w1_photons_plu, w2_photons_plu =\
                        self._get_double_click_windows(self.all_events,2)
        last_anal=self.last_anal
        range=kw.pop('range',(last_anal.w1_start[0], last_anal.w1_stop[0]))
        bins = kw.pop('bins',range[1]-range[0])
        bins_plu= kw.pop('bins_plu',30)

        self._plot_windowed_photons(last_anal.w1, last_anal.w2, 
                bins=bins,
                range=range,
                w1_photons_2=w1_photons_plu, w2_photons_2=w2_photons_plu, 
                bins_2=bins_plu, **kw)

    def plot_plu_vs_hh_events(self,**kw):
        """Plots a histogram of the arrival times of the gated plu events vs the
        double click events as seen by the hydraharp
        keyword args: - all_plu_events, all_hh_clicks 
                        (default to self.all_events, self.all_double_clicks)
                      - bins=700, range=(0,699), bins_2=700, 
                        log_plots=False, save_plots=False, save_path=''   
        """
        all_plu_events=kw.pop('all_plu_events',self.all_events)
        all_hh_clicks=kw.pop('all_hh_clicks',self.all_double_clicks)
        
        w1_photons_plu, w2_photons_plu =\
                self._get_double_click_windows(all_plu_events,2)
        w1_photons_hh, w2_photons_hh =\
                self._get_double_click_windows(all_hh_clicks)

        self._plot_windowed_photons(w1_photons_hh,w2_photons_hh,
                                   w1_photons_2=w1_photons_plu,
                                   w2_photons_2=w2_photons_plu,
                                   **kw)

    def plot_plu_good_vs_bad_events(self,good_ssro=[1,2],**kw):
        """Plots a histogram of the arrival times of the gated plu events,
        with one axis good events, and on the other the bad events
        as defined by the good_ssro as a list of integers 0-3 [00,01,10,11]
        - **kw:  bins=700, range=(0,699), bins_2=700, 
                        log_plots=False, save_plots=False, save_path=''
        """
        
        all_plu_events=kw.pop('all_plu_events',self.all_events)
        all_ssro1=kw.pop('all_ssro1', self.all_ssro1)
        all_ssro2=kw.pop('all_ssro2', self.all_ssro2)


        good_events, bad_events=\
                self.get_good_bad_events(all_plu_events,all_ssro1,all_ssro2, good_ssro)
        print 'good:', len(good_events)/3, 'bad:', len(bad_events)/3
        w1_photons_good, w2_photons_good = \
                self._get_double_click_windows(good_events,2)
        w1_photons_bad, w2_photons_bad = \
                self._get_double_click_windows(bad_events, 2)

        self._plot_windowed_photons(w1_photons_good,w2_photons_good,
                                   w1_photons_2=w1_photons_bad,
                                   w2_photons_2=w2_photons_bad,
                                   **kw)

    def get_good_bad_events(self,all_plu_events,all_ssro1,all_ssro2,good_ssro):
        """Divides hhplu data into good and bad by the ssro correlation they have
        as defined by the good_ssro as a list of integers 0-3 [00,01,10,11]"""
        
        bad_events= np.zeros((0,4),dtype=np.uint32)
        good_events= np.zeros((0,4),dtype=np.uint32)
        for i in range(len(all_ssro1)):
            corr=sscorr.correlations(np.array([all_ssro1[i]]),
                    np.array([all_ssro2[i]]),np.array([1]))[0]
            if (1 in [corr[good] for good in good_ssro]):
                good_events=np.vstack((good_events,all_plu_events[3*i]))
                good_events=np.vstack((good_events,all_plu_events[3*i+1]))
                good_events=np.vstack((good_events,all_plu_events[3*i+2]))
            else:
                bad_events=np.vstack((bad_events,all_plu_events[3*i]))
                bad_events=np.vstack((bad_events,all_plu_events[3*i+1]))
                bad_events=np.vstack((bad_events,all_plu_events[3*i+2]))

        return good_events, bad_events


    def plot_dts(self,hrange=(-200,200), bins=700, good_ssro=[1,2], save_fig=False, save_path='', **kw):
        """Plots a histogram of the dt's (w2 arrival time - w1 arrival time) 
        of the gated plu events, with one axis good events, and on the other 
        the bad events as defined by the good_ssro as a list of 
        integers 0-3 [00,01,10,11]"""
        
        all_plu_events=kw.pop('all_plu_events',self.all_events)
        all_ssro1=kw.pop('all_ssro1', self.all_ssro1)
        all_ssro2=kw.pop('all_ssro2', self.all_ssro2)

        good_events, bad_events=self.get_good_bad_events(all_plu_events,all_ssro1,all_ssro2, good_ssro)
        dts_good=np.array([])
        for i in range(int(len(good_events)/3)):
            dt=int(good_events[3*i+1][1])-int(good_events[3*i][1])
            dts_good=np.append(dts_good,dt)
        
        dts_bad=np.array([])
        for i in range(int(len(bad_events)/3)):
            dt=int(bad_events[3*i+1][1])-int(bad_events[3*i][1])
            dts_bad=np.append(dts_bad,dt)


        
        plt.figure()
        ax = plt.subplot(111)
        plt.title('dts for all plu events')
        plt.hist(dts_good,range=hrange, bins=bins, facecolor='none', edgecolor= 'g' ,
                    hatch = '\\')
        ax.set_xlabel('Bins')
        ax2=ax.twinx()
        plt.hist(dts_bad,range=hrange, bins=bins, facecolor='none', edgecolor= 'r' ,
                    hatch = '//')

        if save_fig:
            fig1.savefig(save_path+'photons_window1.png')

        

    def _get_double_click_windows(self, double_clicks, delete_marker=None):
        if delete_marker!=None:
            double_clicks=hht3.delete_markers(double_clicks,delete_marker)
        w1_photons=np.zeros((0,4),dtype=np.uint32)
        w2_photons=np.zeros((0,4),dtype=np.uint32)

        for i in np.arange(len(double_clicks)):
            if (i%2)==0:
                w1_photons=np.vstack((w1_photons,double_clicks[i]))
            if (i%2)==1:
                w2_photons=np.vstack((w2_photons,double_clicks[i]))
        
        return w1_photons, w2_photons 

    def _get_time_histogram(self,w1_photons,w2_photons, bins, range):
        
        w1_ch0,w1_ch1=hht3.get_click_times(w1_photons)
        w2_ch0,w2_ch1=hht3.get_click_times(w2_photons)

        #print 'w1_ch0',len(w1_ch0),'w1_ch1',len(w1_ch1)
        #print  'w2_ch0',len(w2_ch0), 'w2_ch0',len(w2_ch1)

        w1_h0y, _tmp = np.histogram(w1_ch0, bins=bins, range=range)
        w1_h1y, _tmp = np.histogram(w1_ch1, bins=bins, range=range)

        w2_h0y, _tmp = np.histogram(w2_ch0, bins=bins, range=range)
        w2_h1y, x = np.histogram(w2_ch1, bins=bins, range=range)

        return  x, w1_h0y, w1_h1y, w2_h0y, w2_h1y

    def _plot_windowed_photons(self,w1_photons,w2_photons,bins=700,range=(0,699),
            w1_photons_2=None, w2_photons_2=None,bins_2=700, **kw):
        
        two_axis_plot= w1_photons_2!=None
       
        print 'number of photons in Window 1', len(w1_photons)
        print 'number of photons in Window 2', len(w2_photons)
        #print 'number of w1 ch0' , len(w1_photons[np.where(w1_photons[:,2]==0)])
        #print 'number of w1 ch1' , len(w1_photons[np.where(w1_photons[:,2]==1)])
        #print 'number of w2 ch0' , len(w2_photons[np.where(w2_photons[:,2]==0)])
        #print 'number of w2 ch1' , len(w2_photons[np.where(w2_photons[:,2]==1)]) 
        if(two_axis_plot):
            print '2 number of photons in Window 1', len(w1_photons_2)
            print '2 number of photons in Window 2', len(w2_photons_2)
            #print '2 number of w1 ch0' , len(w1_photons_2[np.where(w1_photons_2[:,2]==0)])
            #print '2 number of w1 ch1' , len(w1_photons_2[np.where(w1_photons_2[:,2]==1)])
            #print '2 number of w2 ch0' , len(w2_photons_2[np.where(w2_photons_2[:,2]==0)])
            #print '2 number of w2 ch1' , len(w2_photons_2[np.where(w2_photons_2[:,2]==1)]) 

        x, w1_h0y, w1_h1y, w2_h0y, w2_h1y = \
                self._get_time_histogram(w1_photons,w2_photons, bins, range)
        
        if(two_axis_plot):
            x2, w1_h0y_2, w1_h1y_2, w2_h0y_2, w2_h1y_2 =\
                    self._get_time_histogram(w1_photons_2,w2_photons_2, bins_2, range)
            self._plot_histograms(x, w1_h0y, w1_h1y, w2_h0y, w2_h1y, 
                       x2, w1_h0y_2, w1_h1y_2, w2_h0y_2, w2_h1y_2, **kw)
        else:
            self._plot_histograms(x, w1_h0y, w1_h1y, w2_h0y, w2_h1y, **kw)

    def _plot_histograms(self, x, w1_h0y, w1_h1y, w2_h0y, w2_h1y, 
                       x2=None, w1_h0y_2=None, w1_h1y_2=None, w2_h0y_2=None, w2_h1y_2=None,
                       c1='g', c2= 'r',log_plots=False,
                       save_fig=False,save_path=''):
        """This plots 2 figures x 2 subplots x 2 y-axes in barplots, 
        as given by w1,w2_h0,h1_(x,y)_1,2 variables."""
        two_axis_plot= x2!=None
        if save_path=='':
            save_path=self.savedir
    
        fig1 = plt.figure()  

        ax = plt.subplot(121)
        plt.title('window 1, channel 0')
        plot_bar(x[:-1], w1_h0y, width=(max(x)-min(x))/len(x),  facecolor='none', 
                edgecolor= c1,label='channel 0', hatch = '/', log=log_plots)
        ax.set_xlabel('Bins')
        if two_axis_plot:
            ax2=ax.twinx()
            plot_bar(x2[:-1], w1_h0y_2, width=(max(x2)-min(x2))/len(x2),  
                    facecolor='none', edgecolor= c2 ,label='channel 0_2',
                    hatch = '\\', log=log_plots)

        ax = plt.subplot(122)
        plt.title('window 1, channel 1')
        plot_bar(x[:-1], w1_h1y, width=(max(x)-min(x))/len(x),  facecolor='none', 
                edgecolor= c1 ,label='channel 1', hatch = '/', log=log_plots)
        ax.set_xlabel('Bins')
        if two_axis_plot:
            ax2=ax.twinx()
            plot_bar(x2[:-1], w1_h1y_2, width=(max(x2)-min(x2))/len(x2),  
                    facecolor='none', edgecolor= c2, label='channel 1_2', 
                    hatch = '\\', log=log_plots)

        if save_fig:
            fig1.savefig(os.path.join(save_path,'photons_window1.png'))


        fig2 = plt.figure()

        ax = plt.subplot(121)
        plt.title('window 2, channel 0')
        plot_bar(x[:-1], w2_h0y, width=(max(x)-min(x))/len(x),  facecolor='none', 
                edgecolor= c1,label='channel 0', hatch = '/', log=log_plots)
        ax.set_xlabel('Bins')
        if two_axis_plot:
            ax2=ax.twinx()
            plot_bar(x2[:-1], w2_h0y_2, width=(max(x2)-min(x2))/len(x2),  
                    facecolor='none', edgecolor= c2 ,label='channel 0_2', 
                    hatch = '\\', log=log_plots)

        ax = plt.subplot(122)
        plt.title('window 2, channel 1')
        plot_bar(x[:-1], w2_h1y, width=(max(x)-min(x))/len(x),  facecolor='none', 
                edgecolor= c1 ,label='channel 1', hatch = '/', log=log_plots)
        ax.set_xlabel('Bins')
        if two_axis_plot:
            ax2=ax.twinx()
            plot_bar(x2[:-1], w2_h1y_2, width=(max(x2)-min(x2))/len(x2), 
                    facecolor='none', edgecolor = c2,label='channel 1_2',
                    hatch = '\\', log=log_plots)

        if save_fig:
            fig2.savefig(os.path.join(save_path,'photons_window2.png'))
       
def plot_bar(x,y,**kw):
    if len(x)>0 and len(y)>0 and np.sum(y)>0:
        #print len(x), len(y), np.sum(y)
        plt.bar(x,y,**kw)
 
def load_previous_analysis(filename):
    f=open(filename,'rb')
    a = pickle.load(f)
    f.close()
    return a
       


""" Below are some typical uses of the above class."""

if __name__ == '__main__':
    
    ###option 1: The measurement is defined by a folder, 
                #we find all hhpdata-files in the folder an analyse them
    datadir=r'D:\Analysis\2012-09_ldetesting\autodata\20120905'
    a = LDEAnalysis()
    a.analyse_lde_from_dir(datadir, w_start = (234,229), w_length=150, w_dt=-1)

    ###option 2, measurement defined by runs, subruns:
    run1=r'D:\Analysis\2012-09_ldetesting\autodata\20120905\013113_LDE_Entanglement_ZZ'
    run2=r'D:\Analysis\2012-09_ldetesting\autodata\20120905\013113_LDE_Entanglement_ZZ'
    runs = [run1,run2]
    subruns = [np.arange(14),[1,3,4]]
    a = LDEAnalysis()
    a.analyse_lde(runs, subruns,savedir=r'D:\Analysis\2012-09_lde',
                    w_start = (234,229), w_length=150, w_dt=-1)

    #we can add some more data to the existing anaysis:
    another_run=r'D:\Analysis\2012-09_lde\2012-09-05\093113_LDE_Entanglement_ZZ'
    subruns=[np.arange(4)]
    a.analyse_lde([another_run], subruns,savedir=r'D:\Analysis\2012-09_lde',
                    w_start = (234,229), w_length=150, w_dt=-1)
    
    ### save the current analysis (saves all_photons, all_total_corrs, all gate_phases etc) 
    a.save(r'D:\Analysis\test_lde_analysis.npz')
    
    ### we can also work/add to an existing previously saved_analysis: 
    #a = load_previous_analysis('D:\Analysis\test_lde_analysis.npz')
    
    ###Plot the results and save the images
    a.plot(F0LT2 = 0.805, F0LT1 = 0.905, F1LT2 = 0.998, F1LT1 = 0.9937)
    
    ### do some more filtering
    gate_filtered_corr, gate_filtered_corr_psi1, gate_filtered_corr_psi2 = \
            a.filter_on_gatephase()
    cr_filtered_corr, cr_filtered_corr_psi1, cr_filtered_corr_psi2 = \
           a.filter_on_CR(40,10)
    
    ### plot the filtered total correlations, and leave out the psi-corrs.
    a.plot(total_corr=gate_filtered_corr,
    F0LT2 = 0.805, F0LT1 = 0.905, F1LT2 = 0.998, F1LT1 = 0.9937,
    xx_measurement=False)

    ###plot some more stuff:
    a.plot_tailcts()
    a.plot_statistics(setup='lt2')
    a.plot_statistics(setup='lt1')
    
    print'All:', a.total_corr, '\n 00:', a.total_corr_00,'\n 01:', a.total_corr_01, \
            '\n 10:', a.total_corr_10,'\n 11:', a.total_corr_11
    
    a.filter_correlations(a.all_laser<.015)
    a.total_corr,a.total_corr_00,a.total_corr_01,a.total_corr_10,a.total_corr_11=a.filter_on_gatephase()
    a.total_corr,a.total_corr_00,a.total_corr_01,a.total_corr_10,a.total_corr_11=a.filter_on_CR(10,1)
    a.plot_correlations()
    ###You can get your analysis object form a saved p[revious analysis
    s= load(r'D:\Analysis\test_lde_analysis.npz')






=======
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cbook
import pickle

from analysis.lib.fitting import fit,common
from analysis.lib.lde import spcorr, sscorr, tpqi
from analysis.lib.pq import hht3

HHPFILEBASE='hhp_data'
CHMAXTIME=2300


class SingleLDEAnalysis:

    datfolder = ''
    rawfolder = ''
 
    def __init__(self):
        pass        

    def prepare_data(self, ch0maxtime=CHMAXTIME, ch1maxtime=CHMAXTIME, DATIDX = 0):
        self.hhdata, self.w1, self.w2 = sscorr.load_hhdata(self.rawfolder, ch0maxtime=ch0maxtime,
                ch1maxtime=ch1maxtime, do_filter_crap=False)
        self.hhp, self.hhp1, self.hhp2 = sscorr.get_hhpludata(self.hhdata)

       
        self.ssro1, self.ssro2, self.cr1, self.cr2, self.gate_phase = sscorr.load_adwin_data(
                self.datfolder, DATIDX=DATIDX)
        if (len(self.hhp)>0):
            self.noof_hh_ssros = len(self.hhp[np.logical_and(self.hhp[:,3]==1, 
                self.hhp[:,2]==2)])
        else:
            self.noof_hh_ssros=0
        self.noof_adwin_ssros = len(self.ssro2)

        if self.noof_hh_ssros > self.noof_adwin_ssros:
            print 'HH SSROs do not match Adwin!'
            raise(Exception('hharp vs ssro mismatch'))

        np.savez(os.path.join(self.datfolder,HHPFILEBASE+('-%.'+str(sscorr.DATIDXDIGITS)+'d') % DATIDX),\
                hhdata=self.hhdata, hhp=self.hhp, w1=self.w1, w2=self.w2)
        
    def get_adwin_data(self, DATIDX = 0):
        self.ssro1, self.ssro2, self.cr1, self.cr2, self.gate_phase = sscorr.load_adwin_data(
                self.datfolder, DATIDX=DATIDX)  
        return True
    
    def get_dict_from_npz(self,d):
        e={}
        for key in d.keys():
            e[key]=d[key]
        return e
        
    def get_adwin_statistics(self,DATIDX):
        fn = os.path.join(self.datfolder, 
                sscorr.ADWINLT1FOLDER+('-%.'+str(sscorr.DATIDXDIGITS)+'d') % DATIDX, 
                sscorr.ADWINDATAFN)
        d = np.load(fn)
        self.adwin_lt1_stats=self.get_dict_from_npz(d)
        d.close()
        
        sn = os.path.join(self.datfolder, 
                sscorr.ADWINLT2FOLDER+('-%.'+str(sscorr.DATIDXDIGITS)+'d') % DATIDX, 
                sscorr.ADWINDATAFN)
        s = np.load(sn)
        self.adwin_lt2_stats=self.get_dict_from_npz(s)
        s.close()


    def correlations(self, *arg, **kw):
        if (len(self.hhp)>0):
            self.w = sscorr.windowdata(self.hhp, w1_start=self.w1_start, 
                    w1_stop=self.w1_stop, w2_start=self.w2_start, 
                    w2_stop=self.w2_stop, dt_max=self.dt_max,dt_min=self.dt_min,ch1_offset=self.ch1_offset,
                    single_photon_range_start=self.single_photon_range_start,
                    single_photon_range_stop=self.single_photon_range_stop,
                    *arg, **kw)
            self.uncond_corr, self.corr_00, self.corr_01, self.corr_10, self.corr_11 =\
                    sscorr.correlations(self.ssro1, self.ssro2, self.w)
        else:
            self.uncond_corr=np.zeros(4)

            self.corr_00=np.zeros(4)
            self.corr_01=np.zeros(4)
            self.corr_10=np.zeros(4)
            self.corr_11=np.zeros(4)
            self.w=[]

        return True
    
    def tail_counts(self):
        """calculates the tail counts - needs get_adwin_statistics to be run first"""
        w1tailcts_ch0=len(np.where(
            np.logical_and(self.w1[:,2]==0,np.logical_and(self.w1[:,1]>self.w1_start[0],self.w1[:,1]<self.w1_stop[0])))[0])
        w1tailcts_ch1=len(np.where(
            np.logical_and(self.w1[:,2]==1,np.logical_and(self.w1[:,1]>self.w1_start[1],self.w1[:,1]<self.w1_stop[1])))[0])
        w2tailcts_ch0=len(np.where(
            np.logical_and(self.w2[:,2]==0,np.logical_and(self.w2[:,1]>self.w2_start[0],self.w2[:,1]<self.w2_stop[0])))[0])
        w2tailcts_ch1=len(np.where(
            np.logical_and(self.w2[:,2]==1,np.logical_and(self.w2[:,1]>self.w2_start[1],self.w2[:,1]<self.w2_stop[1])))[0])
        seq_starts=self.adwin_lt2_stats['get_noof_seq_starts']
        if float(seq_starts)>0:
            self.total_tail=1.*(w1tailcts_ch0+w1tailcts_ch1+w2tailcts_ch0+w2tailcts_ch1)/(float(seq_starts)*300.)
        else:
            self.total_tail=0

    def laser_counts(self, laser_length=16):
        """calculates the laser counts - needs get_adwin_statistics to be run first"""
        w1tailcts_ch0=len(np.where(
            np.logical_and(self.w1[:,2]==0,np.logical_and(self.w1[:,1]>self.w1_start[0]-laser_length,self.w1[:,1]<self.w1_start[0])))[0])
        w1tailcts_ch1=len(np.where(
            np.logical_and(self.w1[:,2]==1,np.logical_and(self.w1[:,1]>self.w1_start[1]-laser_length,self.w1[:,1]<self.w1_start[1])))[0])
        w2tailcts_ch0=len(np.where(
            np.logical_and(self.w2[:,2]==0,np.logical_and(self.w2[:,1]>self.w2_start[0]-laser_length,self.w2[:,1]<self.w2_start[0])))[0])
        w2tailcts_ch1=len(np.where(
            np.logical_and(self.w2[:,2]==1,np.logical_and(self.w2[:,1]>self.w2_start[1]-laser_length,self.w2[:,1]<self.w2_start[1])))[0])
        seq_starts=self.adwin_lt2_stats['get_noof_seq_starts']
        if float(seq_starts)>0:
            self.total_laser=1.*(w1tailcts_ch0+w1tailcts_ch1+w2tailcts_ch0+w2tailcts_ch1)/(float(seq_starts)*300.)
        else:
            self.total_laser=0
        

###NOTE does not consider gate phase yet
    def get_event_photons(self, windowdata, hhp):
        plu_idx=np.where(hhp[:,3]==1)[0][np.where(windowdata>0)]
        event_photons= np.zeros((0,4),dtype=np.uint32)
        
        for i in plu_idx:
            event_photons=np.vstack((event_photons,hhp[i-2]))
            event_photons=np.vstack((event_photons,hhp[i-1]))
            event_photons=np.vstack((event_photons,hhp[i]))
        return event_photons

         
    def get_g2_coincidences(self, deltas, **kw):
        # plotsavename = 'g2_ch0-%d-%d_ch1-%d-%d.png'
        data = kw.get('hhdata',self.hhdata)
        ch0window=kw.get('ch0window',(self.w1_start[0], self.w1_stop[0]))
        ch1window=kw.get('ch1window',(self.w1_start[1], self.w1_stop[1]))

        data = hht3.filter_timewindow(data, 0, 
                mintime=ch0window[0], maxtime=ch0window[1])
        data = hht3.filter_timewindow(data, 1, 
                mintime=ch1window[0], maxtime=ch1window[1])

        d,c = tpqi.coincidences(data, deltas)
        # tpqi.plot_g2(d,c,savepath=os.path.join(self.savedir,plotsavename))
        return c


    def print_SSRO(self, DATIDX=0):
        self.ssro1, self.ssro2, self.cr1, self.cr2 = \
                sscorr.load_adwin_data(self.datfolder, DATIDX=DATIDX)
        print self.ssro2

    def print_plu_vs_hh(self):
        if (len(self.hhp)>0):
            self.noof_hh_ssros= len(self.hhp[np.logical_and(self.hhp[:,3]==1, 
                    self.hhp[:,2]==2)])
        else:
            self.noof_hh_ssros=0
        if len(self.w)>0:    
            self.noof_valid_ssros=len(np.where(self.w>0)[0])
        else:
            self.noof_valid_ssros=0
        print 'Plu gave', self.noof_hh_ssros, 'entganglement events'
        print 'of which',self.noof_valid_ssros ,\
                'are double click events according to the HH'
        if self.noof_valid_ssros >  self.noof_hh_ssros:
            print 'w',self.w
            print '-------------------------------'
            if (len(self.hhp)>0):
                print 'hhp', self.hhp[np.logical_and(self.hhp[:,3]==1, 
                            self.hhp[:,2]==2)]

            raise Exception('more valid than total plu events!')

    def save(self,savefile):
        np.savez(savefile, corr=self.uncond_corr, 
               corr_psi1=self.corr_psi1,
               corr_psi2=self.corr_psi2,
               windows=self.w)
      

class LDEAnalysis:
    """This class defines an analysis of an LDE experiment consisting 
    of various runs, subruns etc. The main methoud lde_analysis() can be called 
    with a number of runs, subruns, and the resulting analysis these will be
    added to previous runs/subruns. The constructor can be given an existing
    LDEAnalysis datafile to add more runs to"""
    
    
    def __init__(self):
        
        self.all_events=np.zeros((0,4),dtype=np.uint32)
        self.all_double_clicks=np.zeros((0,4),dtype=np.uint32)
        self.all_photon_hist_x=None
        self.all_photon_hist_w1_h0y=None
        self.all_photon_hist_w1_h1y=None
        self.all_photon_hist_w2_h0y=None
        self.all_photon_hist_w2_h1y=None
               
        self.all_window_data=[]
        self.all_gate_phases=[]
        self.all_ssro1=[]
        self.all_ssro2=[]
        self.all_CR1=[]
        self.all_CR2=[]
        
        self.all_tails=[]
        self.all_laser=[]
        self.all_statistics_lt1=[]
        self.all_statistics_lt2=[]
        self.all_statistics_plu=[]
        self.all_double_click_statistics=[]
        
        self.save_statistics_lt1=['get_noof_CR_checks', 'get_below_threshold_events']
        self.save_statistics_lt2=['get_noof_seq_starts', 'get_below_threshold_events', 
                                    'get_noof_CR_checks']
    
        self.total_corr=np.zeros(4,dtype=np.uint32)
        self.total_corr_00=np.zeros(4,dtype=np.uint32)
        self.total_corr_01=np.zeros(4,dtype=np.uint32)
        self.total_corr_10=np.zeros(4,dtype=np.uint32)
        self.total_corr_11=np.zeros(4,dtype=np.uint32)
        self.savename=''
        self.savedir=''
        self.all_runs=[]
        self.all_subruns=[]

     
    def analyse_lde_from_dir(self, startpath, **kw): 
        """This convenience method calls most of the interesting analysis 
        for an LDE experiment defined 'hhp-data-' files in a directory tree"""

        identifier = 'rawdata'
        runs=[]
        subruns=[]
        savedir=kw.pop('savedir',os.path.join(startpath,'analysis'))
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        for (path, dirs, files) in os.walk(startpath):
            for drn in dirs:
                if drn.find(identifier)==0: 
                    print drn
                    print path
                    idx=int(drn[-sscorr.DATIDXDIGITS:])
                    try:
                        run_nr = runs.index(path)
                        subruns[run_nr].append(idx)
                    except ValueError:
                        runs.append(path)
                        #run_nr=len(runs)-1
                        subruns.append([idx])
        print runs
        print subruns
        return self.analyse_lde(runs,subruns,savedir=savedir, **kw)
            
    def analyse_lde_from_runfile(self, runfile, **kw):
        """This convenience method calls most of the interesting analysis 
        for an LDE experiment defined by a runfile generated by the analysis deamon"""
        
        d=np.load(runfile)
        runs=d['runs']
        subruns=d['subruns']
        d.close()
        savedir=kw.pop('savedir',
                os.path.join(os.path.dirname(runfile),'analysis'))
        if not os.path.isdir(savedir):
            os.mkdir(savedir)
        return self.analyse_lde(runs,subruns,savedir=savedir, **kw)
    
    def reanalyse_lde(self, w_start = (637,666), 
            w_length=150, w_dt=-1, ch1_offset=-29, w_dt_min=-1, **kw):
        
        anal = SingleLDEAnalysis()

        anal.w1_start = w_start
        anal.w2_start = kw.get('w2_start', w_start)
        w1_length=w_length
        w2_length= kw.get('w2_length', w_length)
        anal.w1_stop = (anal.w1_start[0]+w1_length,anal.w1_start[1]+w1_length)
        anal.w2_stop = (anal.w2_start[0]+w2_length,anal.w2_start[1]+w2_length)
        anal.dt_max=w_dt #NOTE with this number you can adjust the time difference between the photons
        anal.dt_min=w_dt_min
        anal.ch1_offset=ch1_offset#-4
        
        filter_ap_length=kw.get('filter_ap', 0)
        anal.single_photon_range_start=(anal.w1_start[0]-filter_ap_length,anal.w1_start[1]-filter_ap_length)
        anal.single_photon_range_stop=anal.w1_stop 
        anal.hhp=self.all_events
        anal.ssro1=self.all_ssro1
        anal.ssro2=self.all_ssro2
        anal.correlations()
        apply_to_self=kw.pop('apply_to_self',True)
        if apply_to_self:
            self.total_corr,self.total_corr_00,self.total_corr_01,self.total_corr_10, self.total_corr_11 = \
                    anal.uncond_corr, anal.corr_00, anal.corr_01, anal.corr_10, anal.corr_11
        return anal.uncond_corr, anal.corr_00, anal.corr_01, anal.corr_10, anal.corr_11
        

    def analyse_lde(self, runs, subruns, savedir='', w_start = (637,666), 
            w_length=150, w_dt=-1, ch1_offset=-29, w_dt_min=-1, **kw):
        """
        This convenience method calls most of the interesting analysis 
        for an LDE experiment defined by an array of directories of LDE 'runs', 
        and a 2D-array with for each run the integers defining which subruns 
        to use. 
        The results are saved in the savedir, which defaults to the current 
        working directory.
        
        params: 
        - w_start(ch0,ch1) in bins
        - w_length in bins
        - w_dt max distance between clicks in bins
        - kw: * w2_start, w2_length = w_start, w_length
              * filter on afterpulsing, given as the length in bins before w_start
              to look for laser photons: filter_ap=0
              * 

        returns:
        SingleLDEAnalysis object of last subrun
        """
        do_analyse_g2=kw.get('analyse_g2',True)
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        self.savedir=savedir
        anal = SingleLDEAnalysis()

        anal.w1_start = w_start
        anal.w2_start = kw.get('w2_start', w_start)
        w1_length=w_length
        w2_length= kw.get('w2_length', w_length)
        anal.w1_stop = (anal.w1_start[0]+w1_length,anal.w1_start[1]+w1_length)
        anal.w2_stop = (anal.w2_start[0]+w2_length,anal.w2_start[1]+w2_length)
        anal.dt_max=w_dt #NOTE with this number you can adjust the time difference between the photons
        anal.dt_min=w_dt_min
        anal.ch1_offset=ch1_offset#-4
        
        filter_ap_length=kw.get('filter_ap', 0)
        anal.single_photon_range_start=(anal.w1_start[0]-filter_ap_length,anal.w1_start[1]-filter_ap_length)
        anal.single_photon_range_stop=anal.w1_stop
        

        if self.all_photon_hist_x == None:
            self.all_photon_hist_x=np.zeros(CHMAXTIME,dtype=np.uint32)
            self.all_photon_hist_w1_h0y=np.zeros(CHMAXTIME,dtype=np.uint32)
            self.all_photon_hist_w1_h1y=np.zeros(CHMAXTIME,dtype=np.uint32)
            self.all_photon_hist_w2_h0y=np.zeros(CHMAXTIME,dtype=np.uint32)
            self.all_photon_hist_w2_h1y=np.zeros(CHMAXTIME,dtype=np.uint32)


        self.savename = 'window1_ch0_%d-%d_ch1_%d-%d__window2_ch0_%d-%d_ch1_%d-%d-%d-dt_max-%d-dt_min-run1-4' % \
                    (anal.w1_start[0], anal.w1_stop[0], anal.w1_start[1], 
                            anal.w1_stop[1],anal.w2_start[0], anal.w2_stop[0], 
                            anal.w2_start[1], anal.w2_stop[1], anal.dt_max, 
                            anal.dt_min)
        
        
        
        self.g2_deltas = range(-2,3)
        self.g2_coincidences_all = [np.array([],dtype=int) \
                for d in self.g2_deltas]
        self.g2_coincidences_corr = [np.array([],dtype=int) \
                for d in self.g2_deltas]
        self.g2_coincidences_fullwindow = [np.array([],dtype=int) \
                for d in self.g2_deltas]
        self.empty_subruns=0
        for _i, datfolder in np.ndenumerate(runs):
            self.all_runs.append(datfolder)
            i=_i[0]
            self.all_subruns.append(subruns[i])
            for DATIDX in subruns[i]:
                DATIDXSTR=('-%.'+str(sscorr.DATIDXDIGITS)+'d') % DATIDX
                anal.datfolder = datfolder
                anal.rawfolder = os.path.join(datfolder, 'rawdata'+DATIDXSTR)
                hhpdata=os.path.join(datfolder, HHPFILEBASE+DATIDXSTR+'.npz')
                corrsavefile = os.path.join(datfolder, 
                        'corr_data_'+self.savename+DATIDXSTR)
                            
                ### data preparation or loading
                if not os.path.exists(hhpdata):
                    #print hhpdata
                    print 'run', datfolder, ', subrun', DATIDXSTR,\
                            'not yet prepared, starting preparation'
                    anal.prepare_data(DATIDX=DATIDX)  
                else:
                    print 'loading from hhp-data'
                    d=np.load(hhpdata)
                    anal.hhp=d['hhp']
                    anal.hhdata=d['hhdata']
                    anal.w1=d['w1']
                    anal.w2=d['w2']
                    d.close()
                if len(anal.hhdata)==0:
                    print 'run', datfolder, ', subrun', DATIDXSTR, 'EMPTY'
                    self.empty_subruns+=1
                    anal.hhdata=np.array([np.zeros(4,dtype=np.uint32)])
                    anal.hhp=np.array([np.zeros(4,dtype=np.uint32)])
                    anal.w1=np.array([np.zeros(4,dtype=np.uint32)])
                    anal.w2=np.array([np.zeros(4,dtype=np.uint32)])
                 
                ### adwin data 
                anal.get_adwin_data(DATIDX=DATIDX)
                
                ###statistics
              
                anal.get_adwin_statistics(DATIDX)            
                statistics_lt1 = [anal.adwin_lt1_stats[key] for key in self.save_statistics_lt1]
                statistics_lt2 = [anal.adwin_lt2_stats[key] for key in self.save_statistics_lt2]
                self.all_statistics_lt1.append(statistics_lt1)
                self.all_statistics_lt2.append(statistics_lt2)
                #print statistics_lt1
                
                ### g2
                if do_analyse_g2:
                    c_all = anal.get_g2_coincidences(self.g2_deltas, 
                            ch0window=(0,500), ch1window=(0,500))
                    c_corr = anal.get_g2_coincidences(self.g2_deltas)
                    c_fullwindow = anal.get_g2_coincidences(self.g2_deltas,
                            ch0window=(anal.w1_start[0],500), ch1window=(anal.w1_start[1],500))
                
                    self.g2_coincidences_all = tpqi.add_coincidences(
                            self.g2_coincidences_all, c_all)
                    self.g2_coincidences_corr = tpqi.add_coincidences(
                            self.g2_coincidences_corr, c_corr)
                    self.g2_coincidences_fullwindow = tpqi.add_coincidences(
                            self.g2_coincidences_fullwindow, c_fullwindow)
                
                ### the correlations            
                anal.correlations()
                self.total_corr+=anal.uncond_corr ##NOTE has to be brought into right order  
#                self.total_corr_psi1+=anal.corr_psi1
#                self.total_corr_psi2+=anal.corr_psi2
                self.total_corr_00+=anal.corr_00
                self.total_corr_01+=anal.corr_01
                self.total_corr_10+=anal.corr_10
                self.total_corr_11+=anal.corr_11

                ###save single analysis
                #anal.save(corrsavefile)

                if len(anal.w)>0:
                    if (len(anal.w[np.where(anal.w>0)])>0):
                        event_photons=anal.get_event_photons(anal.w,anal.hhp)
                       # print 'event photons', event_photons
                        self.all_events=np.vstack((self.all_events,event_photons))
                    x, w1_h0y, w1_h1y, w2_h0y, w2_h1y =\
                            self._get_time_histogram(anal.w1,anal.w2, CHMAXTIME,
                                    (0,CHMAXTIME-1))
                    self.all_photon_hist_x=x
                    self.all_photon_hist_w1_h0y+=w1_h0y
                    self.all_photon_hist_w1_h1y+=w1_h1y
                    self.all_photon_hist_w2_h0y+=w2_h0y
                    self.all_photon_hist_w2_h1y+=w2_h1y

                ###double_clicks check
                raw_double_clicks=sscorr.get_double_clicks(anal.hhdata,ch0_start=anal.w1_start[0], 
                        ch0_stop=anal.w1_start[0]+w1_length, \
                        ch1_start=anal.w1_start[1], 
                        ch1_stop=anal.w1_start[1]+w1_length)
                self.all_double_clicks=np.vstack((self.all_double_clicks,raw_double_clicks))
                double_click_statistics=sscorr.get_channel_statistics(raw_double_clicks)
                self.all_double_click_statistics.append(double_click_statistics)
                anal.print_plu_vs_hh()
                print 'The HH detected in total', len(raw_double_clicks)/2,\
                        'valid double click events'
                self.all_statistics_plu.append([anal.noof_hh_ssros, 
                    anal.noof_valid_ssros,len(raw_double_clicks)/2]) 
                ###tail_cts per run
                anal.tail_counts()
                self.all_tails = np.append(self.all_tails, anal.total_tail)
                anal.laser_counts()
                self.all_laser = np.append(self.all_laser, anal.total_laser)

                ### adding to all_arrays
                #NOTE when using indexing arrays, make sure that
                # a) the indexing array is not empty (or only consists only of False elements)
                # b) the array to be indexed is not empty either
                self.all_window_data =  np.append(self.all_window_data,anal.w[np.where(anal.w>0)] 
                      if len(anal.w) > 0 else [])
                self.all_gate_phases = np.append(self.all_gate_phases, 
                        anal.gate_phase[np.where(anal.w>0)] if len(anal.gate_phase) > 0 else [])
                self.all_ssro1 = np.append(self.all_ssro1,
                        anal.ssro1[np.where(anal.w>0)] if len(anal.ssro1) > 0 else [] )
                self.all_ssro2 = np.append(self.all_ssro2,
                        anal.ssro2[np.where(anal.w>0)] if len(anal.ssro2) > 0 else [] )
                self.all_CR1 = np.append(self.all_CR1,
                        anal.cr1[np.where(anal.w>0)] if len(anal.cr1) > 0 else [] )
                self.all_CR2 = np.append(self.all_CR2,
                        anal.cr2[np.where(anal.w>0)] if len(anal.cr2) > 0 else [] )
                
                
        print self.total_corr
        self.last_anal=anal
        return anal
        
    def save(self, filename=''):        
        if filename=='': 
            filename=os.path.join(self.savedir,'analysis_'+self.savename+'.pkl')
        #print filename
        f=open(filename,'wb')
        pickle.dump(self,f)
        f.close()


    def filter_on_gatephase(self, **kw):
        """returns the total correlations filtered on good gatephase. 
        kw: apply_to_self : boolean apply the filter to the total 
        correlations total_corr_ij"""
        if len(np.where(self.all_gate_phases>0)[0])==0:
            return [],[],[]
        gate_phases=self.all_gate_phases[self.all_gate_phases!=0]
        fltr = gate_phases>0
        
        return self.filter_correlations(fltr, **kw)

    def filter_on_CR(self, lt1_min, lt2_min, **kw):
        """returns the total correlations filtered on CR checks after RO. 
        Minimum counts during CR check for lt1 and lt2 should be given as arguments.
        kw: apply_to_self : boolean apply the filter to the total 
        correlations total_corr_ij"""

        if len(np.where(np.logical_and(self.all_CR1>lt1_min, self.all_CR2>lt2_min))[0])==0:
            return [],[],[]
        fltr = np.logical_and(self.all_CR1>lt1_min, self.all_CR2>lt2_min)
        return self.filter_correlations(fltr, **kw)
        
    def filter_correlations(self, filter, **kw):
        """returns the total correlations filtered on any Boolean array of length len(self.ssro1)
        kw: apply_to_self : boolean apply the filter to the total 
        correlations total_corr_ij"""
        apply_to_self=kw.pop('apply_to_self',True)
        if apply_to_self:
            self.total_corr,self.total_corr_00,self.total_corr_01,self.total_corr_10, self.total_corr_11 = \
            sscorr.correlations(self.all_ssro1[filter],
                self.all_ssro2[filter], 
                self.all_window_data[filter])
        return sscorr.correlations(self.all_ssro1[filter],
                self.all_ssro2[filter], 
                self.all_window_data[filter])

                
    def plot_correlations(self, F0LT2 = 0.805, F0LT1 = 0.905, F1LT2 = 0.998, F1LT1 = 0.9937,
            save_plots=True, **kw):

        """ Plots the correlations"""
        total_corr=kw.get('total_corr', self.total_corr) 
        total_corr_psi1=kw.get('total_corr_psi1', self.total_corr_00 + self.total_corr_11) 
        total_corr_psi2=kw.get('total_corr_psi2', self.total_corr_01 + self.total_corr_10)
        
        total_corr_err=sscorr.get_correlation_errors(total_corr)
        total_corr_psi1_err=sscorr.get_correlation_errors(total_corr_psi1)
        total_corr_psi2_err=sscorr.get_correlation_errors(total_corr_psi2)
      
        corrected_corr, corr_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(total_corr, 
                F0LT1, F0LT2, F1LT1, F1LT2, verbose = True, return_error_bars=True)
        corrected_corr_psi1, corr_psi1_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(total_corr_psi1, 
                F0LT1, F0LT2, F1LT1, F1LT2, verbose = True, return_error_bars=True)
        corrected_corr_psi2, corr_psi2_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(total_corr_psi2, 
                F0LT1, F0LT2, F1LT1, F1LT2, verbose = True, return_error_bars=True)

        sscorr.plot_uncorrected(total_corr, total_corr_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr'))
        sscorr.plot_uncorrected(total_corr_psi1, total_corr_psi1_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr_psi1'))
        sscorr.plot_uncorrected(total_corr_psi2, total_corr_psi2_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr_psi2'))
                 
        sscorr.plot_corrected(corrected_corr, corr_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr'))
        sscorr.plot_corrected(corrected_corr_psi1, corr_psi1_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr_psi1'))
        sscorr.plot_corrected(corrected_corr_psi2, corr_psi2_err,
                sav_fig=save_plots,save_path=\
                 os.path.join(self.savedir,self.savename+'total_corr_psi2'))

        ### g2 plots
        if save_plots:
            tpqi.plot_g2(self.g2_deltas, self.g2_coincidences_all,
                delta_separation_time=200, 
                savepath=os.path.join(self.savedir, 'g2_all.png'))
            tpqi.plot_g2(self.g2_deltas, self.g2_coincidences_corr,
                delta_separation_time=200, 
                savepath=os.path.join(self.savedir, 'g2_correlationwindow.png'))
            tpqi.plot_g2(self.g2_deltas, self.g2_coincidences_fullwindow,
                delta_separation_time=200, 
                savepath=os.path.join(self.savedir, 'g2_fulltail.png'))
                
    def get_run_subrun_from_datapoint_idx(self,idx):
        """returns the measurement run (directory string) and subrun index 
        corresponding to a total subrun number in the analysis"""
        number=int(idx)
        totsubs=0
        run=0
        for subruns in self.all_subruns:
            if totsubs+len(subruns)>=number:
                break
            run+=1
            totsubs+=len(subruns)
        #print 'run', run, 'number', number, 'totsubs', totsubs
        return list(cbook.flatten(self.all_runs))[run], list(cbook.flatten(self.all_subruns))[idx]
    
    def plot_tailcts(self,save_plots=True):
        """plots total subrun number versus tailcounts and laser counts"""
        fig=plt.figure()
        plt.title('Combined counts per seq-start*300 from both \
                channels+windows')
        ax=plt.subplot(211)
        ax.plot(np.arange(len(list(cbook.flatten(self.all_subruns)))),self.all_tails)
        ax.set_xlabel('Subrun total #')
        ax.set_ylabel('Tailcounts')
        ax=plt.subplot(212)
        ax.plot(np.arange(len(list(cbook.flatten(self.all_subruns)))),self.all_laser)
        ax.set_ylabel('Lasercounts')
        
        if save_plots:
            fig.savefig(os.path.join(self.savedir,self.savename)+'tailcounts.png')
        
    def plot_statistics(self,setup='lt2',save_plots=True,**kw):
        """plots total subrun number versus adwin statistics saved 
        from in the lde_analysis
        pars:
        - setup: h=give setup to plot statistics for
        - known kw: plot_statistics: list of statistics to plot, 
                    defaults to all saved satistics for given setup.
                    note: statistics in plot_statistics must be 
                    contained in self.save_statistics for the given setup.
        """
        if setup=='lt1':
            stats=self.all_statistics_lt1
            saved_stats=self.save_statistics_lt1
            plot_statistics=kw.get('plot_statistics',self.save_statistics_lt1)
        elif setup=='lt2':
            stats=self.all_statistics_lt2
            saved_stats=self.save_statistics_lt2
            plot_statistics=kw.get('plot_statistics',self.save_statistics_lt2)
        else: 
            print 'unknown setup'
            return
        num_stat=len(plot_statistics)
        i=0
        fig=plt.figure()
        for stat in plot_statistics:
            try:
                stat_column=saved_stats.index(stat)
            except ValueError:
                print 'statistic', stat, 'not found in analysed statistcs'
                continue
            i+=1
            ax=plt.subplot(num_stat,1,i)
            ax.plot(np.arange(len(list(cbook.flatten(self.all_subruns)))),
                    [subrun[stat_column] for subrun in stats])
            ax.set_xlabel('Subrun total #')
            ax.set_ylabel(stat)
        if save_plots:
            fig.savefig(os.path.join(self.savedir,self.savename)+'statistics_'+setup+'.png')
    
    def plot_plu_vs_tail_all(self,range=None,bins=None,bins_plu=30,**kw):
        """Plots a histogram of the arrival times of the gated plu events vs the
        raw tail data of the last subrun
        keyword args:   range=last_anal.w1_(start,stop), bins=full_range,
                        bins_plu=30, log_plots=False, save_plots=False, save_path=''   
        """
        w1_photons_plu, w2_photons_plu =\
                        self._get_double_click_windows(self.all_events,2)
        last_anal=self.last_anal
        if range == None:
            range=(last_anal.w1_start[0], last_anal.w1_stop[0])
        if bins == None:
            bins = range[1]-range[0]
        #print 'length x before:', len(self.all_photon_hist_x)
        #print 'length x:', len(self.all_photon_hist_x[range[0]:range[1]])
        #print 'range:', range
        #print 'length y:', len(self.all_photon_hist_w1_h1y[range[0]:range[1]])
        self._plot_windowed_photons(w1_photons_plu, w2_photons_plu,
                    bins=bins_plu,
                    range=range,
                    x2=self.all_photon_hist_x[range[0]:range[1]+1],
                    w1_h0y_2=self.all_photon_hist_w1_h0y[range[0]:range[1]], 
                    w1_h1y_2=self.all_photon_hist_w1_h1y[range[0]:range[1]], 
                    w2_h0y_2=self.all_photon_hist_w2_h0y[range[0]:range[1]], 
                    w2_h1y_2=self.all_photon_hist_w2_h1y[range[0]:range[1]],
                    c1='r', c2='g',
                     **kw)

    def plot_plu_vs_tail_last_subrun(self,**kw):
        """Plots a histogram of the arrival times of the gated plu events vs the
        raw tail data of the last subrun
        keyword args:   range=last_anal.w1_(start,stop),bins=full_range,
                        bins_plu=30, log_plots=False, save_plots=False, save_path=''   
        """
        w1_photons_plu, w2_photons_plu =\
                        self._get_double_click_windows(self.all_events,2)
        last_anal=self.last_anal
        range=kw.pop('range',(last_anal.w1_start[0], last_anal.w1_stop[0]))
        bins = kw.pop('bins',range[1]-range[0])
        bins_plu= kw.pop('bins_plu',30)

        self._plot_windowed_photons(last_anal.w1, last_anal.w2, 
                bins=bins,
                range=range,
                w1_photons_2=w1_photons_plu, w2_photons_2=w2_photons_plu, 
                bins_2=bins_plu, **kw)

    def plot_plu_vs_hh_events(self,**kw):
        """Plots a histogram of the arrival times of the gated plu events vs the
        double click events as seen by the hydraharp
        keyword args: - all_plu_events, all_hh_clicks 
                        (default to self.all_events, self.all_double_clicks)
                      - bins=700, range=(0,699), bins_2=700, 
                        log_plots=False, save_plots=False, save_path=''   
        """
        all_plu_events=kw.pop('all_plu_events',self.all_events)
        all_hh_clicks=kw.pop('all_hh_clicks',self.all_double_clicks)
        
        w1_photons_plu, w2_photons_plu =\
                self._get_double_click_windows(all_plu_events,2)
        w1_photons_hh, w2_photons_hh =\
                self._get_double_click_windows(all_hh_clicks)

        self._plot_windowed_photons(w1_photons_hh,w2_photons_hh,
                                   w1_photons_2=w1_photons_plu,
                                   w2_photons_2=w2_photons_plu,
                                   **kw)

    def plot_plu_good_vs_bad_events(self,good_ssro=[1,2],**kw):
        """Plots a histogram of the arrival times of the gated plu events,
        with one axis good events, and on the other the bad events
        as defined by the good_ssro as a list of integers 0-3 [00,01,10,11]
        - **kw:  bins=700, range=(0,699), bins_2=700, 
                        log_plots=False, save_plots=False, save_path=''
        """
        
        all_plu_events=kw.pop('all_plu_events',self.all_events)
        all_ssro1=kw.pop('all_ssro1', self.all_ssro1)
        all_ssro2=kw.pop('all_ssro2', self.all_ssro2)


        good_events, bad_events=\
                self.get_good_bad_events(all_plu_events,all_ssro1,all_ssro2, good_ssro)
        print 'good:', len(good_events)/3, 'bad:', len(bad_events)/3
        w1_photons_good, w2_photons_good = \
                self._get_double_click_windows(good_events,2)
        w1_photons_bad, w2_photons_bad = \
                self._get_double_click_windows(bad_events, 2)

        self._plot_windowed_photons(w1_photons_good,w2_photons_good,
                                   w1_photons_2=w1_photons_bad,
                                   w2_photons_2=w2_photons_bad,
                                   **kw)

    def get_good_bad_events(self,all_plu_events,all_ssro1,all_ssro2,good_ssro):
        """Divides hhplu data into good and bad by the ssro correlation they have
        as defined by the good_ssro as a list of integers 0-3 [00,01,10,11]"""
        
        bad_events= np.zeros((0,4),dtype=np.uint32)
        good_events= np.zeros((0,4),dtype=np.uint32)
        for i in range(len(all_ssro1)):
            corr=sscorr.correlations(np.array([all_ssro1[i]]),
                    np.array([all_ssro2[i]]),np.array([1]))[0]
            if (1 in [corr[good] for good in good_ssro]):
                good_events=np.vstack((good_events,all_plu_events[3*i]))
                good_events=np.vstack((good_events,all_plu_events[3*i+1]))
                good_events=np.vstack((good_events,all_plu_events[3*i+2]))
            else:
                bad_events=np.vstack((bad_events,all_plu_events[3*i]))
                bad_events=np.vstack((bad_events,all_plu_events[3*i+1]))
                bad_events=np.vstack((bad_events,all_plu_events[3*i+2]))

        return good_events, bad_events


    def plot_dts(self,hrange=(-200,200), bins=700, good_ssro=[1,2], save_fig=False, save_path='', **kw):
        """Plots a histogram of the dt's (w2 arrival time - w1 arrival time) 
        of the gated plu events, with one axis good events, and on the other 
        the bad events as defined by the good_ssro as a list of 
        integers 0-3 [00,01,10,11]"""
        
        all_plu_events=kw.pop('all_plu_events',self.all_events)
        all_ssro1=kw.pop('all_ssro1', self.all_ssro1)
        all_ssro2=kw.pop('all_ssro2', self.all_ssro2)

        good_events, bad_events=self.get_good_bad_events(all_plu_events,all_ssro1,all_ssro2, good_ssro)
        dts_good=np.array([])
        for i in range(int(len(good_events)/3)):
            dt=int(good_events[3*i+1][1])-int(good_events[3*i][1])
            dts_good=np.append(dts_good,dt)
        
        dts_bad=np.array([])
        for i in range(int(len(bad_events)/3)):
            dt=int(bad_events[3*i+1][1])-int(bad_events[3*i][1])
            dts_bad=np.append(dts_bad,dt)


        
        plt.figure()
        ax = plt.subplot(111)
        plt.title('dts for all plu events')
        plt.hist(dts_good,range=hrange, bins=bins, facecolor='none', edgecolor= 'g' ,
                    hatch = '\\')
        ax.set_xlabel('Bins')
        ax2=ax.twinx()
        plt.hist(dts_bad,range=hrange, bins=bins, facecolor='none', edgecolor= 'r' ,
                    hatch = '//')

        if save_fig:
            fig1.savefig(save_path+'photons_window1.png')

        

    def _get_double_click_windows(self, double_clicks, delete_marker=None):
        if delete_marker!=None:
            double_clicks=hht3.delete_markers(double_clicks,delete_marker)
        w1_photons=np.zeros((0,4),dtype=np.uint32)
        w2_photons=np.zeros((0,4),dtype=np.uint32)

        for i in np.arange(len(double_clicks)):
            if (i%2)==0:
                w1_photons=np.vstack((w1_photons,double_clicks[i]))
            if (i%2)==1:
                w2_photons=np.vstack((w2_photons,double_clicks[i]))
        
        return w1_photons, w2_photons 

    def _get_time_histogram(self,w1_photons,w2_photons, bins, range):
        
        w1_ch0,w1_ch1=hht3.get_click_times(w1_photons)
        w2_ch0,w2_ch1=hht3.get_click_times(w2_photons)

        #print 'w1_ch0',len(w1_ch0),'w1_ch1',len(w1_ch1)
        #print  'w2_ch0',len(w2_ch0), 'w2_ch0',len(w2_ch1)

        w1_h0y, _tmp = np.histogram(w1_ch0, bins=bins, range=range)
        w1_h1y, _tmp = np.histogram(w1_ch1, bins=bins, range=range)

        w2_h0y, _tmp = np.histogram(w2_ch0, bins=bins, range=range)
        w2_h1y, x = np.histogram(w2_ch1, bins=bins, range=range)

        return  x, w1_h0y, w1_h1y, w2_h0y, w2_h1y

    def _plot_windowed_photons(self,w1_photons,w2_photons,bins=700,range=(0,699),
            w1_photons_2=None, w2_photons_2=None,bins_2=700, **kw):
        
        two_axis_plot= w1_photons_2!=None
       
        print 'number of photons in Window 1', len(w1_photons)
        print 'number of photons in Window 2', len(w2_photons)
        #print 'number of w1 ch0' , len(w1_photons[np.where(w1_photons[:,2]==0)])
        #print 'number of w1 ch1' , len(w1_photons[np.where(w1_photons[:,2]==1)])
        #print 'number of w2 ch0' , len(w2_photons[np.where(w2_photons[:,2]==0)])
        #print 'number of w2 ch1' , len(w2_photons[np.where(w2_photons[:,2]==1)]) 
        if(two_axis_plot):
            print '2 number of photons in Window 1', len(w1_photons_2)
            print '2 number of photons in Window 2', len(w2_photons_2)
            #print '2 number of w1 ch0' , len(w1_photons_2[np.where(w1_photons_2[:,2]==0)])
            #print '2 number of w1 ch1' , len(w1_photons_2[np.where(w1_photons_2[:,2]==1)])
            #print '2 number of w2 ch0' , len(w2_photons_2[np.where(w2_photons_2[:,2]==0)])
            #print '2 number of w2 ch1' , len(w2_photons_2[np.where(w2_photons_2[:,2]==1)]) 

        x, w1_h0y, w1_h1y, w2_h0y, w2_h1y = \
                self._get_time_histogram(w1_photons,w2_photons, bins, range)
        
        if(two_axis_plot):
            x2, w1_h0y_2, w1_h1y_2, w2_h0y_2, w2_h1y_2 =\
                    self._get_time_histogram(w1_photons_2,w2_photons_2, bins_2, range)
            self._plot_histograms(x, w1_h0y, w1_h1y, w2_h0y, w2_h1y, 
                       x2, w1_h0y_2, w1_h1y_2, w2_h0y_2, w2_h1y_2, **kw)
        else:
            self._plot_histograms(x, w1_h0y, w1_h1y, w2_h0y, w2_h1y, **kw)

    def _plot_histograms(self, x, w1_h0y, w1_h1y, w2_h0y, w2_h1y, 
                       x2=None, w1_h0y_2=None, w1_h1y_2=None, w2_h0y_2=None, w2_h1y_2=None,
                       c1='g', c2= 'r',log_plots=False,
                       save_fig=False,save_path=''):
        """This plots 2 figures x 2 subplots x 2 y-axes in barplots, 
        as given by w1,w2_h0,h1_(x,y)_1,2 variables."""
        two_axis_plot= x2!=None
        if save_path=='':
            save_path=self.savedir
    
        fig1 = plt.figure()  

        ax = plt.subplot(121)
        plt.title('window 1, channel 0')
        plot_bar(x[:-1], w1_h0y, width=(max(x)-min(x))/len(x),  facecolor='none', 
                edgecolor= c1,label='channel 0', hatch = '/', log=log_plots)
        ax.set_xlabel('Bins')
        if two_axis_plot:
            ax2=ax.twinx()
            plot_bar(x2[:-1], w1_h0y_2, width=(max(x2)-min(x2))/len(x2),  
                    facecolor='none', edgecolor= c2 ,label='channel 0_2',
                    hatch = '\\', log=log_plots)

        ax = plt.subplot(122)
        plt.title('window 1, channel 1')
        plot_bar(x[:-1], w1_h1y, width=(max(x)-min(x))/len(x),  facecolor='none', 
                edgecolor= c1 ,label='channel 1', hatch = '/', log=log_plots)
        ax.set_xlabel('Bins')
        if two_axis_plot:
            ax2=ax.twinx()
            plot_bar(x2[:-1], w1_h1y_2, width=(max(x2)-min(x2))/len(x2),  
                    facecolor='none', edgecolor= c2, label='channel 1_2', 
                    hatch = '\\', log=log_plots)

        if save_fig:
            fig1.savefig(os.path.join(save_path,'photons_window1.png'))


        fig2 = plt.figure()

        ax = plt.subplot(121)
        plt.title('window 2, channel 0')
        plot_bar(x[:-1], w2_h0y, width=(max(x)-min(x))/len(x),  facecolor='none', 
                edgecolor= c1,label='channel 0', hatch = '/', log=log_plots)
        ax.set_xlabel('Bins')
        if two_axis_plot:
            ax2=ax.twinx()
            plot_bar(x2[:-1], w2_h0y_2, width=(max(x2)-min(x2))/len(x2),  
                    facecolor='none', edgecolor= c2 ,label='channel 0_2', 
                    hatch = '\\', log=log_plots)

        ax = plt.subplot(122)
        plt.title('window 2, channel 1')
        plot_bar(x[:-1], w2_h1y, width=(max(x)-min(x))/len(x),  facecolor='none', 
                edgecolor= c1 ,label='channel 1', hatch = '/', log=log_plots)
        ax.set_xlabel('Bins')
        if two_axis_plot:
            ax2=ax.twinx()
            plot_bar(x2[:-1], w2_h1y_2, width=(max(x2)-min(x2))/len(x2), 
                    facecolor='none', edgecolor = c2,label='channel 1_2',
                    hatch = '\\', log=log_plots)

        if save_fig:
            fig2.savefig(os.path.join(save_path,'photons_window2.png'))
       
def plot_bar(x,y,**kw):
    if len(x)>0 and len(y)>0 and np.sum(y)>0:
        #print len(x), len(y), np.sum(y)
        plt.bar(x,y,**kw)
 
def load_previous_analysis(filename):
    f=open(filename,'rb')
    a = pickle.load(f)
    f.close()
    return a
       


""" Below are some typical uses of the above class."""

if __name__ == '__main__':
    
    ###option 1: The measurement is defined by a folder, 
                #we find all hhpdata-files in the folder an analyse them
    datadir=r'D:\Analysis\2012-09_ldetesting\autodata\20120905'
    a = LDEAnalysis()
    a.analyse_lde_from_dir(datadir, w_start = (234,229), w_length=150, w_dt=-1)

    ###option 2, measurement defined by runs, subruns:
    run1=r'D:\Analysis\2012-09_ldetesting\autodata\20120905\013113_LDE_Entanglement_ZZ'
    run2=r'D:\Analysis\2012-09_ldetesting\autodata\20120905\013113_LDE_Entanglement_ZZ'
    runs = [run1,run2]
    subruns = [np.arange(14),[1,3,4]]
    a = LDEAnalysis()
    a.analyse_lde(runs, subruns,savedir=r'D:\Analysis\2012-09_lde',
                    w_start = (234,229), w_length=150, w_dt=-1)

    #we can add some more data to the existing anaysis:
    another_run=r'D:\Analysis\2012-09_lde\2012-09-05\093113_LDE_Entanglement_ZZ'
    subruns=[np.arange(4)]
    a.analyse_lde([another_run], subruns,savedir=r'D:\Analysis\2012-09_lde',
                    w_start = (234,229), w_length=150, w_dt=-1)
    
    ### save the current analysis (saves all_photons, all_total_corrs, all gate_phases etc) 
    a.save(r'D:\Analysis\test_lde_analysis.npz')
    
    ### we can also work/add to an existing previously saved_analysis: 
    #a = load_previous_analysis('D:\Analysis\test_lde_analysis.npz')
    
    ###Plot the results and save the images
    a.plot(F0LT2 = 0.805, F0LT1 = 0.905, F1LT2 = 0.998, F1LT1 = 0.9937)
    
    ### do some more filtering
    gate_filtered_corr, gate_filtered_corr_psi1, gate_filtered_corr_psi2 = \
            a.filter_on_gatephase()
    cr_filtered_corr, cr_filtered_corr_psi1, cr_filtered_corr_psi2 = \
           a.filter_on_CR(40,10)
    
    ### plot the filtered total correlations, and leave out the psi-corrs.
    a.plot(total_corr=gate_filtered_corr,
    F0LT2 = 0.805, F0LT1 = 0.905, F1LT2 = 0.998, F1LT1 = 0.9937,
    xx_measurement=False)

    ###plot some more stuff:
    a.plot_tailcts()
    a.plot_statistics(setup='lt2')
    a.plot_statistics(setup='lt1')
    
    print'All:', a.total_corr, '\n 00:', a.total_corr_00,'\n 01:', a.total_corr_01, \
            '\n 10:', a.total_corr_10,'\n 11:', a.total_corr_11
    
    a.filter_correlations(a.all_laser<.015)
    a.total_corr,a.total_corr_00,a.total_corr_01,a.total_corr_10,a.total_corr_11=a.filter_on_gatephase()
    a.total_corr,a.total_corr_00,a.total_corr_01,a.total_corr_10,a.total_corr_11=a.filter_on_CR(10,1)
    a.plot_correlations()
    ###You can get your analysis object form a saved p[revious analysis
    s= load(r'D:\Analysis\test_lde_analysis.npz')






>>>>>>> 998bfb9d754ee59f58c639b71f91f58a0a5b6921
