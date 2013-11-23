from analysis.lib.lde import sscorr, lde_analysis
import os
import numpy as np
import matplotlib.pyplot as plt
import pickle


def ro_c_F_S(corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=True):
    if ro_correct:
        roc=sscorr.ssro_correct_twoqubit_state_photon_numbers(corr, F0a, F0b, F1a, F1b, 
                        dF0a=dF0a, dF0b=dF0b, dF1a=dF1a, dF1b=dF1b, verbose=False)
        val=roc[0][0]*roc[3][0]
        if val > 0:
            return np.sqrt(val)
        else:
            return 0.0
    else:
        N=float(corr.sum())
        val=(corr[0]/N*corr[3]/N)
        if val > 0:
            return np.sqrt(val)
        else:
            return 0.0

def ro_c_F_even(corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=True):
    if ro_correct:
        roc=sscorr.ssro_correct_twoqubit_state_photon_numbers(corr, F0a, F0b, F1a, F1b, 
                        dF0a=dF0a, dF0b=dF0b, dF1a=dF1a, dF1b=dF1b, verbose=False)
        return (roc[0][0]+roc[3][0])
    else:
        N=float(corr.sum())
        return (corr[0]/N+corr[3]/N)
    
def ro_c_F_odd(corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=True):
    if ro_correct:
        roc=sscorr.ssro_correct_twoqubit_state_photon_numbers(corr, F0a, F0b, F1a, F1b, 
                        dF0a=dF0a, dF0b=dF0b, dF1a=dF1a, dF1b=dF1b, verbose=False)
        return (roc[1][0]+roc[2][0])
    else:
        N=float(corr.sum())
        return (corr[1]/N+corr[2]/N)

def ro_c_dF(corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=True):
    if ro_correct:
        return sscorr.get_fidelity_error(corr, F0a, F0b, F1a, F1b, 
                        dF0a=dF0a, dF0b=dF0b, dF1a=dF1a, dF1b=dF1b)
    else:
        N11,N10,N01,N00=[float(i) for i in corr]
        return np.sqrt((N00+N11)*(N01+N10)/(N00+N01+N10+N11)**3)
    
def ro_c_dF_S(corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=True):
    if ro_correct:
        return sscorr.get_fidelity_error_sqrt_ZZ(corr, F0a, F0b, F1a, F1b, 
                        dF0a=dF0a, dF0b=dF0b, dF1a=dF1a, dF1b=dF1b)
    else:
        N11,N10,N01,N00=[float(i) for i in corr]
        val=1/4.*((N00*(N00+N01+N10)*(N00+N01+N10-N11)**2+4*N00*N10*N11*(N00+N01+N11)+\
                   4*N00*N01*N11*(N00+N10+N11)+N11*(N01+N10+N11)*(-N00+N01+N10+N11)**2)/(N00+N01+N10+N11)**5)
        if val > 0:
            return np.sqrt(val)
        else:
            return 0.0
            
def get_fidelity(ZZ_corr,XX_corr,XmX_corr,ro_correct=True, psi1=True, F0a=0.905, F0b=0.805, 
                    F1a=0.9937, F1b=0.998, dF0a=0.01, dF0b=0.01, dF1a=0.01, dF1b=0.01):

    ZZ=   ro_c_F_odd(ZZ_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    dZZ=     ro_c_dF(ZZ_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    ZZS=    ro_c_F_S(ZZ_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    dZZS=  ro_c_dF_S(ZZ_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    XX=  ro_c_F_even(XX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)\
                        if psi1 else ro_c_F_odd(XX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    dXX=     ro_c_dF(XX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    XmX= ro_c_F_odd(XmX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)\
                        if psi1 else ro_c_F_even(XmX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    dXmX=   ro_c_dF(XmX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    
    XXavg = XmX#(XX/dXX**2 + XmX/dXmX**2)/(1/dXX**2+1/dXmX**2)
    dXXavg= dXmX#np.sqrt(1/(1/dXX**2+1/dXmX**2))
    
    F=ZZ/2. - ZZS + (XXavg - 1/2.)
    dF= np.sqrt(1/4.*dZZ**2 + dZZS**2 + dXXavg**2)
    
    return F,dF,{'ZZ': ZZ, 'dZZ': dZZ, 'ZZS': ZZS, 'dZZS': dZZS, 
                    'XX': XX, 'dXX': dXX, 'XmX': XmX, 'dXmX':dXmX, 
                    'XXavg':XXavg,'dXXavg':dXXavg}
    
class FidelityAnalysis:

    def __init__(self,**kw):
        self.basepath=kw.pop('basepath',r'D:\analysis\data\lde')
        self.folder_ZZ=os.path.join(self.basepath,'ZZ')
        self.folder_XX=os.path.join(self.basepath,'XX')
        self.folder_XmX=os.path.join(self.basepath,'X-X')
        self.F0b = 0.805
        self.F0a = 0.905 
        self.F1b = 0.998
        self.F1a = 0.9937
        
        self.dF0b = 0.01
        self.dF0a = 0.01 
        self.dF1b = 0.01
        self.dF1a = 0.01

    
    def analyse_fidelity(self,state,ro_correct=True,w_start = (637,666), w_length=270, w_dt=-1, **kw):
        
        self.ZZ_an=lde_analysis.LDEAnalysis()
        self.XX_an=lde_analysis.LDEAnalysis()
        self.XmX_an=lde_analysis.LDEAnalysis()

        self.ZZ_an.analyse_lde_from_dir(self.folder_ZZ, w_start = w_start, 
        w_length=w_length, w_dt=w_dt, analyse_g2=False, **kw)
        self.ZZ_an.filter_on_gatephase()
        self.XX_an.analyse_lde_from_dir(self.folder_XX, w_start = w_start, 
        w_length=w_length, w_dt=w_dt, analyse_g2=False, **kw)
        self.XX_an.filter_on_gatephase()
        self.XmX_an.analyse_lde_from_dir(self.folder_XmX, w_start = w_start, 
        w_length=w_length, w_dt=w_dt, analyse_g2=False, **kw)
        self.XmX_an.filter_on_gatephase()      
        
        self.state=state
        self.w_start=w_start
        self.w_length=w_length
        self.w_dt=w_dt
        self.state=state
        if state=='psi1':
            psi1=True
            self.ZZ_corr =   self.ZZ_an.total_corr_00 +  self.ZZ_an.total_corr_11
            self.XX_corr =   self.XX_an.total_corr_00 +  self.XX_an.total_corr_11
            self.XmX_corr = self.XmX_an.total_corr_00 + self.XmX_an.total_corr_11
        elif state == 'psi2':
            psi1=False
            self.ZZ_corr=   self.ZZ_an.total_corr_01 +  self.ZZ_an.total_corr_10
            self.XX_corr=   self.XX_an.total_corr_01 +  self.XX_an.total_corr_10
            self.XmX_corr= self.XmX_an.total_corr_01 + self.XmX_an.total_corr_10       
        else:
            raise(Exception('Unknown state' + state))

        self.F, self.dF,self.F_dict = \
                get_fidelity(self.ZZ_corr,self.XX_corr,self.XmX_corr,
                            F0a=self.F0a, F0b=self.F0b, F1a=self.F1a, F1b=self.F1b, 
                            dF0a=self.dF0a, dF0b=self.dF0b, dF1a=self.dF1a, dF1b=self.dF1b,
                            ro_correct=ro_correct, psi1=psi1)
                
        self.dZZ=1/4.*self.F_dict['dZZ']**2 + self.F_dict['dZZS']**2 
        self.dXX=self.F_dict['dXXavg']**2
            
        self.F2stdev = (self.F-0.5)/self.dF

        self.print_fidelity()
       

    def reanalyse_fidelity(self,state,ro_correct=True, w_start = (637,666), w_length=150, w_dt=50, **kw):
    
        ZZ_total,ZZ_00,ZZ_01,ZZ_10,ZZ_11 = self.ZZ_an.reanalyse_lde(w_start = w_start, 
        w_length=w_length, w_dt=w_dt, apply_to_self=False, **kw)
        
        XX_total,XX_00,XX_01,XX_10,XX_11 = self.XX_an.reanalyse_lde(w_start = w_start, 
        w_length=w_length, w_dt=w_dt, apply_to_self=False, **kw)
        
        XmX_total,XmX_00,XmX_01,XmX_10,XmX_11 = self.XmX_an.reanalyse_lde(w_start = w_start, 
        w_length=w_length, w_dt=w_dt, apply_to_self=False, **kw)
        
        self.state=state
        self.w_start=w_start
        self.w_length=w_length
        self.w_dt=w_dt
        
        if state=='psi1':
            psi1=True
            self.ZZ_corr =   ZZ_00 +  ZZ_11
            self.XX_corr =   XX_00 +  XX_11
            self.XmX_corr = XmX_00 + XmX_11
        elif state == 'psi2':
            psi1=False
            self.ZZ_corr=   ZZ_01 +  ZZ_10
            self.XX_corr=   XX_01 +  XX_10
            self.XmX_corr= XmX_01 + XmX_10       
        else:
            raise(Exception('Unknown state' + state))

        self.F, self.dF,self.F_dict = \
                get_fidelity(self.ZZ_corr,self.XX_corr,self.XmX_corr,
                            F0a=self.F0a, F0b=self.F0b, F1a=self.F1a, F1b=self.F1b, 
                            dF0a=self.dF0a, dF0b=self.dF0b, dF1a=self.dF1a, dF1b=self.dF1b,
                            ro_correct=ro_correct, psi1=psi1)
                
        self.dZZ=1/4.*self.F_dict['dZZ']**2 + self.F_dict['dZZS']**2 
        self.dXX=self.F_dict['dXXavg']**2
            
        self.F2stdev = (self.F-0.5)/self.dF
        
        self.print_fidelity()

    def print_fidelity(self):
        
        print 'F:',self.F
        print 'sigma_F:', self.dF
        print 'total N:', self.ZZ_corr.sum()+self.XX_corr.sum()+self.XmX_corr.sum()

        print 'sigma_ZZ_contrib:', self.dZZ
        print 'sigma_XX_contrib:', self.dXX

        print '(F-0.5)/sigma:', self.F2stdev
    
    def plot(self,save_path=r'D:\analysis\output'):
        
        ZZ_corr_err=sscorr.get_correlation_errors(self.ZZ_corr)
        XX_corr_err=sscorr.get_correlation_errors(self.XX_corr)
        XmX_corr_err=sscorr.get_correlation_errors(self.XmX_corr)
        
        fig=plt.figure()
        plt.subplot(231)
        sscorr.plot_uncorrected(self.ZZ_corr, ZZ_corr_err)
        plt.subplot(232)
        sscorr.plot_uncorrected(self.XX_corr, XX_corr_err)
        plt.subplot(233)
        sscorr.plot_uncorrected(self.XmX_corr, XmX_corr_err)
        
        ZZ_c_corr, ZZ_c_corr_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(self.ZZ_corr, 
                        self.F0a, self.F0b, self.F1a, self.F1b, 
                        dF0a=self.dF0a, dF0b=self.dF0b, dF1a=self.dF1a, dF1b=self.dF1b,
                        verbose=False, return_error_bars=True)
        XX_c_corr, XX_c_corr_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(self.XX_corr, 
                        self.F0a, self.F0b, self.F1a, self.F1b, 
                        dF0a=self.dF0a, dF0b=self.dF0b, dF1a=self.dF1a, dF1b=self.dF1b,
                        verbose=False, return_error_bars=True)
        XmX_c_corr, XmX_c_corr_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(self.XmX_corr, 
                        self.F0a, self.F0b, self.F1a, self.F1b, 
                        dF0a=self.dF0a, dF0b=self.dF0b, dF1a=self.dF1a, dF1b=self.dF1b,
                        verbose=False, return_error_bars=True)
        plt.subplot(234)
        sscorr.plot_corrected(ZZ_c_corr, ZZ_c_corr_err)
        plt.subplot(235)
        sscorr.plot_corrected(XX_c_corr, XX_c_corr_err)
        plt.subplot(236)
        sscorr.plot_corrected(XmX_c_corr, XmX_c_corr_err)
        
        plt.suptitle('Correlations for state '+self.state+' with w_start ' + str(self.w_start) +\
                        ', w_length ' +str(self.w_length)+ ', w_dt ' + str(self.w_dt) + \
                        '\n and Fidelity F='+str(self.F)+'+/-'+str(self.dF))
        fig.set_size_inches(12,12)
        fig.savefig(os.path.join(save_path,'fidelity'+'_'+self.state+'.pdf'))
        return fig
    
    def save(self, save_path=r'D:\analysis\output'):        
        filename=os.path.join(save_path,'fidelity_analysis_'+self.state+'.pkl')
        #print filename
        f=open(filename,'wb')
        pickle.dump(self,f)
        f.close()
        
def load_previous_analysis(filename):
    f=open(filename,'rb')
    a = pickle.load(f)
    f.close()
    return a
>>>>>>> 998bfb9d754ee59f58c639b71f91f58a0a5b6921
=======
from analysis.lib.lde import sscorr, lde_analysis
import os
import numpy as np
import matplotlib.pyplot as plt
import pickle


def ro_c_F_S(corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=True):
    if ro_correct:
        roc=sscorr.ssro_correct_twoqubit_state_photon_numbers(corr, F0a, F0b, F1a, F1b, 
                        dF0a=dF0a, dF0b=dF0b, dF1a=dF1a, dF1b=dF1b, verbose=False)
        val=roc[0][0]*roc[3][0]
        if val > 0:
            return np.sqrt(val)
        else:
            return 0.0
    else:
        N=float(corr.sum())
        val=(corr[0]/N*corr[3]/N)
        if val > 0:
            return np.sqrt(val)
        else:
            return 0.0

def ro_c_F_even(corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=True):
    if ro_correct:
        roc=sscorr.ssro_correct_twoqubit_state_photon_numbers(corr, F0a, F0b, F1a, F1b, 
                        dF0a=dF0a, dF0b=dF0b, dF1a=dF1a, dF1b=dF1b, verbose=False)
        return (roc[0][0]+roc[3][0])
    else:
        N=float(corr.sum())
        return (corr[0]/N+corr[3]/N)
    
def ro_c_F_odd(corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=True):
    if ro_correct:
        roc=sscorr.ssro_correct_twoqubit_state_photon_numbers(corr, F0a, F0b, F1a, F1b, 
                        dF0a=dF0a, dF0b=dF0b, dF1a=dF1a, dF1b=dF1b, verbose=False)
        return (roc[1][0]+roc[2][0])
    else:
        N=float(corr.sum())
        return (corr[1]/N+corr[2]/N)

def ro_c_dF(corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=True):
    if ro_correct:
        return sscorr.get_fidelity_error(corr, F0a, F0b, F1a, F1b, 
                        dF0a=dF0a, dF0b=dF0b, dF1a=dF1a, dF1b=dF1b)
    else:
        N11,N10,N01,N00=[float(i) for i in corr]
        return np.sqrt((N00+N11)*(N01+N10)/(N00+N01+N10+N11)**3)
    
def ro_c_dF_S(corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=True):
    if ro_correct:
        return sscorr.get_fidelity_error_sqrt_ZZ(corr, F0a, F0b, F1a, F1b, 
                        dF0a=dF0a, dF0b=dF0b, dF1a=dF1a, dF1b=dF1b)
    else:
        N11,N10,N01,N00=[float(i) for i in corr]
        val=1/4.*((N00*(N00+N01+N10)*(N00+N01+N10-N11)**2+4*N00*N10*N11*(N00+N01+N11)+\
                   4*N00*N01*N11*(N00+N10+N11)+N11*(N01+N10+N11)*(-N00+N01+N10+N11)**2)/(N00+N01+N10+N11)**5)
        if val > 0:
            return np.sqrt(val)
        else:
            return 0.0
            
def get_fidelity(ZZ_corr,XX_corr,XmX_corr,ro_correct=True, psi1=True, F0a=0.905, F0b=0.805, 
                    F1a=0.9937, F1b=0.998, dF0a=0.01, dF0b=0.01, dF1a=0.01, dF1b=0.01):

    ZZ=   ro_c_F_odd(ZZ_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    dZZ=     ro_c_dF(ZZ_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    ZZS=    ro_c_F_S(ZZ_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    dZZS=  ro_c_dF_S(ZZ_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    XX=  ro_c_F_even(XX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)\
                        if psi1 else ro_c_F_odd(XX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    dXX=     ro_c_dF(XX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    XmX= ro_c_F_odd(XmX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)\
                        if psi1 else ro_c_F_even(XmX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    dXmX=   ro_c_dF(XmX_corr, F0a, F0b, F1a, F1b, dF0a, dF0b, dF1a, dF1b, ro_correct=ro_correct)
    
    XXavg = XmX#(XX/dXX**2 + XmX/dXmX**2)/(1/dXX**2+1/dXmX**2)
    dXXavg= dXmX#np.sqrt(1/(1/dXX**2+1/dXmX**2))
    
    F=ZZ/2. - ZZS + (XXavg - 1/2.)
    dF= np.sqrt(1/4.*dZZ**2 + dZZS**2 + dXXavg**2)
    
    return F,dF,{'ZZ': ZZ, 'dZZ': dZZ, 'ZZS': ZZS, 'dZZS': dZZS, 
                    'XX': XX, 'dXX': dXX, 'XmX': XmX, 'dXmX':dXmX, 
                    'XXavg':XXavg,'dXXavg':dXXavg}
    
class FidelityAnalysis:

    def __init__(self,**kw):
        self.basepath=kw.pop('basepath',r'D:\analysis\data\lde')
        self.folder_ZZ=os.path.join(self.basepath,'ZZ')
        self.folder_XX=os.path.join(self.basepath,'XX')
        self.folder_XmX=os.path.join(self.basepath,'X-X')
        self.F0b = 0.822
        self.F0a = 0.9210 
        self.F1b = 0.9892
        self.F1a = 0.9965
        
        self.dF0b = 0.0075
        self.dF0a = 0.0029 
        self.dF1b = 0.01
        self.dF1a = 0.0009

    
    def analyse_fidelity(self,state,ro_correct=True,w_start = (637,666), w_length=270, w_dt=-1, **kw):
        
        self.ZZ_an=lde_analysis.LDEAnalysis()
        self.XX_an=lde_analysis.LDEAnalysis()
        self.XmX_an=lde_analysis.LDEAnalysis()

        self.ZZ_an.analyse_lde_from_dir(self.folder_ZZ, w_start = w_start, 
        w_length=w_length, w_dt=w_dt, analyse_g2=False, **kw)
        self.ZZ_an.filter_on_gatephase()
        self.XX_an.analyse_lde_from_dir(self.folder_XX, w_start = w_start, 
        w_length=w_length, w_dt=w_dt, analyse_g2=False, **kw)
        self.XX_an.filter_on_gatephase()
        self.XmX_an.analyse_lde_from_dir(self.folder_XmX, w_start = w_start, 
        w_length=w_length, w_dt=w_dt, analyse_g2=False, **kw)
        self.XmX_an.filter_on_gatephase()      
        
        self.state=state
        self.w_start=w_start
        self.w_length=w_length
        self.w_dt=w_dt
        self.state=state
        if state=='psi1':
            psi1=True
            self.ZZ_corr =   self.ZZ_an.total_corr_00 +  self.ZZ_an.total_corr_11
            self.XX_corr =   self.XX_an.total_corr_00 +  self.XX_an.total_corr_11
            self.XmX_corr = self.XmX_an.total_corr_00 + self.XmX_an.total_corr_11
        elif state == 'psi2':
            psi1=False
            self.ZZ_corr=   self.ZZ_an.total_corr_01 +  self.ZZ_an.total_corr_10
            self.XX_corr=   self.XX_an.total_corr_01 +  self.XX_an.total_corr_10
            self.XmX_corr= self.XmX_an.total_corr_01 + self.XmX_an.total_corr_10       
        else:
            raise(Exception('Unknown state' + state))

        self.F, self.dF,self.F_dict = \
                get_fidelity(self.ZZ_corr,self.XX_corr,self.XmX_corr,
                            F0a=self.F0a, F0b=self.F0b, F1a=self.F1a, F1b=self.F1b, 
                            dF0a=self.dF0a, dF0b=self.dF0b, dF1a=self.dF1a, dF1b=self.dF1b,
                            ro_correct=ro_correct, psi1=psi1)
                
        self.dZZ=1/4.*self.F_dict['dZZ']**2 + self.F_dict['dZZS']**2 
        self.dXX=self.F_dict['dXXavg']**2
            
        self.F2stdev = (self.F-0.5)/self.dF

        self.print_fidelity()
       

    def reanalyse_fidelity(self,state,ro_correct=True, w_start = (637,666), w_length=150, w_dt=50, **kw):
    
        ZZ_total,ZZ_00,ZZ_01,ZZ_10,ZZ_11 = self.ZZ_an.reanalyse_lde(w_start = w_start, 
        w_length=w_length, w_dt=w_dt, apply_to_self=False, **kw)
        
        XX_total,XX_00,XX_01,XX_10,XX_11 = self.XX_an.reanalyse_lde(w_start = w_start, 
        w_length=w_length, w_dt=w_dt, apply_to_self=False, **kw)
        
        XmX_total,XmX_00,XmX_01,XmX_10,XmX_11 = self.XmX_an.reanalyse_lde(w_start = w_start, 
        w_length=w_length, w_dt=w_dt, apply_to_self=False, **kw)
        
        self.state=state
        self.w_start=w_start
        self.w_length=w_length
        self.w_dt=w_dt
        
        if state=='psi1':
            psi1=True
            self.ZZ_corr =   ZZ_00 +  ZZ_11
            self.XX_corr =   XX_00 +  XX_11
            self.XmX_corr = XmX_00 + XmX_11
        elif state == 'psi2':
            psi1=False
            self.ZZ_corr=  ZZ_01 +  ZZ_10
            self.XX_corr=  XX_01 +  XX_10
            self.XmX_corr= XmX_01 + XmX_10       
        else:
            raise(Exception('Unknown state' + state))

        self.F, self.dF,self.F_dict = \
                get_fidelity(self.ZZ_corr,self.XX_corr,self.XmX_corr,
                            F0a=self.F0a, F0b=self.F0b, F1a=self.F1a, F1b=self.F1b, 
                            dF0a=self.dF0a, dF0b=self.dF0b, dF1a=self.dF1a, dF1b=self.dF1b,
                            ro_correct=ro_correct, psi1=psi1)
                
        self.dZZ=1/4.*self.F_dict['dZZ']**2 + self.F_dict['dZZS']**2 
        self.dXX=self.F_dict['dXXavg']**2
            
        self.F2stdev = (self.F-0.5)/self.dF
        
        self.print_fidelity()

    def print_fidelity(self):
        
        print 'F:',self.F
        print 'sigma_F:', self.dF
        print 'total N:', self.ZZ_corr.sum()+self.XX_corr.sum()+self.XmX_corr.sum()
        print self.ZZ_corr,self.XX_corr,self.XmX_corr
        print 'sigma_ZZ_contrib:', self.dZZ
        print 'sigma_XX_contrib:', self.dXX

        print '(F-0.5)/sigma:', self.F2stdev
    
    def plot(self,save_path=r'D:\analysis\output'):
        
        ZZ_corr_err=sscorr.get_correlation_errors(self.ZZ_corr)
        XX_corr_err=sscorr.get_correlation_errors(self.XX_corr)
        XmX_corr_err=sscorr.get_correlation_errors(self.XmX_corr)
        
        fig=plt.figure()
        plt.subplot(231)
        sscorr.plot_uncorrected(self.ZZ_corr, ZZ_corr_err)
        plt.subplot(232)
        sscorr.plot_uncorrected(self.XX_corr, XX_corr_err)
        plt.subplot(233)
        sscorr.plot_uncorrected(self.XmX_corr, XmX_corr_err)
        
        ZZ_c_corr, ZZ_c_corr_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(self.ZZ_corr, 
                        self.F0a, self.F0b, self.F1a, self.F1b, 
                        dF0a=self.dF0a, dF0b=self.dF0b, dF1a=self.dF1a, dF1b=self.dF1b,
                        verbose=False, return_error_bars=True)
        XX_c_corr, XX_c_corr_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(self.XX_corr, 
                        self.F0a, self.F0b, self.F1a, self.F1b, 
                        dF0a=self.dF0a, dF0b=self.dF0b, dF1a=self.dF1a, dF1b=self.dF1b,
                        verbose=False, return_error_bars=True)
        XmX_c_corr, XmX_c_corr_err=sscorr.ssro_correct_twoqubit_state_photon_numbers(self.XmX_corr, 
                        self.F0a, self.F0b, self.F1a, self.F1b, 
                        dF0a=self.dF0a, dF0b=self.dF0b, dF1a=self.dF1a, dF1b=self.dF1b,
                        verbose=False, return_error_bars=True)
        plt.subplot(234)
        sscorr.plot_corrected(ZZ_c_corr, ZZ_c_corr_err)
        plt.subplot(235)
        sscorr.plot_corrected(XX_c_corr, XX_c_corr_err)
        plt.subplot(236)
        sscorr.plot_corrected(XmX_c_corr, XmX_c_corr_err)
        
        plt.suptitle('Correlations for state '+self.state+' with w_start ' + str(self.w_start) +\
                        ', w_length ' +str(self.w_length)+ ', w_dt ' + str(self.w_dt) + \
                        '\n and Fidelity F='+str(self.F)+'+/-'+str(self.dF))
        fig.set_size_inches(12,12)
        fig.savefig(os.path.join(save_path,'fidelity'+'_'+self.state+'.pdf'))
        return fig
    
    def save(self, save_path=r'D:\analysis\output'):        
        filename=os.path.join(save_path,'fidelity_analysis_'+self.state+'.pkl')
        #print filename
        f=open(filename,'wb')
        pickle.dump(self,f)
        f.close()
        
def load_previous_analysis(filename):
    f=open(filename,'rb')
    a = pickle.load(f)
    f.close()
    return a
>>>>>>> e5231864cbe6d2afad8f94e692ac4014d2bdb65f
