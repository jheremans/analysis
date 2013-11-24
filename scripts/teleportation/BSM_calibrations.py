import os, sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pprint

from analysis.lib.fitting import fit
from analysis.lib.tools import plot
from measurement.lib.tools import toolbox

from analysis.lib.m2.ssro import mbi
from analysis.lib.m2.ssro import ssro
reload(ssro)

# adapt
name = 'hans-sil1'
pi2_4mhz_value = 1. - 0.473

def stage_1_calibrations():
    # ssro first
    fig, ax = plt.subplots(1,1, figsize=(4,4))
  
    print 80*'='
    print 'Slow pi pulse for MBI'
    print 80*'='
    slow_pi_length = slow_pi(ax)

def stage_2a_calibrations():
    fig, ax = plt.subplots(1,1)

    print 80*'='
    print '4 MHz electron Rabi'
    print 80*'='
    CORPSE_freq = fast_rabi(ax)

def stage_2_calibrations():
    fig, ([ax1, ax2],[ax3, ax4]) = plt.subplots(2,2, figsize=(10,10))
    
    print 80*'='
    print 'fast pi pulse'
    print 80*'='
    fast_pi_amp = fast_pi(ax1)
    
    print 80*'='
    print 'CORPSE pi'
    print 80*'='
    CORPSE_pi_amp = CORPSE_pi(ax2)

    print 80*'='
    print 'pi2pi'
    print 80*'='
    pi2pi_pi_amp = pi_pi2pi(ax3)

    print 80*'='
    print 'pi2pi, m_I=0'
    print 80*'='
    pi2pi_pi_mI0_amp = pi_pi2pi_mI0(ax4)

def stage_3a_calibrations():
    fig, ax = plt.subplots(1,1, figsize = (5,3))
    print 80*'='
    print 'Nitrogen driving frequency'
    print 80*'='
    nitrogen_frequency = N_frq(ax)

def stage_3_calibrations():
    fig, ax = plt.subplots(1,1, figsize = (5,3))
    print 80*'='
    print 'Nitrogen Rabi oscillation'
    print 80*'='
    nitrogen_rabi = N_rabi(ax)  

def stage_4_calibrations():
    fig, ax = plt.subplots(1,1, figsize = (5,3))
    print 80*'='
    print 'CORPSE pi phase shift (look for the maximum)'
    print 80*'='
    CORPSE_pi_phase = CORPSE_phase(ax)    

def stage_5_calibrations():
    # fig, ax = plt.subplots(1,1, figsize = (5,3))
    # print 80*'='
    # print 'UNROT evolution time calibration'
    # print 80*'='
    # UNROT_evolution_time = UNROT_evtime(ax)

    fig, ax = plt.subplots(1,1, figsize = (5,3))
    print 80*'='
    print 'UNROT evolution time calibration (look for the maximum)'
    print 80*'='
    UNROT_evolution_time = UNROT_evtime_large_range(ax)

def stage_6_calibrations():
    # fig, ax = plt.subplots(1,1, figsize = (5,3))
    # print 80*'='
    # print 'UNROT evolution time calibration'
    # print 80*'='
    # UNROT_evolution_time = UNROT_evtime(ax)

    fig, ax = plt.subplots(1,1, figsize = (5,3))
    print 80*'='
    print 'spin echo (LDE-BSM) calibrations (look for the minimum)'
    print 80*'='
    spin_echo_time = spin_echo(ax)

def stage_7_calibrations():
    # fig, ax = plt.subplots(1,1, figsize = (5,3))
    # print 80*'='
    # print 'Hadamard phase calibration (aim for 10 high)'
    # print 80*'='
    # Hadamard_phase_cal = Hadamard_phase(ax, do_print_text=True)

    fig, ax = plt.subplots(1,1, figsize = (5,3))
    print 80*'='
    print 'Hadamard phase calibration (aim for 10 high)'
    print 80*'='
    Hadamard_phase_cal = Hadamard_phase_large_range(ax, do_print_text=True)

def stage_8_calibrations():
    # fig, ax = plt.subplots(1,1, figsize = (5,3))
    # print 80*'='
    # print 'Hadamard evolution time calibration (aim for 01 high)'
    # print 80*'='
    # Hadamard_evolution_time = Hadamard_ev_time(ax)

    fig, ax = plt.subplots(1,1, figsize = (5,3))
    print 80*'='
    print 'Hadamard evolution time calibration (aim for 01 high)'
    print 80*'='
    Hadamard_evolution_time = Hadamard_ev_time_large_range(ax)


### generic calibration functions, don't use those directly
def calibrate_epulse_amplitude(folder, ax, *args, **kw):
    double_ro = kw.pop('double_ro', 'False')

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    a.plot_results_vs_sweepparam(ax=ax, name='adwindata')
    
    x = a.sweep_pts.reshape(-1)[:]
    if double_ro == 'electron':
        y = a.p0[:,0]
    elif double_ro == 'nitrogen':
        y = a.p0[:,1]
    else:
        y = a.p0.reshape(-1)[:]
    
    x0 = fit.Parameter(args[0], 'x0')
    of = fit.Parameter(args[1], 'of')
    a = fit.Parameter(args[2], 'a')
    fitfunc_str = '(1-of) + a (x-x0)**2'
    
    def fitfunc_parabolic(x):
        return (1.-of()) + a() * (x-x0())**2
    
    fit_result = fit.fit1d(x,y, None, p0=[of, a, x0], fitfunc=fitfunc_parabolic,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def fit_correlation_parabolic(folder, ax, *args, **kw):
    which_correlation = kw.pop('which_correlation', 2)

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_correlations(name = 'adwindata')
    a.plot_results_vs_sweepparam(ax=ax, mode = 'correlations')
        
    x = a.sweep_pts.reshape(-1)[:]
    y = a.normalized_correlations[:,which_correlation]
    
    x0 = fit.Parameter(args[0], 'x0')
    of = fit.Parameter(args[1], 'of')
    a = fit.Parameter(args[2], 'a')
    fitfunc_str = '(1-of) + a (x-x0)**2'
    
    def fitfunc_parabolic(x):
        return (1.-of()) + a() * (x-x0())**2
    
    fit_result = fit.fit1d(x,y, None, p0=[of, a, x0], fitfunc=fitfunc_parabolic,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def fit_correlation_oscillation(folder, ax, *args, **kw):
    which_correlation = kw.pop('which_correlation', 2)

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_correlations(name = 'adwindata')
    a.plot_results_vs_sweepparam(ax=ax, mode = 'correlations')
        
    x = a.sweep_pts.reshape(-1)[:]
    y = a.normalized_correlations[:,which_correlation]
    
    x0 = fit.Parameter(args[0], 'x0')
    f = fit.Parameter(args[1], 'f')
    A = fit.Parameter(args[2], 'A')  
    o = fit.Parameter(0.25, 'o')
    fitfunc_str = '(1 - A) + A * cos(2pi f (x - x0))'

    def fitfunc(x) : 
        return o() + A() * np.cos(2*np.pi*f()*(x - x0()))

    
    fit_result = fit.fit1d(x,y, None, p0=[f, A, x0,o], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def fit_population_vs_detuning(folder, ax, *args):
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC() # a.get_N_ROC(0.97,0.03,0.96,0.01,0.93,0.01)#(0.99, 0.02, 0.94, 0.005)
    ax = a.plot_results_vs_sweepparam(ret='ax', name='adwindata')

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

    guess_a = 1
    guess_A = -0.8
    guess_F = 0.005
    guess_x0 = args[0]

    fit_result=fit.fit1d(x, y,
            rabi.fit_population_vs_detuning,
            guess_a, guess_A, guess_F, guess_x0,
            fixed=[], do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(a.sweep_pts[0],a.sweep_pts[-1],201), 
        ax=ax, plot_data=False)

    return fit_result

def fit_linear(folder, ax, value, *args):
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    a.plot_results_vs_sweepparam(ax=ax, name='adwindata')
    
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

    a = fit.Parameter(args[0], 'a')
    b = fit.Parameter(args[0], 'b')
    fitfunc_str = 'a x + b'

    def fitfunc(x):
        return a() * x + b()

    fit_result = fit.fit1d(x,y, None, p0=[a, b], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def calibrate_epulse_rabi(folder, ax, *args, **kws):
    fit_x0 = kws.pop('fit_x0', True)
    fit_k = kws.pop('fit_k', True)
    double_ro = kws.pop('double_ro', 'False')
    guess_x0 = kws.pop('guess_x0', 0)

    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    a.plot_results_vs_sweepparam(ax=ax, name='adwindata')
    
    x = a.sweep_pts.reshape(-1)[:]
    if double_ro == 'electron':
        y = a.p0[:,0]
    elif double_ro == 'nitrogen':
        y = a.p0[:,1]
    else:
        y = a.p0.reshape(-1)[:]

    f = fit.Parameter(args[0], 'f')
    A = fit.Parameter(args[1], 'A')
    x0 = fit.Parameter(guess_x0, 'x0')
    k = fit.Parameter(0, 'k')
    p0 = [f, A]
    if fit_x0:
        p0.append(x0)
    if fit_k:
        p0.append(k)

    fitfunc_str = '(1 - A) + A * exp(-kx) * cos(2pi f (x - x0))'

    def fitfunc(x) : 
        return (1.-A()) + A() * np.exp(-k()*x) * \
            np.cos(2*np.pi*f()*(x - x0()))

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def calibrate_Npulse_rabi(folder, ax, *args):
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    a.plot_results_vs_sweepparam(ax=ax, name='adwindata')
    
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0[:,0]

    f = fit.Parameter(args[0], 'f')
    A = fit.Parameter(args[1], 'A')
    o = fit.Parameter(0.5, 'o')
    x0 = fit.Parameter(0, 'x0')
    p0 = [f, A, x0, o]

    fitfunc_str = 'o + A * cos(2pi f x - 2pi phi)'

    def fitfunc(x):
        return o() + A() * np.cos(2*np.pi*(f()*(x - x0())))

    fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

def epulse_fidelity(folder, ax, *args):
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    a.plot_results_vs_sweepparam(ax=ax, name='adwindata')
    
    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]
    
    of = fit.Parameter(args[0], 'of')
    fitfunc_str = '(1-of)'
    
    def fitfunc_fid(x):
        return (1.-of())
    
    fit_result = fit.fit1d(x,y, None, p0=[of], fitfunc=fitfunc_fid,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
    plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
        
    return fit_result

### end generic calibration functions

### actual functions to call, for specific calibrations    
def slow_pi(ax=None):
    folder = toolbox.latest_data('cal_slow_pi_'+name)
    if ax==None:
        fig,ax = plt.subplots(1,1)
    fit_result = calibrate_epulse_rabi(folder, ax, 0.015, 1.0, 0., fit_x0 = False, fit_k = False)

    f = fit_result['params_dict']['f']
    u_f = fit_result['error_dict']['f']
    #ax.text(500, 0.4, 'pi = (%.0f +/- %.0f) ns' % (0.5/f, 0.5/f**2 * u_f),
    #    va='bottom', ha='left')
    ax.text(0.015, 0.4, 'pi = (%.4f +/- %.4f) V' % (0.5/f, 0.5/f**2 * u_f),
        va='bottom', ha='left')
    
    return (0.5/f, 0.5/f**2 * u_f)

def fast_rabi(ax=None):
    folder = toolbox.latest_data('cal_fast_rabi'+name)
    if ax==None:
        fig,ax = plt.subplots(1,1)
    fit_result = calibrate_epulse_rabi(folder, ax, 1./125, 0.5, fit_k=False)

    f = fit_result['params_dict']['f']
    u_f = fit_result['error_dict']['f']
    ax.text(100, 0.9, '$f_r$ = (%.3f +/- %.3f) MHz' % (f*1e3, u_f*1e3),
        va='bottom', ha='left')

    return (f*1e3, u_f*1e3)

def fast_pi(ax=None, do_print_text=True):
    folder = toolbox.latest_data('cal_fast_pi_'+name)

    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.8,  1, 0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    
    if do_print_text:
        ax.text(0.65, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A 


def fast_pi2(ax=None):
    folder = toolbox.latest_data('cal_fast_pi_over_2') 
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = fit_linear(folder, ax, -1,  1)
    a = fit_result['params_dict']['a']
    b = fit_result['params_dict']['b']
    u_a = fit_result['error_dict']['a']
    u_b = fit_result['error_dict']['b']
    A = (0.5 - b) / a
    u_A = np.sqrt(A**2) * np.sqrt ( (u_a/a)**2 + (u_b/b)**2 )
    ax.text(0.65, 0.8, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A 


def CORPSE_pi(ax=None, do_print_text=True):
    folder = toolbox.latest_data('CORPSEPiCalibration') # _'+name)
    
    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.53,  1, 0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    
    if do_print_text:
        ax.text(0.53, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A 

def N_frq(ax=None):
    folder = toolbox.latest_data('NMR_frq_scan')
    if ax==None:
        fig,ax = plt.subplots(1,1)
    fit_result = fit_population_vs_detuning(folder, ax, 7.135)

    x0 = fit_result['params_dict']['x0']
    u_x0 = fit_result['error_dict']['x0']
    ax.text(7.1, 0.3, '$N_{frq}$ = (%.3f +/- %.3f) MHz' % (x0, u_x0),
        va='bottom', ha='left')

    return (f*1e3, u_f*1e3)


def N_rabi(ax=None):
    folder = toolbox.latest_data('Nuclear_rabi')
    if ax==None:
        fig,ax = plt.subplots(1,1)
    fit_result = calibrate_Npulse_rabi(folder, ax, 1./100, 0.5)

    f = fit_result['params_dict']['f']
    u_f = fit_result['error_dict']['f']
    ax.text(100, 0.9, '$f_r$ = (%.3f +/- %.3f) MHz' % (f*1e3, u_f*1e3),
        va='bottom', ha='left')
    ax.text(100, 0.7, 'pi pulse = (%.3f +/- %.3f) us' % (1./(2*f), 1./(2*f)*u_f/f),
        va='bottom', ha='left')

    return (f*1e3, u_f*1e3)


def CORPSE_phase(ax=None, do_print_text=True):
    folder = toolbox.latest_data('CalibrateCORPSEPiPhase') # _'+name)
    
    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 107, 1, 0., double_ro='electron')
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    
    if do_print_text:
        ax.text(100, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A 


def UNROT_evtime(ax=None, do_print_text=True):
    folder = toolbox.latest_data('Calibrate_UNROT_X_timing') 
    
    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 51, 0, 1., double_ro='nitrogen')
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    
    if do_print_text:
        ax.text(51.03, 0.5, 'A = (%.3f +/- %.3f) us' % (A, u_A))

    return A, u_A 

def UNROT_evtime_large_range(ax=None, do_print_text=True):
    folder = toolbox.latest_data('Calibrate_UNROT_X_timing') 
    
    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
    guess_x0 =  52.5   
    fit_result = calibrate_epulse_rabi(folder, ax, 1./0.45, 0.5, guess_x0=guess_x0, double_ro='nitrogen', fit_x0 = True)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    
    if do_print_text:
        ax.text(guess_x0, 0.5, 'A = (%.3f +/- %.3f) us' % (A, u_A))

    return A, u_A 

def spin_echo(ax=None, do_print_text=True):
    folder = toolbox.latest_data('calibrate_echo') # _'+name)
    
    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 200, 1, 0., double_ro='electron')
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    
    if do_print_text:
        ax.text(100, 0.5, 'echo = (%.3f +/- %.3f) ns' % (A, u_A))

    return A, u_A 



def Hadamard_phase(ax=None, do_print_text=True):
    folder = toolbox.latest_data('TestBSM_LDE_calibrate_H_phase') # _'+name)
    
    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
        
    fit_result = fit_correlation_parabolic(folder, ax, 100, 1, 0., which_correlation=0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    
    if do_print_text:
        ax.text(10, 0.8, 'H phase = (%.3f +/- %.3f) ' % (A, u_A))
    return A, u_A 

def Hadamard_phase_large_range(ax=None, do_print_text=True):
    folder = toolbox.latest_data('TestBSM_LDE_calibrate_H_phase') # _'+name)
    
    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
        
    fit_result = fit_correlation_oscillation(folder, ax, 100, 1./360, 0.25, which_correlation=0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    
    if do_print_text:
        ax.text(10, 0.8, 'H phase = (%.3f +/- %.3f) ' % (A, u_A))
    return A, u_A 

def Hadamard_ev_time(ax=None, do_print_text=True):
    folder = toolbox.latest_data('calibrate_H_ev_time') # _'+name)

    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
        
    fit_result = fit_correlation_parabolic(folder, ax, 35.0, 1, 0., which_correlation=3)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']

    if do_print_text:
        ax.text(35.0, 0.8, 'evo time = (%.3f +/- %.3f) us' % (A, u_A))

    return A, u_A 

def Hadamard_ev_time_large_range(ax=None, do_print_text=True):
    folder = toolbox.latest_data('calibrate_H_ev_time') # _'+name)
    
    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
        
    fit_result = fit_correlation_oscillation(folder, ax, 35.08, 1./0.05/10., 0.25, which_correlation=3)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    
    if do_print_text:
        ax.text(35.07, 0.8, 'H phase = (%.4f +/- %.4f) ' % (A, u_A))
    return A, u_A 

def pi_pi2pi(ax=None, do_print_text=True):
    folder = toolbox.latest_data('cal_pi2pi_pi_'+name)
    
    if ax==None:
        do_print_text = False
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.11, 1, 0)#ax, 0.085, 1, 0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    
    if do_print_text:
        ax.text(0.1, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A

def pi_pi2pi_mI0(ax=None):
    folder = toolbox.latest_data('cal_pi2pi_pi_mI0')
    
    if ax==None:
        fig,ax = plt.subplots(1,1)
        
    fit_result = calibrate_epulse_amplitude(folder, ax, 0.11, 1, 0)
    A = fit_result['params_dict']['x0']
    u_A = fit_result['error_dict']['x0']
    ax.text(0.1, 0.5, 'A = (%.3f +/- %.3f) V' % (A, u_A))

    return A, u_A

# #### Not fully finished yet

# def pi_397ns(ax=None):
#     folder = toolbox.latest_data('cal_397ns_pi_'+name)
    
#     if ax==None:
#         fig,ax = plt.subplots(1,1)
        
#     fit_result = calibrate_epulse_amplitude(folder, ax, 0.16,  1, 0)

# def hard_pi_duration(ax=None):
#     folder = toolbox.latest_data('cal_hard_pi_dur_'+name)
    
#     if ax==None:
#         fig,ax = plt.subplots(1,1)
        
#     fit_result = calibrate_epulse_amplitude(folder, ax, 120,  1, 0)

# def hard_pi_amplitude(ax=None):
#     folder = toolbox.latest_data('cal_hard_pi_amp_'+name)
    
#     if ax==None:
#         fig,ax = plt.subplots(1,1)
        
#     fit_result = calibrate_epulse_amplitude(folder, ax, 0.85,  1, 0)
          
    
# def CORPSE(ax=None):
#     folder = toolbox.latest_data('Test_CORPSE')
    
#     if ax==None:
#         fig,ax = plt.subplots(1,1)
        
#     fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)
    
# def CORPSE_60(ax=None):
#     folder = toolbox.latest_data('cal_CORPSE60')
    
#     if ax==None:
#         fig,ax = plt.subplots(1,1)
        
#     fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)

# def CORPSE_60_duration(ax=None):
#     folder = toolbox.latest_data('CORPSE_vs_60length')
    
#     if ax==None:
#         fig,ax = plt.subplots(1,1)
        
#     fit_result = calibrate_epulse_amplitude(folder, ax, 52,  1, 0)
     
    
# def CORPSE_m300(ax=None):
#     folder = toolbox.latest_data('cal_CORPSEm300')
    
#     if ax==None:
#         fig,ax = plt.subplots(1,1)
        
#     fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)
    
# def CORPSE_420(ax=None):
#     folder = toolbox.latest_data('cal_CORPSE420')
    
#     if ax==None:
#         fig,ax = plt.subplots(1,1)
        
#     fit_result = calibrate_epulse_amplitude(folder, ax, 0.68,  1, 0)
    
# def CORPSE_fidelity(ax=None):
#     folder = toolbox.latest_data('CORPSE')
    
#     if ax==None:
#         fig,ax = plt.subplots(1,1)
        
#     fit_result = epulse_fidelity(folder, ax, 1)
  
                                
# if __name__ == '__main__':
#     # calibrate_all()
#     stage_1_calibrations()