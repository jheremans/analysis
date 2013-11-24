import os, sys
import numpy as np
import h5py
import logging

from matplotlib import pyplot as plt
from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi
from analysis.lib.fitting import fit
from analysis.lib.tools import toolbox
from analysis.lib.m2 import m2
from analysis.lib.m2.ssro import ssro
from analysis.lib.tools import plot

timestamp = None

nr_of_revivals = 11
revivals_nrs = np.arange(nr_of_revivals)


fidelities = np.ones(nr_of_revivals)*0
revivals = np.ones(nr_of_revivals)*0
u_fid = np.ones(nr_of_revivals)*0
u_rev = np.ones(nr_of_revivals)*0

i = 0
r_max = 7

for r in revivals_nrs:
  
    name = 'revival_{}_'.format(r)


    if timestamp != None:
        folder = toolbox.data_from_time(timestamp)
    else:
        folder = toolbox.latest_data(name)


    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC()
    ax = a.plot_results_vs_sweepparam(ret='ax', name='adwindata')

    #a = sequence.SequenceAnalysis(folder)
    #a.get_sweep_pts()
    #a.get_readout_results(name='ssro')
    #a.get_electron_ROC()
    #ax = a.plot_result_vs_sweepparam(ret='ax', name='ssro')

    x = a.sweep_pts.reshape(-1)[:]
    y = a.p0.reshape(-1)[:]

    x0 = fit.Parameter(r*100.4, 'x0')
    a = fit.Parameter(0.5, 'a')
    o = fit.Parameter(0.5, 'o')
    c = fit.Parameter(10, 'c')
    fitfunc_str = ''

    def fitfunc_gauss(x):
        return o() + a() * np.exp( -((x-x0())/ c())**2)

    if r > r_max:
        a = fit.Parameter(0, 'a')
        def fitfunc_linear(x):
            return o() + a() * np.exp( -((x-x0())/ c())**2)
        fit_result = fit.fit1d(x,y, None, p0=[o,a], fixed= [1], fitfunc=fitfunc_linear,
            fitfunc_str=fitfunc_str, do_print = False, ret=True)
    else:
        fit_result = fit.fit1d(x,y, None, p0=[o,x0,a,c], fitfunc=fitfunc_gauss,
            fitfunc_str=fitfunc_str, do_print=False, ret=True)

    #plot.plot_fit1d(fit_result, np.linspace(x[0],x[-1],201), ax=ax,
    #    plot_data=False)


    if r>r_max:
        fidelities[i] = fit_result['params_dict']['o']
        revivals[i] = r * 108.4
        u_fid[i] = fit_result['error_dict']['o']
    else:   
        fidelities[i] = fit_result['params_dict']['a'] + fit_result['params_dict']['o'] 
        revivals[i] = fit_result['params_dict']['x0']
        u_fid[i] = np.sqrt(fit_result['error_dict']['a']**2 + fit_result['error_dict']['o']**2)
        u_rev[i] = fit_result['error_dict']['x0']
    
    i = i + 1

exp_x = revivals
exp_y = fidelities
exp_u_y = u_fid


fig, ax = plt.subplots(1,1)

ax.set_title(''+'\n'+'DynamicalDecoupling_t2')

ax.errorbar(exp_x, exp_y, fmt='o',
                yerr=exp_u_y)
    
ax.set_xlabel('total free evolution time')
ax.set_ylabel(r'$F(|0\rangle)$')

ax.set_ylim(-0.05, 1.05)

fig.savefig(os.path.join(folder, 't2_result_vs_total_fet.png'),format='png')
        
a_exp = fit.Parameter(0.5, 'a_exp')
o_exp = fit.Parameter(0.5, 'o_exp')
c_exp = fit.Parameter(200, 'c_exp')
n = fit.Parameter(3, 'n')
#a2_exp = fit.Parameter(0.18, 'a2_exp')
#c2_exp = fit.Parameter(700, 'c2_exp')
#n2 = fit.Parameter(3, 'n2')

def fitfunc_exp(x):
    return o_exp() + a_exp() * np.exp( -(x/ c_exp())**n()) #* np.exp( -(x/ c2_exp())**n2())

fitfunc_str_exp = 'o + a * e^(-(x/c)^n)'

fit_result_exp = fit.fit1d(exp_x,exp_y, None, p0=[o_exp,a_exp,c_exp,n], 
        fixed = [0,1,3], fitfunc=fitfunc_exp,
        fitfunc_str=fitfunc_str_exp, do_print=True, ret=True)

plot.plot_fit1d(fit_result_exp, np.linspace(exp_x[0],exp_x[-1],201), ax=ax,
        plot_data=True)

t2 = fit_result_exp['params_dict']['c_exp']
u_t2 = fit_result_exp['error_dict']['c_exp']

ax.text(200, 0.5, 't2 = (%.3f +/- %.3f) us' % (t2, u_t2))