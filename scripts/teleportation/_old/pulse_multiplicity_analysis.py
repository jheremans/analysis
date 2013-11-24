from analysis.lib import fitting
from analysis.lib.m2.ssro import mbi
from analysis.lib.fitting import fit
from measurement.lib.tools import toolbox
from analysis.lib.m2 import m2
from analysis.lib.m2.ssro import ssro
from analysis.lib.tools import plot

timestamp = None # '20130604232420'

if timestamp != None:
    folder = toolbox.data_from_time(timestamp)
else:
    folder = toolbox.latest_data('sweep_pipulse_cnt')

k = fit.Parameter(0.98, 'k')
pop = fit.Parameter(1., 'pop')
p0 = [k, pop]
fitfunc_str = '0.5 + (p0 - 0.5) exp(-2(1-k)x)'
def fitfunc(x):
    return 0.5 + (pop() - 0.5) * np.exp(-2*(1-k())*x)

a = mbi.MBIAnalysis(folder)
a.get_sweep_pts()
a.get_readout_results(name='adwindata')
a.get_electron_ROC()
#a.get_N_ROC(1.00, 0.02, 0.94, 0.01, 0.96, 0.01)
ax = a.plot_results_vs_sweepparam(ret='ax', name='adwindata')

x = a.sweep_pts.reshape(-1)[:]
y = a.p0.reshape(-1)[:]

# first separate even and odd number of pulses
evenidx = where(x%2==0)[0]
oddidx = where(x%2!=0)[0]

even_result = fit.fit1d(x[evenidx], y[evenidx], None, p0=[k], fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
plot.plot_fit1d(even_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)


odd_result = fit.fit1d(x[oddidx], y[oddidx], None, p0=p0, fitfunc=fitfunc,
        fitfunc_str=fitfunc_str, do_print=True, ret=True)
plot.plot_fit1d(odd_result, np.linspace(x[0],x[-1],201), ax=ax,
        plot_data=False, print_info=False)
