
import numpy as np
import pylab as plt
from  matplotlib import rc
from analysis.lib.spin import spin_control as sc

f_normal_mszero = r'D:\measuring\data\20130717\152845_Seg_RO__LT2_SIL10_normalRO_normalRO_5nW_initmszero/'
f_normal_msone = r'D:\measuring\data\20130717\152304_Seg_RO__LT2_SIL10_normalRO_normalRO_5nW_initmsone/'
f_seg_mszero = r'D:\measuring\data\20130717\150422_Seg_RO__LT2_SIL10_segRO_normalRO_5nW/'
f_seg_msone = r'D:\measuring\data\20130717\151626_Seg_RO__LT2_SIL10_segRO_normalRO_5nW_initmsone/'
ssro_folder=r'D:\measuring\data\20130717\115946_Adwin_SSRO__LT2_SIL10_Ey_A2_5nW/'

def get_double_RO_data(folder='',plot=True):
    L=9

    Ndet = np.zeros(L)
    fid = np.zeros(L)
    ufid = np.zeros(L)
    for n in (np.arange(L)):
        if (n<10):
            nn = '0'+str(n)
        else:
            nn=str(n)

        a = np.load (folder+'Seg_RO-0'+nn+'_segment_number.npz')
        Ndet [n] = np.sum(a['segment_number'][0:len(a['segment_number'])-1])
        
        b = np.load(folder+'Seg_RO-0'+nn+'_Spin_RO.npz')
        reps=len(b['counts'])
        final_ROclicks_norm=sum(b['counts'])/float(reps)
        f,uf = sc.get_electron_ROC(final_ROclicks_norm,float(reps),ssro_folder)
        fid [n] = f
        ufid=uf
        x=b['sweep_axis']
        xname=b['sweep_par_name']
        b.close()
        a.close()
    if plot:    

        plt.figure()
        plt.plot (x, Ndet/float(reps), 'ob')
        plt.ylim([0,1])
        plt.ylabel ('read-out fidelity')
        plt.xlabel (xname)
        plt.show()

        plt.figure()
        plt.plot (x, fid, 'ob')
        plt.ylim([0,1])
        plt.ylabel ('state fidelity after RO')
        plt.xlabel (xname)
        plt.show()
    return fid,ufid,x    

f_n_0,uf_n_0,x = get_double_RO_data(f_normal_mszero,plot=False)
f_n_1,uf_n_1,x = get_double_RO_data(f_normal_msone,plot=False)
f_s_0,uf_s_0,x = get_double_RO_data(f_seg_mszero,plot=False)
f_s_1,uf_s_1,x = get_double_RO_data(f_seg_msone,plot=False)
print x
print f_n_0
colors={
'fs':'RoyalBlue',
'fn':'Crimson'
}

plt.figure()
plt.errorbar (x, f_n_1,uf_n_1, label='input |1>',marker='o', linestyle='None',mfc=colors['fn'], mec=colors['fn'],ecolor=colors['fn'],color=colors['fn'],elinewidth=15)
plt.errorbar(x, f_n_0,uf_n_0, label='input |0>',marker='o', linestyle='None',mfc=colors['fs'], mec=colors['fs'],ecolor=colors['fs'],color=colors['fs'],elinewidth=15)
plt.ylim([-0.05,1.15])
plt.ylabel ('P(ms=0)')
plt.xlabel ('RO duration [us]')
plt.title('Normal Readout')
plt.legend(loc=1)
plt.show()

plt.figure()
plt.errorbar (x, f_s_1,uf_s_1, label='input |1>',marker='o', linestyle='None',mfc=colors['fn'], mec=colors['fn'],ecolor=colors['fn'],color=colors['fn'],elinewidth=15)
plt.errorbar(x, f_s_0,uf_s_0, label='input |0>',marker='o', linestyle='None',mfc=colors['fs'], mec=colors['fs'],ecolor=colors['fs'],color=colors['fs'],elinewidth=15)
plt.ylim([-0.05,1.15])
plt.ylabel ('P(ms=0)')
plt.xlabel ('RO duration [us]')
plt.title('Segmented Readout')
plt.legend(loc=1)
plt.show()

plt.figure()
plt.bar(1.125,f_s_1[L-1],width=0.75,color=colors['fn'])
plt.bar(0.125,f_s_0[L-1],width=0.75,color=colors['fs'])
plt.xticks([0.5,1.5],['prepare |0>','prepare |1>'])
plt.ylim([-0.05,1])
plt.ylabel('P (ms=0)')
plt.title('Segmented Readout')

plt.figure()
plt.bar(1.125,f_n_1[L-1],width=0.75,color=colors['fn'])
plt.bar(0.125,f_n_0[L-1],width=0.75,color=colors['fs'])
plt.xticks([0.5,1.5],['prepare |0>','prepare |1>'])
plt.ylim([-0.05,1])
plt.ylabel('P (ms=0)')
plt.title('Normal Readout')