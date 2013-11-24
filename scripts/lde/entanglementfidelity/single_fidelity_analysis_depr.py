import os
from analysis.lde import fidelities, lde_analysis

        
folder_ZZ=r'D:\bjhensen\data\ZZ'
folder_XX=r'D:\bjhensen\data\XX'
folder_XmX=r'D:\bjhensen\data\X-X'

ro_correct=True
state='psi1' #or 'psi2'
analysis_kw= {'w_length':100,
              'w_dt': 50,
              'w_start': (634,662),
              'analyse_g2': False,
              }


ZZ_an=lde_analysis.LDEAnalysis()
XX_an=lde_analysis.LDEAnalysis()
XmX_an=lde_analysis.LDEAnalysis()

ZZ_an.analyse_lde_from_dir(folder_ZZ, **analysis_kw)
ZZ_an.filter_on_gatephase()
XX_an.analyse_lde_from_dir(folder_XX, **analysis_kw)
XX_an.filter_on_gatephase()
XmX_an.analyse_lde_from_dir(folder_XmX, **analysis_kw)
XmX_an.filter_on_gatephase()      

if state=='psi1':
    psi1=True
    ZZ_corr =   ZZ_an.total_corr_00 +  ZZ_an.total_corr_11
    XX_corr =   XX_an.total_corr_00 +  XX_an.total_corr_11
    XmX_corr = XmX_an.total_corr_00 + XmX_an.total_corr_11
elif state == 'psi2':
    psi1=False
    ZZ_corr=   ZZ_an.total_corr_01 +  ZZ_an.total_corr_10
    XX_corr=   XX_an.total_corr_01 +  XX_an.total_corr_10
    XmX_corr= XmX_an.total_corr_01 + XmX_an.total_corr_10       
else:
    raise(Exception('Unknown state' + state))

F, dF,_dict = \
        fidelities.get_fidelity(ZZ_corr,XX_corr,XmX_corr, ro_correct=ro_correct, psi1=psi1)
        
dZZ=1/4.*_dict['dZZ']**2 + _dict['dZZS']**2 
dXX=_dict['dXXavg']**2
    
F2stdev = (F-0.5)/dF


print 'F:',F
print 'sigma_F:', dF

print 'sigma_ZZ_contrib:', dZZ
print 'sigma_XX_contrib:', dXX

print '(F-0.5)/sigma:', F2stdev

