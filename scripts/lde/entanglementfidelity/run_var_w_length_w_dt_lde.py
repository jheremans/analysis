import os
folder=r'D:\bjhensen\data\ZZ'
w_lengths=range(60,290,20)
w_dts=range(10,100,10)

corr=zeros((len(w_lengths),len(w_dts),4))
corr_00=zeros((len(w_lengths),len(w_dts),4))
corr_01=zeros((len(w_lengths),len(w_dts),4))
corr_10=zeros((len(w_lengths),len(w_dts),4))
corr_11=zeros((len(w_lengths),len(w_dts),4))



for i,w_length in enumerate(w_lengths):
    for j,w_dt in enumerate(w_dts):
        a=lde_analysis.LDEAnalysis()
        a.analyse_lde_from_dir(folder, w_start=(637,666), w_length=w_length, w_dt=w_dt, analyse_g2=False)
        a.total_corr,a.total_corr_00,a.total_corr_01,a.total_corr_10,a.total_corr_11=a.filter_on_gatephase()
        corr[i,j]=a.total_corr
        corr_00[i,j]=a.total_corr_00
        corr_01[i,j]=a.total_corr_01
        corr_10[i,j]=a.total_corr_10
        corr_11[i,j]=a.total_corr_11
        
savez(os.path.join(folder,'all_ro_cor.npz'), 
    corr=corr, corr_00=corr_00, corr_01=corr_01, corr_10=corr_10, corr_11=corr_11,
    w_lengths=w_lengths, w_dts=w_dts)