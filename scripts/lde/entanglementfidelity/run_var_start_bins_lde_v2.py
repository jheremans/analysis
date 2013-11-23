import os
folders=[r'D:\bjhensen\data\ZZ',r'D:\bjhensen\data\XX',r'D:\bjhensen\data\X-X']
for folder in folders:
    w_ch0s=range(631,642)
    w_ch1s=range(661,672)

    corr=zeros((len(w_ch0s),len(w_ch1s),4))
    corr_00=zeros((len(w_ch0s),len(w_ch1s),4))
    corr_01=zeros((len(w_ch0s),len(w_ch1s),4))
    corr_10=zeros((len(w_ch0s),len(w_ch1s),4))
    corr_11=zeros((len(w_ch0s),len(w_ch1s),4))

    a=lde_analysis.LDEAnalysis()
    a.analyse_lde_from_dir(folder, w_start=(630,660), w_length=70, w_dt=30, analyse_g2=False)

    for i,w_ch0 in enumerate(w_ch0s):
        print i+1,'/',len(w_ch0s)
        for j,w_ch1 in enumerate(w_ch1s):
            corr[i,j], corr_00[i,j], corr_01[i,j],corr_10[i,j],corr_11[i,j]=\
                    a.reanalyse_lde(w_start=(w_ch0,w_ch1), w_length=50, w_dt=10)
            
    savez(os.path.join(folder,'all_ro_cor_w_start.npz'), 
        corr=corr, corr_00=corr_00, corr_01=corr_01, corr_10=corr_10, corr_11=corr_11,
        w_ch0s=w_ch0s, w_ch1s=w_ch1s)