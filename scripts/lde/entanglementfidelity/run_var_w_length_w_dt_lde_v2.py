<<<<<<< HEAD
import os

base = r'D:\measuring\analysis\data\lde'

folders=[os.path.join(base,'ZZ'),os.path.join(base,'XX'),os.path.join(base,'X-X')]
for folder in folders:
    w_lengths=range(10,100,20)
    w_dts=range(10,100,10)

    corr=zeros((len(w_lengths),len(w_dts),4))
    corr_00=zeros((len(w_lengths),len(w_dts),4))
    corr_01=zeros((len(w_lengths),len(w_dts),4))
    corr_10=zeros((len(w_lengths),len(w_dts),4))
    corr_11=zeros((len(w_lengths),len(w_dts),4))

    a=lde_analysis.LDEAnalysis()
    a.analyse_lde_from_dir(folder, w_start=(637,666), w_length=max(w_lengths), w_dt=max(w_dts), analyse_g2=False)

    for i,w_length in enumerate(w_lengths):
        print i+1,'/',len(w_lengths)
        for j,w_dt in enumerate(w_dts):
            corr[i,j], corr_00[i,j], corr_01[i,j],corr_10[i,j],corr_11[i,j]=\
                    a.reanalyse_lde(w_start=(637,666), w_length=w_length, w_dt=w_dt)
            
    savez(os.path.join(folder,'all_ro_cor_fast.npz'), 
        corr=corr, corr_00=corr_00, corr_01=corr_01, corr_10=corr_10, corr_11=corr_11,
        w_lengths=w_lengths, w_dts=w_dts)
=======
import os
folders=[r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde/ZZ',r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde/XX',r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/lde/X-X']
for folder in folders:
    w_lengths=range(10,100,20)
    w_dts=range(10,100,10)

    corr=zeros((len(w_lengths),len(w_dts),4))
    corr_00=zeros((len(w_lengths),len(w_dts),4))
    corr_01=zeros((len(w_lengths),len(w_dts),4))
    corr_10=zeros((len(w_lengths),len(w_dts),4))
    corr_11=zeros((len(w_lengths),len(w_dts),4))

    a=lde_analysis.LDEAnalysis()
    a.analyse_lde_from_dir(folder, w_start=(637,666), w_length=max(w_lengths), w_dt=max(w_dts), analyse_g2=False)

    for i,w_length in enumerate(w_lengths):
        print i+1,'/',len(w_lengths)
        for j,w_dt in enumerate(w_dts):
            corr[i,j], corr_00[i,j], corr_01[i,j],corr_10[i,j],corr_11[i,j]=\
                    a.reanalyse_lde(w_start=(637,666), w_length=w_length, w_dt=w_dt)
            
    savez(os.path.join(folder,'all_ro_cor_fast.npz'), 
        corr=corr, corr_00=corr_00, corr_01=corr_01, corr_10=corr_10, corr_11=corr_11,
        w_lengths=w_lengths, w_dts=w_dts)
>>>>>>> 998bfb9d754ee59f58c639b71f91f58a0a5b6921
