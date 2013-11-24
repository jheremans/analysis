
from numpy import *
import pylab as plt
from analysis.lib.fitting import fit, common
from analysis.lib.tools import plot


def loadData (fname = 'Z:/Diamond/Autobackup/LT2/data/20121211/160345_MBI_RO_N_vs_e_RO_steps_RO_no_shel_init+1_RO+1/SSRO_readout.dat'):

    data = {}
    i = 0
    with open(fname) as inf:
        N = []
        RO = []
        ROcorr = []
        errRO = []
        for line in inf:
            li = line.strip()
            if not li.startswith("#"):
                parts = line.split()
                N.append(float(parts[0]))
                RO.append(float(parts [1]))
                ROcorr.append(float(parts[2]))
                errRO.append(float(parts[3]))
                i = i+1


        data = {'Nreadouts': N, 'rawRO': RO, 'corrRO': ROcorr, 'errRO': errRO}
    return data


def plotNuclearState (N, matrix3):
    #input:
    #N; different readout steps
    #matrix3 is a matrix with 3 rows, corresponding respectively to values for mI = -, mI = 0, mI = +1
    
    Am1 = squeeze(asarray(matrix3[0, :]))
    A0 = squeeze(asarray(matrix3[1, :]))
    Ap1 = squeeze(asarray(matrix3[2, :]))
    
    fig = figure()
    plt.plot (N, Am1, 'b>', label='mI = -1')
    plt.plot (N, A0, 'ro', label='mI = 0')
    plt.plot (N, Ap1, 'kD', label = 'mI = +1')
    plt.rcParams.update({'font.size': 18})
    plt.xlabel ('number of read-out steps', fontsize = 18)
    plt.ylabel ('readout corrected nuclear spin fidelity', fontsize = 18)   
    plt.legend (loc =1, prop={'size':12})
    plt.show()
    return fig

def fitNflips (N, cM, fixAmpl = False, fixOff=False, show_plots = True):

    cM_m1 = squeeze(asarray(cM[0, :]))
    cM_0 = squeeze(asarray(cM[1, :]))
    cM_p1 = squeeze(asarray(cM[2, :]))

    A_guess = 1.
    off_guess= 0.
    p_guess=0.05
    
    fixP = []
    if fixAmpl == True:
        fixP.append(0)
    if fixOff == True:
        fixP.append(1)

    fit_m1 = fit.fit1d (asarray(N), cM_m1, common.fit_nuclearSpinFlips_1, A_guess, off_guess, p_guess, fixed = fixP,  do_plot = show_plots, ret = True)
    fit_0 = fit.fit1d (asarray(N), cM_0, common.fit_nuclearSpinFlips_2, A_guess, off_guess, p_guess,  fixed = fixP, do_plot = show_plots, ret=True)
    fit_p1 = fit.fit1d (asarray(N), cM_p1, common.fit_nuclearSpinFlips_3, A_guess, off_guess, p_guess,  fixed = fixP, do_plot = show_plots, ret=True)
    results = [fit_m1, fit_0, fit_p1]
    print fit_m1[0]['fitdata']
    plt.figure()
    plt.plot (N, fit_m1[0]['fitdata'], 'b', label = 'mI=-1')
    plt.plot (N, fit_0[0]['fitdata'], 'r', label = 'mI=0')
    plt.plot (N, fit_p1[0]['fitdata'], 'k', label = 'mI=+1')
    plt.plot (N, cM_m1, '>b')
    plt.plot (N, cM_0, 'ro')
    plt.plot (N, cM_p1, 'kD')
    plt.rcParams.update({'font.size': 18})
    plt.xlabel ('number of read-out steps', fontsize = 16)
    plt.ylabel ('readout corrected nuclear spin fidelity', fontsize = 16)   
    plt.legend (loc =1, prop={'size':12})
    plt.show()

    return results


def correction_pseudoinverse (N, matrix3, corrMatrix):

    #Here we force the resulting quantum state to be such that |a|^2+|b|^2+|c|^2=1,
    #basically getting a system ofg 3 equations in 3 unknowns, that we solve by pseudo-inverse
    #(which is equivalent to a least-squares minimization)
    
    M = []
    for i in arange(len(N)):
        a0 = matrix3[0,i] - corrMatrix[2,0]
        a1 = matrix3[1,i] - corrMatrix[2,1]
        a2 = matrix3[2,i] - corrMatrix[2,2]
        b0 = corrMatrix[0,0] - corrMatrix[2,0]
        b1 = corrMatrix[0,1] - corrMatrix[2,1]
        b2 = corrMatrix[0,2] - corrMatrix[2,2]
        c0 = corrMatrix[1,0] - corrMatrix[2,0]
        c1 = corrMatrix[1,1] - corrMatrix[2,1]
        c2 = corrMatrix[1,2] - corrMatrix[2,2]
        A = matrix([[a0],[a1],[a2]])
        B = matrix([[b0, c0], [b1, c1], [b2, c2]])
        x = linalg.pinv(B)*A
        x = squeeze(asarray(x))
        x = hstack ([x, 1-float(x[0])-float(x[1])])
        M.append(x)
    #M = transpose(M, [0, 1], [1, 0])
    M = vstack(M)
    M = transpose(M)
    return M

def correction (N, matrix3, corrMatrix):

    #corrects the read-out values just inverting the matrix with the calibration data
    F = corrMatrix.I
    M = F*matrix3
    return M

def plot_total_probability (N, correctedMatrix):
    
    total = squeeze(asarray( correctedMatrix[0,:]+correctedMatrix[1,:]+correctedMatrix[2,:]))

    plt.figure()
    plt.plot (N, total, 'bo')
    plt.rcParams.update({'font.size': 18})
    plt.xlabel ('number of read-out steps', fontsize = 16)
    plt.ylabel ('readout corrected nuclear spin fidelity', fontsize = 16)   
    plt.ylim ([0, 1.5])
    plt.show()
    


#calibration table for 11-12-2012
#U0 = matrix ([[0.09, 0.724, 0.692],[0.756, 0.06, 0.71],[0.734, 0.72, 0.076]]) #non-corrected
#U = matrix([[0.1017, 0.9283, 0.8866],[0.97, 0.0626,0.91],[0.9413, 0.92308,0.083]]) #electron SSRO corrected

#calibration table for 12-12-2012
#U0 = matrix([[1-0.117, 1-0.939, 1-0.963],[1-0.924, 1-0.127, 1-0.921],[1-0.887, 1-0.900, 1-0.065]])
#U0 = U0.T

#calibration table for 14-12-2012
#U0 = matrix([[1-0.088, 1-0.705, 1-0.710],[1-0.75, 1-0.10, 1-0.733],[1-0.72, 1-0.728, 1-0.082]])
#U0 = U0.T

# FIXME : I have the feeling that if you only put 'raw counts' in this matrix,
# you overcompensate for electron readout errors.
#
#calibration table for 16-12-2012
U0 = matrix([[1-0.083, 1-0.927, 1-0.938],[1-0.918, 1-0.082, 1-0.920],[1-0.91, 1-0.9, 1-0.068]])
#U0 = U0.T

#14-12-2012
#A2 and Ey

FN_m1 = r'D:\measuring\data\20121216\144632_MBI_pumping_A2_Ey_RO_mI-1_init-1'
FN_0 = r'D:\measuring\data\20121216\150333_MBI_pumping_A2_Ey_RO_mI0_init-1'
FN_p1 = r'D:\measuring\data\20121216\152128_MBI_pumping_A2_Ey_RO_mI+1_init-1'

#14-12-2012
#A2 and pi pulse
FN_m1 = r'D:\measuring\data\20121216\173627_MBI_pumping_A2_pi_pulse_RO_mI-1_init-1'
FN_0 = r'D:\measuring\data\20121216\172313_MBI_pumping_A2_pi_pulse_RO_mI0_init-1'
#FN_p1 = r'D:\measuring\data\20121216\164529_MBI_pumping_A2_pi_pulse_RO_mI+1_init-1'
FN_p1=r'D:\measuring\data\20121216\170904_MBI_pumping_A2_pi_pulse_RO_mI+1_init-1'
data_m1 = loadData(FN_m1+'\SSRO_readout.dat')
data_0 = loadData(FN_0+'\SSRO_readout.dat')
data_p1 = loadData(FN_p1+'\SSRO_readout.dat')

Nm1 = data_m1['Nreadouts']
N0 = data_0['Nreadouts']
Np1 = data_p1['Nreadouts']
rawRO_m1 = data_m1['corrRO']
rawRO_0 = data_0['corrRO']
rawRO_p1 = data_p1['corrRO']

uni = ones(len(Nm1))
rawRO = matrix([uni-asarray(rawRO_m1), uni-asarray(rawRO_0), uni-asarray(rawRO_p1)])
#corrRO = matrix([asarray(data_m1['corrRO']), asarray(data_0['corrRO']), asarray(data_p1['corrRO'])])
#corrNuclear = correction(Nm1, rawRO, U0)


plotNuclearState (Nm1,  correction(Nm1, rawRO, U0))

results = fitNflips (Nm1, correction(Nm1, rawRO, U0), fixAmpl= True, fixOff=False)
fit_m1 = results[0]
fit_0 = results[1]
fit_p1 = results[2]


#print "RO -1, paramters: " + str(fit_m1['params_dict'])
#print "RO -1, errors: "+ str(fit_m1['error_dict'])
#print
#print "RO 0, paramters: " + str(fit_0['params_dict'])
#print "RO 0, errors: "+ str(fit_0['error_dict'])
#print
#print "RO +1, paramters: " + str(fit_p1['params_dict'])
#print "RO +1, errors: "+ str(fit_p1['error_dict'])


fName = "C:\Documents and Settings\localadmin\My Documents\corrected.txt"
M = squeeze(asarray(correction (Nm1, rawRO, U0)))
savetxt(fName, M)


plot_total_probability (Nm1, correction(Nm1, rawRO, U0))





    

