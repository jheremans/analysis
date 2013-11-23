"""
Module that helps reconstructing the density matrix from measurement
data and calculating fidelities.
"""

import numpy as np
import sympy
import error
from analysis.lib.math import error
from sympy.matrices import Matrix

def dm(statevec):
    '''
    returns the density matrix for the input state vector.
    '''
    return np.outer(statevec, statevec.conj())

def measured_single_qubit_dm(mx, my, mz, u_mx=0., u_my=0., u_mz=0.):
    '''
    returns the density matrix for measurement results on a single qubit
    in the x,y,z bases. we take the measurement result to be the probability
    to measure the qubit state |0>.
    note: measuring in the x-basis is achieved by a +pi/2 rotation around -y,
    a measurement in the y basis by a +pi/2 rotation around +x.
    returns the density matrix and the uncertainty thereof.
    '''
    a = mz; u_a = u_mz
    x = mx - .5; u_x = u_mx
    y = .5 - my; u_y = u_my

    return (np.array([[a, x+1j*y],[x-1j*y, 1.-a]]), 
            np.array([[u_a, u_x+1j*u_y],[u_x-1j*u_y, u_a]])) #1-u_a

def fidelity(statevec, dm, u_dm, ustatevec=np.array([0,0])):
    '''
    returns the fidelity (and its uncertainty) of the measured density
    matrix with a given state vector.
    '''
   
    f = error.Formula()
    beta,a,x,b,alpha = sympy.symbols('beta,a,x,b,alpha')
    v = Matrix([alpha, beta])
    rho = Matrix([[x,a+1j*b],[a-1j*b, 1-x]])    
    f.formula = (v.conjugate().transpose() * rho * v)[0]
    
    f.values[alpha] = statevec[0]
    f.values[beta] = statevec[1]
    f.values[a]=float(np.real(dm[0,1]))
    f.values[x]=float(dm[0,0])
    f.values[b]=float(np.imag(dm[0,1]))

    f.uncertainties[alpha]=ustatevec[0]
    f.uncertainties[beta]=ustatevec[1]
    f.uncertainties[x]=u_dm[0,0]
    f.uncertainties[a]=float(np.real(u_dm[0,1]))
    f.uncertainties[b]=float(np.imag(u_dm[0,1]))
    
    _fid,_ufid = f.num_eval()
    fid = float(_fid.as_real_imag()[0])
    ufid = float(_ufid.as_real_imag()[0])
   
    return (fid,ufid)

