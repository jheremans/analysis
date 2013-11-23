import sympy
from sympy.utilities.lambdify import implemented_function, lambdify
from sympy import DeferredVector
from sympy.physics.quantum import TensorProduct
from analysis.lib.math import error
reload(error)

class CorrelationsROC:
    def __init__(self):
        self.F0_e_ssro = 1.
        self.u_F0_e_ssro = 0.
        self.F1_e_ssro = 1.
        self.u_F1_e_ssro = 0.
        
        self.F0_N_ssro = 1.
        self.u_F0_N_ssro = 0.
        self.F1_N_ssro = 1.
        self.u_F1_N_ssro = 0.
        
        self.P_min1 = 0.95
        self.u_P_min1 = 0
        self.P_0 = 0.04
        self.u_P_0 = 0
        self.P_plus1 = 0.01
        self.u_P_plus1 = 0
        self.F0_RO_pulse = 0.98
        self.u_F0_RO_pulse = 0.0105
        self.F1_RO_pulse = 0.99
        self.u_F1_RO_pulse = 0.01
    
    def _setup(self):
        self._F0_N_ssro, self._F1_N_ssro, self._F0_e_ssro, self._F1_e_ssro, self._F0_RO_pulse, self._F1_RO_pulse, self._P_min1, self._P_0 = \
            sympy.symbols('F0_N_ssro F1_N_ssro F0_e_ssro F1_e_ssro F0_RO_pulse F1_RO_pulse P_min1 P_0')
        self._p_correlations = sympy.DeferredVector('p_correlations')
               
        # The total error matrix is the tensor product of the nitrogen and electron error matrices.
        # The nitrogen error matrix is built up of three matrices RO_err, CNOT_err and Init_err. 
        # RO_err reflects fidelities of the electron RO, 
        # CNOT_err reflects the fidelities in the pi (F0) and 2pi (F1) parts of the pi-2pi pulse (as CNOT gate). 
        # Init_err includes populations in the three nitrogen lines after initialization in -1.
        #The inverse matrices are calculated before the tensorproduct, since this is much faster.
        
        self.error_matrix_N = (sympy.Matrix([[self._F0_N_ssro, 1.-self._F1_N_ssro],[1.-self._F0_N_ssro, self._F1_N_ssro]]) * \
            sympy.Matrix([[self._F0_RO_pulse, 1.-self._F1_RO_pulse, 0.],[1.-self._F0_RO_pulse, self._F1_RO_pulse, 1.]]) * \
                sympy.Matrix([[self._P_min1,self._P_0],[self._P_0,self._P_min1],[1.-self._P_0-self._P_min1,1.-self._P_0-self._P_min1]]))
 
        self.error_matrix_e = (sympy.Matrix([[self._F0_e_ssro, 1.-self._F1_e_ssro],[1.-self._F0_e_ssro, self._F1_e_ssro]]))
        
        self.correction_matrix_N = self.error_matrix_N.inv()
        self.correction_matrix_e = self.error_matrix_e.inv()  
        self.correction_matrix = TensorProduct(self.correction_matrix_N,self.correction_matrix_e)

        corr_vec =  self.correction_matrix * \
            sympy.Matrix([self._p_correlations[0], self._p_correlations[1], self._p_correlations[2], self._p_correlations[3]])
        corr_p_correlations = corr_vec
                
        self.p0_formula = error.Formula()
        self.p0_formula.formula = corr_p_correlations
                      
    def num_evaluation(self, p_correlations, u_p_correlations=None):
        self._setup()
        
        #Here all values and uncertainties of the fidelities are listed.
        self.p0_formula.values[self._F0_e_ssro] = self.F0_e_ssro
        self.p0_formula.uncertainties[self._F0_e_ssro] = self.u_F0_e_ssro
        self.p0_formula.values[self._F1_e_ssro] = self.F1_e_ssro
        self.p0_formula.uncertainties[self._F1_e_ssro] = self.u_F1_e_ssro

        self.p0_formula.values[self._F0_N_ssro] = self.F0_N_ssro
        self.p0_formula.uncertainties[self._F0_N_ssro] = self.u_F0_N_ssro
        self.p0_formula.values[self._F1_N_ssro] = self.F1_N_ssro
        self.p0_formula.uncertainties[self._F1_N_ssro] = self.u_F1_N_ssro
                
        self.p0_formula.values[self._F0_RO_pulse] = self.F0_RO_pulse
        self.p0_formula.uncertainties[self._F0_RO_pulse] = self.u_F0_RO_pulse
        self.p0_formula.values[self._F1_RO_pulse] = self.F1_RO_pulse
        self.p0_formula.uncertainties[self._F1_RO_pulse] = self.u_F1_RO_pulse

        self.p0_formula.values[self._P_min1] = self.P_min1
        self.p0_formula.uncertainties[self._P_min1] = self.u_P_min1
        self.p0_formula.values[self._P_0] = self.P_0
        self.p0_formula.uncertainties[self._P_0] = self.u_P_0

        return self.p0_formula.num_eval_correlation(self._p_correlations, p_correlations.T , u_p_correlations.T)
    
