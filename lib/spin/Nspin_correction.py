import sympy
from sympy.utilities.lambdify import implemented_function, lambdify
from analysis.lib.math import error

class NuclearSpinROC:
    def __init__(self):
        self.F0_ssro = 1.
        self.u_F0_ssro = 0.
        self.F1_ssro = 1.
        self.u_F1_ssro = 0.
        
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
        self._F0_ssro, self._F1_ssro, self._F0_RO_pulse, self._F1_RO_pulse, self._P_min1, self._P_0 = \
           sympy.symbols('F0_ssro F1_ssro F0_RO_pulse F1_RO_pulse P_min1 P_0')
        self._p0 = sympy.symbols('p0')
               
        # The error matrix is built up of three matrices RO_err, CNOT_err and 
        # Init_err. The first reflects fidelities of the electron RO, 
        # the second the fidelities in the pi (F0) and 2pi (F1) parts of the 
        # pi-2pi pulse (as CNOT gate). The third matrix includes populations 
        # in the three nitrogen lines after initialization in -1.
        self.error_matrix =sympy.Matrix(
                                        [[self._F0_ssro, 1.-self._F1_ssro],
                                         [1.-self._F0_ssro, self._F1_ssro]])* \
                                sympy.Matrix([[self._F0_RO_pulse, 1.-self._F1_RO_pulse, 0.],
                                              [1.-self._F0_RO_pulse, self._F1_RO_pulse, 1.]])* \
                                    sympy.Matrix([[self._P_min1,self._P_0],
                                                  [self._P_0,self._P_min1],
                                                  [1.-self._P_0-self._P_min1,1.-self._P_0-self._P_min1]])
        
        self.correction_matrix = self.error_matrix.inv()
                                
        corr_vec =  self.correction_matrix * sympy.Matrix([self._p0, 1-self._p0])
        corr_p0 = corr_vec[0]
        
        self.p0_formula = error.Formula()
        self.p0_formula.formula = corr_p0
                      
    def num_evaluation(self, p0, u_p0=None):
        self._setup()
        
        self.p0_formula.values[self._F0_ssro] = self.F0_ssro
        self.p0_formula.uncertainties[self._F0_ssro] = self.u_F0_ssro
        self.p0_formula.values[self._F1_ssro] = self.F1_ssro
        self.p0_formula.uncertainties[self._F1_ssro] = self.u_F1_ssro
        
        self.p0_formula.values[self._F0_RO_pulse] = self.F0_RO_pulse
        self.p0_formula.uncertainties[self._F0_RO_pulse] = self.u_F0_RO_pulse
        self.p0_formula.values[self._F1_RO_pulse] = self.F1_RO_pulse
        self.p0_formula.uncertainties[self._F1_RO_pulse] = self.u_F1_RO_pulse

        self.p0_formula.values[self._P_min1] = self.P_min1
        self.p0_formula.uncertainties[self._P_min1] = self.u_P_min1
        self.p0_formula.values[self._P_0] = self.P_0
        self.p0_formula.uncertainties[self._P_0] = self.u_P_0
                
        self.p0_formula.num_eval(self._p0, p0, u_p0)

        return self.p0_formula.num_eval(self._p0, p0, u_p0)
    
    