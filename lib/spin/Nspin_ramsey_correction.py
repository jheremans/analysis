import sympy
from analysis.lib.math import error

class NuclearSpinROC(error.SingleQubitROC):

    def __init__(self):
        error.SingleQubitROC.__init__(self)

        self.F0_ssro = 1.
        self.u_F0_ssro = 0.
        self.F1_ssro = 1.
        self.u_F1_ssro = 0.
        self.F_init = 1.
        self.u_F_init = 0.
        self.F_RO_pulse = 1.
        self.u_F_RO_pulse = 0.
    
    def _setup(self):
        error.SingleQubitROC._setup(self)

        self.F0_formula = error.Formula()
        self._F0_ssro, self._F1_ssro, self._F_init, self._F_RO_pulse = \
                sympy.symbols('F0_ssro F1_ssro F_init F_RO_pulse')
        
        self.F0_formula.formula = self._F0_ssro * self._F_RO_pulse / self._F_init
        
        self.F1 = self.F1_ssro
        self.u_F1 = self.u_F1_ssro

        self.F0_formula.values[self._F0_ssro] = self.F0_ssro
        self.F0_formula.uncertainties[self._F0_ssro] = self.u_F0_ssro
        self.F0_formula.values[self._F_init] = self.F_init
        self.F0_formula.uncertainties[self._F_init] = self.u_F_init
        self.F0_formula.values[self._F_RO_pulse] = self.F_RO_pulse
        self.F0_formula.uncertainties[self._F_RO_pulse] = self.u_F_RO_pulse

        self.F0 = float(self.F0_formula.value())
        self.u_F0 = float(self.F0_formula.uncertainty())

    #def num_eval(self, *arg, **kw):
    #    return error.SingleQubitROC.num_eval(self, *arg, **kw)
