import numpy as np
from matplotlib import pyplot as plt

import sympy
from sympy.utilities.lambdify import implemented_function, lambdify

class Formula:

    def __init__(self):
        self.formula = None
        self.values = {}
        self.uncertainties = {}

    def set_symbol_value(self, symbol, value):
        self.values[symbol] = value

    def symbol_value(self, symbol):
        return self.values[symbol]

    def set_symbol_uncertainty(self, symbol, value):
        self.uncertainties[symbol] = value

    def symbol_uncertainty(self, symbol):
        return self.uncertainties[symbol]

    def value(self):
        return self.formula.subs(self.values)
        
    def value_correlation(self):
        # The formula gives a sympy matrix, which lambdify cannot work with properly
        # We turn it into tuple here; lambdify accepts this.
        return tuple(self.formula.subs(self.values))

    def uncertainty_squared(self):
        # usquared = 0.
        for i,s in enumerate(self.uncertainties):
            if i == 0:
                usquared = (self.formula.diff(s).subs(self.values)*\
                    self.uncertainties[s])**2
            else:
                usquared += (self.formula.diff(s).subs(self.values)*\
                    self.uncertainties[s])**2
        
        return usquared

    def uncertainty_squared_correlation(self):
        #Since formula gives a sympy matrix, we cannot square this term by term with **.
        #We therefore square each term one by one, and thus create usquared..
        usquared = {}
        for i,s in enumerate(self.uncertainties):
            if i == 0:
                for j in range(4):
                    usquared[j] = (self.formula.diff(s).subs(self.values)[j]*\
                        self.uncertainties[s])**2
            else:
                for j in range(4):
                    usquared[j] += (self.formula.diff(s).subs(self.values)[j]*\
                        self.uncertainties[s])**2
        
        return usquared

    def uncertainty(self):
        return sympy.sqrt(self.uncertainty_squared())

    def num_eval(self, symbol=None, values=None, uncertainties=None):
        if symbol == None:
            return self.value(), self.uncertainty()
        
        valfunc = lambdify(symbol, self.value(), 'numpy')
        if uncertainties == None:        
            uncertaintyfunc = lambdify(symbol, self.uncertainty(), 'numpy')
            return valfunc(values), uncertaintyfunc(values)
        else:
            u_sym = sympy.symbols('u_sym')
            usquared = self.uncertainty_squared() + \
                    (self.formula.diff(symbol).subs(self.values) * u_sym)**2
            u = sympy.sqrt(usquared)
            uncertaintyfunc = lambdify((symbol, u_sym), u, 'numpy')
            return valfunc(values), uncertaintyfunc(values, uncertainties)
    
    def num_eval_correlation(self, symbol=None, values=None, uncertainties=None):
        if symbol == None:
            return self.value(), self.uncertainty()
        
        valfunc = lambdify(symbol, self.value_correlation(), 'numpy')
        if uncertainties == None:        
            uncertaintyfunc = lambdify(symbol, self.uncertainty(), 'numpy')
            return valfunc(values), uncertaintyfunc(values)
        else:
            #We define u_sym as a 'DeferredVector' to be able to use it in lambdify as one input.
            u_sym = sympy.DeferredVector('u_sym')
            #To be able to square term by term in a sympy matrix, square terms separately:
            usquared = {}
            for j in range(4):
                usquared[j] = self.uncertainty_squared_correlation()[j] + \
                    (self.formula.diff(symbol).subs(self.values)[j] * u_sym[j])**2
            u = (sympy.sqrt(usquared[0]),sympy.sqrt(usquared[1]),sympy.sqrt(usquared[2]),sympy.sqrt(usquared[3]))
            uncertaintyfunc = lambdify((symbol, u_sym), u, 'numpy')
            return valfunc(values), uncertaintyfunc(values, uncertainties)       

class SingleQubitROC:

    def __init__(self):
        self.F0 = 1.
        self.F1 = 1.
        self.u_F0 = 0.
        self.u_F1 = 0.        

    def _setup(self):
        self._F0, self._F1 = sympy.symbols('F0 F1')
        self._p0 = sympy.symbols('p0')

        self.error_matrix = sympy.Matrix(
                [[self._F0, 1-self._F1],
                    [1-self._F0, self._F1]])

        self.correction_matrix = self.error_matrix.inv()
        
        corr_vec = self.correction_matrix * sympy.Matrix([self._p0, 1-self._p0])
        corr_p0 = corr_vec[0]

        self.p0_formula = Formula()
        self.p0_formula.formula = corr_p0
    
    
    def num_eval(self, p0, u_p0=None):
        self._setup()

        self.p0_formula.values[self._F0] = self.F0
        self.p0_formula.values[self._F1] = self.F1
        self.p0_formula.uncertainties[self._F0] = self.u_F0
        self.p0_formula.uncertainties[self._F1] = self.u_F1
        
        return self.p0_formula.num_eval(self._p0, p0, u_p0)


def simple_demo():
    f = Formula()
    a,b,x = sympy.symbols('a,b,x')
    f.formula = a*x + b
    
    f.values[a] = 2
    f.values[b] = 1
    f.uncertainties[a] = 0.2

    xvals = np.arange(10)
    yvals,yerrs = f.num_eval(x, xvals, xvals*0.1)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(xvals, yvals, yerr=yerrs, fmt='o')

def roc_demo():
    roc = SingleQubitROC()
   
    roc.F0 = 0.8
    roc.F1 = 0.9
    roc.u_F0 = 0.001
    roc.u_F1 = 0.001

    xvals = np.linspace(0.1, 0.8, 21)
    xerrs = xvals * 0.05
    yvals, yerrs = roc.num_eval(xvals, xerrs)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(xvals, yvals, yerr=yerrs, fmt='o')


if __name__ == '__main__':
    roc_demo()

