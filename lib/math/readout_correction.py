import numpy as np
import error

def single_qubit_outcome_with_ROC(zero_events, all_events, F0, u_F0, F1, u_F1):
	'''
	Calculate the readout corrected result for P(|0>).
	expects the number of 0-events, and the total number of events, as well
	as the readout fidelities for 0 and 1.

	returns P(0), u_P(0)
	'''
	if type(zero_events) == int or type(zero_events) == float:
		frac0 = float(zero_events)/all_events
		u_frac0 = np.sqrt(frac0*(1.-frac0)/all_events)
	else:
		frac0 = zero_events.astype(float)/float(all_events)
		u_frac0 = np.sqrt(frac0*(1.-frac0)/all_events)

	roc = error.SingleQubitROC()
	roc.F0, roc.u_F0, roc.F1, roc.u_F1 = F0, u_F0, F1, u_F1

	P0, u_P0 = roc.num_eval(frac0, u_frac0)

	return P0, u_P0

def single_qubit_outcome_with_ROC_from_fraction(frac0, u_frac0, F0, u_F0, F1, u_F1):
	'''
	Calculate the readout corrected result for P(|0>).
	expects the fraction and uncertanty of 0-events, as well
	as the readout fidelities for 0 and 1.

	returns P(0), u_P(0)
	'''
	roc = error.SingleQubitROC()
	roc.F0, roc.u_F0, roc.F1, roc.u_F1 = F0, u_F0, F1, u_F1

	P0, u_P0 = roc.num_eval(frac0, u_frac0)

	return P0, u_P0

