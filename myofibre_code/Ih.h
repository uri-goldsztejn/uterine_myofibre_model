#pragma once
#include "Current.h"
//from Testrow et al.
class Ih : public Current {
private:
	double y = 0.002604864867063448;		/* dynamical variable for Ih activation */
	double gh = 0.0542; // 0.04065;	/* Maximal conductance for Ih (nS/pF)  */
	double PK = 1.0;	/* Permeability ratios for K+ */
	double PNa = 0.35; /* Permeability ratios for Na+ */

	double yss(double v) {
		return 1.0 / (1.0 + exp((v + 105.39) / 8.6553));	/* steady-state for y*/
	}

	double ya(double v) {
		return 3.5e-6*exp(-0.0497*v);	
	}

	double yb(double v) {
		return 0.04003*exp(0.05211*v);
	}
	double ytc(double v) {
		return 1.0 / (ya(v) + yb(v)); 		/* time constant for y */
	}


public:

	//returns the current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		double  eh = (R*temp / frdy)*log((ko + (PNa / PK)*nao) / (parameters.ki + (PNa / PK)*parameters.nai)); /* reversal potentials (mV) */
		return gh*y*(parameters.getVoltage() - eh); /* Ih*/
	}

	// updates the ion channel's gates for the following iteration
	void updateGates(CellParameters &parameters) {
		y += dt*((yss(parameters.getVoltage()) - y) / (ytc(parameters.getVoltage())));
	}
};