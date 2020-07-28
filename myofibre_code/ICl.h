#pragma once
#include "Current.h"

class ICl : public Current {
private:
	double c = 0.0003764413740731269;		/* dynamical variable for ICl activation */
	double gcl = 0.1875;	/* Maximal conductance for ICl (nS/pF)  */
	double modulation = 1.0; // conductance modulation

	double vFRT(double v) {
		return (v*frdy) / (R*temp);
	}

	double K1cl(double v) {
		return 0.0006*exp(2.53*vFRT(v));
	}
	double K2cl(double v) {
		return 0.1*exp(-5.0*vFRT(v));
	}
	double css(double v, double cai) {
		return 1.0 / (1.0 + K2cl(v)*((K1cl(v) / cai)*(K1cl(v) / cai) + K1cl(v) / cai + 1.0));		/* steady-state for c*/
	}

	double ctc(double v) {
		return -160.0 + 210.0 / (1.0 + exp((v + 4.56) / 11.62)) + 170.0 / (1.0 + exp(-(v + 25.5) / 11.62));	/* time constant for c */
	}

public:

	// initializes the modulation to the user-defined value 
	void setModulation(double mod) {
		modulation = mod;
	}

	//returns the current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return modulation*gcl*c*(parameters.getVoltage() - parameters.ecl()); /* ICl */
	}
	
	// updates the ion channel's gates for the following iteration
	void updateGates(CellParameters &parameters) {
		c += dt*((css(parameters.getVoltage(), parameters.getCai()) - c) / (ctc(parameters.getVoltage())));

	}
};