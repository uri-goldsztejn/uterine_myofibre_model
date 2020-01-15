#pragma once
//from Testrow et al.
#include "Current.h"

class ICaT : public Current {
private:
	double b = 0.508117603077852;	/* dynamical variable for ICaT activation  */
	double g = 0.03582573962705717;		/* dynamical variable for ICaT inactivation */
	double gcat = 0.058;// 0.02319;	/* Maximal conductance for ICaT (nS/pF)  */
	double ecat = 42.0;		/* reversal potential for ICaT (mV) */
	double modulation = 1.0;	// conductance modulation

	double bss(double v) {
		return 1.0 / (1.0 + exp(-(v + 54.23) / 9.88)); /* steady-state for b*/
	}

	double gss(double v) {
		return 0.02 + (1.0 - 0.02) / (1.0 + exp((v + 72.978) / 4.64)); 	/* steady-state for g*/
	}

	double btc(double v) {
		return 0.45 + 3.9 / (1.0 + ((v + 66.0) / 26.0)*((v + 66.0) / 26.0)); 	/* steady-state for b*/
	}

	double gtc(double v) {
		return 150.0*(1.0 - (1.0 / (1.0 + exp((v - 417.43) / 203.18)))*(1.0 / (1.0 + exp(-(v + 61.11) / 8.07)))); 	/* steady-state for g*/
	}

	
public:
	// initializes the modulation to the user-defined value
	void setModulation(double mod) {
		modulation = mod;
	}
	//returns the channel current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return modulation*gcat*(b*b)*g*(parameters.getVoltage() - ecat); /* ICaT */
	}
	// updates the ion channel gates for the following iteration
	void updateGates(CellParameters &parameters) {
		double v = parameters.getVoltage();
		b += dt*((bss(v) - b) / (btc(v)));
		g += dt*((gss(v) - g) / (gtc(v)));
	}
};