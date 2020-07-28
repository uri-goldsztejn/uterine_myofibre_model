#pragma once
#include "Current.h"
//from Testrow et al.
class Ina: public Current{
private:
	double m = 0.1253518889572223;		/* dynamical variable for INa activation */
	double h = 0.404599170710196;	/* dynamical variable for INa inactivation */
	double gna =  0.1152/2; 		/* 0-0.12; maximal conductance for INa (nS/pF)  (In tong it's 0, in the new one it's 0.1152)*/
	//double ena = (R*temp / frdy)*log(nao / nai);	/* reversal potential for sodium (mV) */
	double modulation = 1.0;	// conductance modulation

	double mss(double v) {
		return 1.0 / (1.0 + exp(-(v + 35.9584) / 9.24013)); /* steady-state function for m gate */
	}

	double hss(double v) {
		return 1.0 / (1.0 + exp((v + 57.0) / 8.0));			/* steady-state function for h gate*/
	}

	double mtc(double v) {
		return 0.25 + 7.0 / (1.0 + exp((v + 38.0) / 10.0));		/* time constant function for m gate */
	}
	
	double htc(double v) {
		return 0.9 + 1002.85 / (1.0 + ((v + 47.5) / 1.5)*((v + 47.5) / 1.5)); /* time constant function for h gate */
	}




public:
	// initializes the modulation to the user-defined value
	void setModulation(double mod) {
		modulation = mod;
	}
	//returns the current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return modulation*gna*(m*m*m)*h*(parameters.getVoltage() - parameters.ena()); /* INa */
	}
	// updates the ion channel's gates for the following iteration
	void updateGates(CellParameters &parameters) {
		double v = parameters.getVoltage();
		m += dt*((mss(v) - m) / (mtc(v)));
		h += dt*((hss(v) - h) / (htc(v)));
	}
};