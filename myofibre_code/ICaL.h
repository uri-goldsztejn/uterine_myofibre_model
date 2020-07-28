#pragma once
//from Testrow et al.
#include "Current.h"

class ICaL: public Current {
private:
    double modulation = 1.0;	// conductance modulation
	double d = 0.01036961357784695;		/* dynamical variable for ICaL activation */
	double f1 = 0.9065941499695301;	/* dynamical variable for ICaL fast inactivation */
	double f2 = 0.9065967263076083;		/* dynamical variable for ICaL slow inactivation */
	double gcal = 0.6;// 0.318; 		/* Maximal conductance for ICaL (nS/pF)  */
	double ecal = 45.0;		/* reversal potential for ICaL (mV) */
	double kmca = 0.001;// 0.0006;	/*  half saturation concentration for ICaL Ca inhibition (mM) */

	double fca(double cai) {
		return 1.0 / (1.0 + (cai / kmca)); /* Ca inhibition in ICaL */
	}

	double dss(double v) {
		return 1.0 / (1.0 + exp(-(v + 22.0) / 7.0)); /* steady-state for d*/
	}

	double f1ss(double v) {
		return 1.0 / (1.0 + exp((v + 38.0) / 7.0));	/* steady-state for f1*/
	}

	double f2ss(double v) {
		return f1ss(v);				/* steady-state for f2*/
	}
	double dtc(double v) {
		return 2.29 + 5.7 / (1.0 + ((v + 29.97) / 9.0)*((v + 29.97) / 9.0)); /* time constant for d */
	}

	double f1tc(double v) {
		return 12.0;		/* time constant for f1 */
	}

	double f2tc(double v) {
		return 90.9699*(1.0 - (1.0 / (1.0 + exp((v + 13.9629) / 45.3782)))*(1.0 / (1.0 + exp(-(v + 9.49866) / 3.3945))));	/* time constant for f2 */
	}



public:

	// initializes the modulation to the user-defined value 
	void setModulation(double mod) {
		modulation = mod;
	}

	//returns the current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return modulation*gcal*fca(parameters.getCai())*(d*d)*(0.8*f1 + 0.2*f2)*(parameters.getVoltage() - ecal);	/* ICaL  */
	}

	// updates the ion channel's gates for the following iteration
	void updateGates(CellParameters &parameters) {
		double v = parameters.getVoltage();
		d += dt*((dss(v) - d) / (dtc(v)));
		f1 += dt*((f1ss(v) - f1) / (f1tc(v)));
		f2 += dt*((f2ss(v) - f2) / (f2tc(v)));
	}
};