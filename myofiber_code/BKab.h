#pragma once
#include "Current.h"

class BKab : public Current {
private:
	double xab = 0.002220456569762898;		/* dynamical variable for BKab */
	double gbkab = 0.1;		/* Conductance ratios between BKa, BKab (dimensionless) */
	double modulation = 1.0;// conductance modulation

	double xabss_z(double cai) {
		return -0.681249 / (1.0 + ((cai*1000.0 - 0.218988) / 0.428335)*((cai*1000.0 - 0.218988) / 0.428335)) + 1.40001 / (1.0 + ((cai*1000.0 + 228.71) / 684.946)*((cai*1000.0 + 228.71) / 684.946)); /* gating charge for xab*/
	}

	double xabss_vh(double cai) {
		return 8540.23 / (1.0 + pow(((cai*1000.0 + 0.401189) / 0.00399115), 0.668054)) - 109.275;	/* half-activation constant for xab*/
	}

	double xabss(double v, double cai) {
		return 1.0 / (1.0 + exp(-xabss_z(cai)*frdy*(v - xabss_vh(cai)) / (R*temp)));		/* steady-state for xab */
	}

	double xabtc(double v) {
		return 13.8049 / (1.0 + ((v - 153.019) / 66.4952)*((v - 153.019) / 66.4952));		/* time constant for xab */
	}


	
public:
	// initializes the modulation to the user-defined value 

	void setModulation(double mod) {
		modulation = mod;
	}

	//returns the channel current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return modulation*gkca*gbkab*xab*(parameters.getVoltage() - parameters.ek()); /* BKab */
	}

	// updates the ion channel gates for the following iteration
	void updateGates(CellParameters &parameters) {
		xab += dt*((xabss(parameters.getVoltage(), parameters.getCai()) - xab) / (xabtc(parameters.getVoltage())));
	}
};