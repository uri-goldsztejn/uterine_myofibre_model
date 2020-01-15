#pragma once
#include "Current.h"

class BKa : public Current {
private:
	double xa = 0.0003569126518797985;		// dynamical variable for BKa 
	double gbka = 0.2;		// Conductance ratios between BKa, BKab (dimensionless) 
	double modulation = 1.0; // conductance modulation

	double xass_z(double cai) {
		return -0.749234 / (1.0 + ((cai*1000.0 - 0.0630535) / 0.161942)*((cai*1000.0 - 0.0630535) / 0.161942)) + 8.38384 / (1.0 + ((cai*1000.0 + 1538.29) / 739.057)*((cai*1000.0 + 1538.29) / 739.057));	/* gating charge for xa*/
	}

	double xass_vh(double cai) {
		return 5011.47 / (1.0 + pow(((cai*1000.0 + 0.237503) / 0.000239278), 0.42291)) - 37.5137;	/* half-activation constant for xa*/
	}

	double xass(double v, double cai) {
		return 1.0 / (1.0 + exp(-xass_z(cai)*frdy*(v - xass_vh(cai)) / (R*temp)));			/* steady-state for xa*/
	}

	double xatc(double v) {
		return 2.40914 / (1.0 + ((v - 158.779) / (-52.1497))*((v - 158.779) / (-52.1497)));			/* time constant for xa */
	}

public:
	// initializes the modulation to the user-defined value 
	void setModulation(double mod) {
		modulation = mod;
	}

	//returns the channel current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return modulation*gkca*gbka*xa*(parameters.getVoltage() - parameters.ek()); /* BKa */
	}
	// updates the ion channel gates for the following iteration
	void updateGates(CellParameters &parameters) {
		xa += dt*((xass(parameters.getVoltage(), parameters.getCai()) - xa) / (xatc(parameters.getVoltage())));

	}
};