#pragma once
#include "Current.h"

class IKa : public Current {
private:
	double s = 0.0307583106982354;		/* dynamical variable for IKa activation */
	double x = 0.08785242843398365;		/* dynamical variable for IKa inactivation  */
	double gka = 0.2;// 0.001418;		/* Conductance ratios between ik1, ik2, ika (dimensionless) */
	double gk = 0.8;		/* Maximal conductance for total voltage-gated K currents (IK1, IK2, IKa) (nS/pF)  */
	double modulation = 1.0; // conductance modulation

	double sss(double v) {
		return 1.0 / (1.0 + exp(-(v + 27.79) / 7.57));	/* steady-state for s*/
	}

	double xss(double v) {
		return 0.02 + 0.98 / (1.0 + exp((v + 69.5) / 6.0));	/* steady-state for x*/
	}

	double stc(double v) {
		return 17.0 / (1.0 + ((v + 20.5232) / 35.0)*((v + 20.5232) / 35.0)); 		/* time constant for s */
	}

	double xtc(double v) {
		return 7.5 + 10.0 / (1.0 + ((v + 34.1765) / 120.0)*((v + 34.1765) / 120.0));	/* time constant for x */
	}


public:

	// initializes the modulation to the user-defined value 
	void setModulation(double mod) {
		modulation = mod;
	}

	//returns the channel current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return modulation*gk*gka*s*x*(parameters.getVoltage() - parameters.ek()); /* IKa*/
	}

	// updates the ion channel gates for the following iteration
	void updateGates(CellParameters &parameters) {
		double v = parameters.getVoltage();
		s += dt*((sss(v) - s) / (stc(v)));
		x += dt*((xss(v) - x) / (xtc(v)));
	}
};