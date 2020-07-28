#pragma once
#include "Current.h"
//from Testrow et al.
class IK1 : public Current {
private:

	double q = 0.2060363247740295;	/* dynamical variable for IK1 activation  */
	double r1 = 0.1922244113609531;	/* dynamical variable for IK1 fast inactivation  */
	double r2 = 0.1932803618375963;	/* dynamical variable for IK1 slow inactivation  */
	double gk1 = 0.65;				/* Conductance ratios between ik1, ik2, ika (dimensionless) */
	double gk = 0.8;		/* Maximal conductance for total voltage-gated K currents (IK1, IK2, IKa) (nS/pF)  */
	double modulation = 1.0; // conductance modulation

	double qss(double v) {
		return 0.978613 / (1.0 + exp(-(v + 18.6736) / 26.6606)); /* steady-state for q*/
	}
	double r1ss(double v) {
		return 1.0 / (1.0 + exp((v + 63.0) / 6.3)); /* steady-state for r1*/
	}

	double r2ss(double v) {
		return r1ss(v);								/* steady-state for r2*/
	}

	double qtc(double v) {
		return 500.0 / (1.0 + ((v + 60.71) / 15.79)*((v + 60.71) / 15.79));	/* time constant for q */
	}

	double r1tc(double v) {
		return 5000.0 / (1.0 + ((v + 62.7133) / 35.8611)*((v + 62.7133) / 35.8611));/* time constant for r1 */
	}

	double r2tc(double v) {
		return 30000.0 + 220000.0 / (1.0 + exp((v + 22.0) / 4.0));			/* time constant for r2 */
	}

public:

	// initializes the modulation to the user-defined value 
	void setModulation(double mod) {
		modulation = mod;
	}

	//returns the current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return modulation*gk*gk1*(q*q)*(0.38*r1 + 0.63*r2)*(parameters.getVoltage() - parameters.ek()); /* IK1 */
	}

	// updates the ion channel's gates for the following iteration
	void updateGates(CellParameters &parameters) {
		double v = parameters.getVoltage();
		q += dt*((qss(v) - q) / (qtc(v)));
		r1 += dt*((r1ss(v) - r1) / (r1tc(v)));
		r2 += dt*((r2ss(v) - r2) / (r2tc(v)));
	}
};