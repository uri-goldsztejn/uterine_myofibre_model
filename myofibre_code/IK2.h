#pragma once
#include "Current.h"
//from Testrow et al.
class IK2 : public Current {
private:
	
	double p = 0.1174074734567931;		/* dynamical variable for IK2 activation */
	double k1 = 0.9968385770271651;	/* dynamical variable for IK2 fast inactivation  */
	double k2 = 0.9968408069904307;		/* dynamical variable for IK2 slow inactivation */
	double gk2 = 0.04;// 0.0315;		/* Conductance ratios between ik1, ik2, ika (dimensionless) */
	double gk = 0.8;		/* Maximal conductance for total voltage-gated K currents (IK1, IK2, IKa) (nS/pF)  */
	double modulation = 1.0;	// conductance modulation

	double pss(double v) {
		return 0.948 / (1.0 + exp(-(v + 17.91) / 18.4));	/* steady-state for p*/
	}

	double k1ss(double v) {
		return 1.0 / (1.0 + exp((v + 21.2) / 5.7));		/* steady-state for k1*/
	}

	double k2ss(double v) {
		return k1ss(v);								/* steady-state for k2*/
	}

	double ptc(double v) {
		return 100.0 / (1.0 + ((v + 64.1) / 28.67)*((v + 64.1) / 28.67));	/* time constant for p */
	}

	double k1tc(double v) {
		return 1.0e6*(1.0 - (1.0 / (1.0 + exp((v - 315.0) / 50.0)))*(1.0 / (1.0 + exp(-(v + 74.9) / 8.0)))); /* time constant for k1 */
	}

	double k2tc(double v) {
		return 1000.0*2500.0*(1.0 - (1.0 / (1.0 + exp((v - 132.868) / 25.3992)))*(1.0 / (1.0 + exp(-(v + 24.9203) / 2.67915))));	/* time constant for k2 */
	}


	
public:

	// initializes the modulation to the user-defined value 
	void setModulation(double mod) {
		modulation = mod;
	}

	//returns the current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return modulation*gk*gk2*(p*p)*(0.75*k1 + 0.25*k2)*(parameters.getVoltage() - parameters.ek()); /* IK2 */
	}

	// updates the ion channel's gates for the following iteration
	void updateGates(CellParameters &parameters) {
		double v = parameters.getVoltage();
		p += dt*((pss(v) - p) / (ptc(v)));
		k1 += dt*((k1ss(v) - k1) / (k1tc(v)));
		k2 += dt*((k2ss(v) - k2) / (k2tc(v)));
	}
};