#pragma once
#include "Current.h"

class INSCC : public Current {
private:
	double gns = 0.0123;	/* Maximal conductance for INSCC (nS/pF)  */
	double gl = 0.0;// 0.009685;		/* Maximal conductance for the leak component of INSCC (nS/pF)  */
	double PnsK	= 1.3;		/* Permeability of K+ */
	double PnsNa = 0.9;		/* Permeability of Na+ */
	double PnsCa = 0.89;	/* Permeability of Ca2+ */
	double PnsCs = 1.0;		/* Permeability of Cs+ */
	double gnsCa = 0.5;		/* Conductance ratio of Ca2+ in INSCC */
	double gnsNa = 1.0;		/* Conductance ratio of Na+ in INSCC */
	double gnsK = 1.19;		/* Conductance ratio of K+ in INSCC */
	double gnsCs = 1.6;		/* Conductance ratio of Cs+ in INSCC */

	double vFRT(double v) {
		return (v*frdy) / (R*temp);
	}

	double fmg() {
		return 0.108043 + 0.903902 / (1.0 + pow((mgo / 0.281007), 1.29834)); /* Magnesium inhibition of INSCC */
	}

	double enscc(CellParameters &parameters) {
		return (R*temp / frdy)*log((PnsK*ko + PnsNa*nao + 4 * PnsCa*(1.0 / (1.0 + exp(vFRT(parameters.getVoltage()))))*cao) / (PnsK*parameters.ki + PnsNa*parameters.nai + 4 * PnsCa*(1.0 / (1.0 + exp(vFRT(parameters.getVoltage()))))*parameters.cai));/* reversal potentials (mV) */
	}
	double gs_x(double x) {
		return (0.03 / (1 + (150.0 / (x + 1e-8))*(150.0 / (x + 1e-8)))) / 0.0123;	/* Extracellular ionic concentrations dependency on INSCC conductance */
	}
	double gs_ca(double x) {
		return (0.03 / (1 + (150.0 / (x + 1e-8))*(150.0 / (x + 1e-8)))) / 0.000525; /* Extracellular calcium concentrations dependency on INSCC conductance */
	}

	

	


public:
	double insca(CellParameters &parameters) {
		return fmg()*(gs_ca(cao)*gnsCa)*gns*(parameters.getVoltage() - (enscc(parameters)));/* Ca component of INSCC */
	} // This method is public to calculate the Ca2+ flux
	
	double insna(CellParameters &parameters) {
		return fmg()*(gs_x(nao)*gnsNa)*gns*(parameters.getVoltage() - (enscc(parameters))); /* Na component of INSCC */
	}

	double insk(CellParameters &parameters) {
		return fmg()*(gs_x(ko)*gnsK)*gns*(parameters.getVoltage() - (enscc(parameters))); /* K component of INSCC */
	}

	double il(CellParameters &parameters) {
		return fmg()*(gl)*(parameters.getVoltage() - (enscc(parameters))); /* leak component of INSCC */
	}

	//returns the current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return insca(parameters) + insna(parameters) + insk(parameters) + il(parameters);
	}

	// This function has a void implementation. Used to preserve the architecture of the program
	void updateGates(CellParameters &parameters) {}
};