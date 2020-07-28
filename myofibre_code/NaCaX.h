#pragma once
#include "Current.h"

class NaCaX : public Current {
private:
	double Jnaca = 3.5e-6; // 1.75e-5;	/* maximal calcium flux via Na-Ca exchanger (mM/ms) */
	double Kmallo = 0.0003;	/* half saturation concentration for calcium allosteric factor  (mM) */
	double nallo = 4;		/* Hill coefficient of allosteric regulation (dimensionless) */
	double ksat = 0.27;		/* saturation factor at very negative potential */
	double xgamma = 0.35;		/* partition parameter  */
	double Kmnai = 30.0;		/* dissociation constants for intracellular Na, Ca ions (mM) */
	double Kmcai = 0.007;	/* dissociation constants for intracellular Na, Ca ions (mM) */
	double Kmnao = 87.5;		/* dissociation constants for extracellular Na, Ca ions (mM) */
	double Kmcao = 1.3;		/* dissociation constants for extracellular Na, Ca ions (mM) */


	double vFRT(double v) {
		return (v*frdy) / (R*temp);
	}
	double f1naca(double v) {
		return exp((xgamma - 1.0)*vFRT(v));
	}
	double f2naca(double v) {
		return exp(xgamma*vFRT(v));
	}
	double fallo(double cai) {
		return 1.0 / (1.0 + pow((Kmallo / cai), nallo));
	}
	double naca_Eup(double v, double cai, double nai) {
		return Jnaca*((nai*nai*nai)*cao*f2naca(v) - (nao*nao*nao)*cai*f1naca(v));
	}
	double naca_Ed1(double v) {
		return 1.0 + ksat*f1naca(v);
	}
	double naca_Ed2(double cai, double nai) {
		return Kmcao*(nai*nai*nai) + (Kmnao*Kmnao*Kmnao)*cai + (Kmnai*Kmnai*Kmnai)*cao*(1.0 + cai / Kmcai);
	}
	double naca_Ed3(double cai, double nai) {
		return cao*(nai*nai*nai) + (nao*nao*nao)*cai + (nao*nao*nao)*Kmcai*(1.0 + (nai / Kmnai)*(nai / Kmnai)*(nai / Kmnai));
	}


public:

	//returns the current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return 0.5*jnaca(parameters)*((zca*frdy) / (AV*Cm*buff));	/* Current via Na-Ca exchanger (pA/pF) */
		// the 1/2 factor is added here instead of during the total calculation
	}

	double jnaca(CellParameters &parameters) {
		return fallo(parameters.getCai())*naca_Eup(parameters.getVoltage(), parameters.getCai(), parameters.nai) 
			/ (naca_Ed1(parameters.getVoltage())*(naca_Ed2(parameters.getCai(), parameters.nai) + naca_Ed3(parameters.getCai(), parameters.nai)));/* Calcium flux via Na-Ca exchanger (mM/ms)*/
	} // This method is public to can calculate the Ca flux

	// This function has a void implementation. Used to preserve the architecture of the program
	void updateGates(CellParameters &parameters){}

};