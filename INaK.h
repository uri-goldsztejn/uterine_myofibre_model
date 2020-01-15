#pragma once
#include "Current.h"

class INaK : public Current {
private:
	double ginak = 1.7;// 3.4;		/* Maximal current for INaK (pA/pF)  */
	double nakKmko = 2.0;	/* Half-saturation concentration for ko dependency */
	double nakKmnai = 22.0;	/* Half-saturation concentration for Nai dependency */

	double vFRT(double v) {
		return (v*frdy) / (R*temp);
	}

	double fnak(double v) {
		return 1.0 / (1.0 + 0.1245*exp(-0.1*vFRT(v)) + 2.19e-3*(exp(nao / 49.71))*exp(-1.9*vFRT(v)));/* Voltage-dependency of INaK */
	}

	double knak() {
		return 1.0 / (1.0 + pow((nakKmko / ko), 1.5)); /* ko-dependency of INaK */
	}
	double nnak(double nai) {
		return 1.0 / (1.0 + pow((nakKmnai / nai), 2)); /* nai-dependency of INaK */
	}

public:
	//returns the channel current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return ginak*knak()*nnak(parameters.nai)*fnak(parameters.getVoltage()); // INaK 
	}

	// updates the ion channel gates for the following iteration
	void updateGates(CellParameters &parameters) {}
};