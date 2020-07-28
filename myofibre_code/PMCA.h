#pragma once
#include "Current.h"

class PMCA : public Current {
private:
	double Jpmca = 3.5e-7;	/* Maximal calcium flux  via PMCA (mM/ms)*/
	double Kmpmca = 0.0005;	/* half saturation concentration (mM)*/
	double npmca = 2;		/* Hill coefficient */

public:

	//returns the current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return Jpmca / (1.0 + pow((Kmpmca / parameters.getCai()), npmca)); /* Calcium flux via PMCA (mM/ms) */
	}

	// This function has a void implementation. Used to preserve the architecture of the program
	void updateGates(CellParameters &parameters) {}
};