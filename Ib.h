#pragma once
#include "Current.h"
//from Testrow et al.
class Ib : public Current {
private:
	double gb = 0.004; 			/* background current conductances (nS/pF) */
	
	
public:
	//returns the channel current obtained in the present state
	double getCurrent(CellParameters &parameters) {
		return gb*(parameters.getVoltage() - parameters.ek()); /* Ib */
	}

	// updates the ion channel gates for the following iteration
	void updateGates(CellParameters &parameters) {}
};