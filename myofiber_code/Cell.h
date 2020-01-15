#pragma once

#include <iterator>
#include <vector>

#include "Ina.h"
#include "ICaL.h"
#include "ICaT.h"
#include "Ib.h"
#include "IK1.h"
#include "IK2.h"
#include "BKa.h"
#include "BKab.h"
#include "IKa.h"
#include "Ih.h"
#include "ICl.h"
#include "INSCC.h"
#include "INaK.h"
#include "NaCaX.h"
#include "PMCA.h"

#include "Parameters.h"
#include <string>
#include <ostream>
#include <iostream>
using namespace std;


//This object represents a single cell that will be part of the myofiber
class Cell
{
protected:
	int cellId;			//Serial number
	static int existingCells;
	
	CellParameters parameters;
	
	// These currents were derived in Tong et al. They can be changed for better models in the future
	Ina ina; ICaL iCaL; ICaT iCaT; Ib ib; IK1 iK1; IK2 iK2; BKa bKa; BKab bKab;
	IKa iKa; Ih ih; ICl iC1; INSCC iNSCC; INaK iNaK; NaCaX iNaCaX; PMCA iPMCA;

	vector<Current *> currents = { &ina, &iCaL, &iCaT, &ib, &iK1, &iK2, &bKa, &bKab,
		&iKa, &ih, &iC1, &iNSCC, &iNaK, &iNaCaX, &iPMCA };


public:
	
	//We add a serial number to every cell to keep track of them
	Cell(){ 
		cellId = existingCells++;
		parameters.id = cellId;
	}

	// A copy constructor used to add the cells in the container that represents the myofiber
	Cell(const Cell &oldCell);

	double getTotalCurrent();
	double getCai();
	//Update the cell's state in each iteration
	void updateAllGates();
	void updateMyosinState(double cai); 
	void iterate(long double diffusionV, double iStim = 0);

	// Obtain the different ionic fluxes
	double getCaFlux();
	double getNaFlux();
	double getKFlux();
	double getClFlux();
	double getVoltage();
	void setVoltage(double dv);
	CellParameters getParameters();

	string getIonCurrents();

	// Initializes the modulations entered by the user
	void setModulations(double CaL_mod, double k_mod, double kCa_mod, double Na_mod, double CaT_mod, double MLCK_mod, double cl_mod);


	// These set-get methods are used to manage the phosporilation states and cellular compartments' lengths
	double getAM();
	double getAMP();
	double getLa();
	double getLx();
	double getLs();
	double getLc();

	void setAM(double AM);
	void setAMP(double AMP);
	void setLa(double la);
	void setLx(double lx);
	void setLs(double ls);
	void setLc(double lc);
};