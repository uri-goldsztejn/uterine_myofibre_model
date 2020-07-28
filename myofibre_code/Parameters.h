/* "A myofibre model for the study of uterine excitation-contraction dynamics"

Copyright(c) 2020 Uri Goldsztejn <uri.goldsztejn@wustl.edu>

This program is free software : you can redistribute it and / or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include "Constants.h"

// This class manages the cell's parameters. These parameters are used to calculate the ionic currents at each iteration.
// The parameters are updated in each iteration.

class CellParameters {
private:

	
	
	double radius = 7e-4; //cm Testrow et al.
	double surfaceArea = 2.0*PI*radius*length;
	double areaCrossSection = PI*radius*radius;
	double volume = PI*radius*radius*length; // Cell's volume [cm^3]. The volume is assumed constant	
	double force = 0;
	//constants to update ionic concentrations
	double Pna = 0.35; //permeability of sodium ions
	double Fbca = 20; //buffer fractional strength factor
	double Scm = 0.3; //total concentration of calmodulin sites for ca2+
	double Bf = 0.2; //total concentration of other buffer sites for ca2+
	double Kd = 2.6e-4; //dissociation constant
	double kdb = 5.3981e-4; //dissociation constant in buffer

	

public:
	double v = -53.90915441282156;		// Membrane potential (mV)
	double cai = 0.0001161881607214449;		// Intracellular Ca2+ (mM)
	double ki = 140.0;	/* Intracellular K+ (mM) */
	double nai = 4.0;		/* Intracellular Na+ (mM) */
	double cli = 46.0;	/* Intracellular Cl- (mM) */

	// Myosin states
	double M_ = 0.97, MP_ = 0.01, AMP_ = 0.01, AM_ = 0.01;

	// Lengths of each cellular compartment, used to calculate the force developed.
	// The initial values were determined so that the system is in equilibrium.
	double la_ = 47.3542, ls_ = 41.365, lx_ = 31.28, lc_ = 120;

	//constants for the mechanical model
	const double MlckMax = 0.84; // Maximum fraction of MLCK
	const double Kmlck = 721.749e-6; // Half-activation of MLCK [mM]
	const double Kcamlck = 1.080e-3; // Half-activation of MLCK set by [Ca2+]i [mM]
	const double Kcacm = 1.78e-7;
	const double pm = 1.0; // Hill coefficient for MLCK
	

	//Coefficients from Yochum et al for myosin cycling.
	const double nm = 4.7135; // Hill coefficient for K1
	const double K2 = 0.1399; //Myosin dephosphorylation rate constant [ms^-1]
	const double K3 = 14.4496; //Cross-bridge formation rate constant
	const double K4 = 3.6124; // Cross-bridge detachment rate constant [ms^-1]
	const double K5 = 0.1399; //Attached myosin dephosphorylation rate constant [ms^-1]
	const double K7 = 0.1340; //Latch state detachment rate constant [ms^-1]
	double CaMLCK = 0.0004640758; //this is not const because it can be changed at the beginning of the simulation by the user

	const double cms = (1.7057e8)/2.0; 
	double id = 0;
	double getVoltage() {
		return v;
	}

	void setVoltage(double voltage) {
		v = voltage;
	}

	double getCai() {
		return cai;
	}

	double getSurfaceArea() {
		return surfaceArea; /* 2.0*PI*radius*length;*/
	}

	double getAreaCrossSection() {
		return PI*radius*radius;
	}

	double getVolume() {
		return volume;
	}
	double getLength() {
		return lc_;
	}
	void setLength(double l) {
		lc_ = l;
	}
	double getForce() {
		return force;
	}
	void setForce(double Ft) {
		force = Ft;
	}


	// In this implementation only the intracellular calcium concentration is updated. 
	// Uncomment the following lines to update all the intracellular ionic concentrations.
	void updateIonicConcentrations(double Ica, double  Ina, double Ik, double Icl) {
		cai += -dt * Ica;

		/*cai -= dt * (area*Cm / (volume*zca*frdy)*Ica*Fbca*1.0 / (1.0 + Scm*Kd / pow((Kd + cai), 2) + Bf*kdb / pow(kdb + cai, 2)));

		nai -= dt * (area*Cm) / (volume*zna*frdy)*Ina;

		ki -=dt * (area*Cm) / (volume*zk*frdy)*Ik;
		
		cli -= dt * (area*Cm) / (volume*zcl*frdy)*Icl;*/
	}

	// Calculate the Nernst  potentials.
	double ek() {
		return (R*temp / frdy)*log(ko / ki);
	}

	double eca() {
		return (R*temp / frdy)*log(cao / cai);
	}

	double ena() {
		return (R*temp / frdy)*log(nao / nai);
	}
	double ecl() {
		return (R*temp / frdy)*log(cli/clo);
	}



};