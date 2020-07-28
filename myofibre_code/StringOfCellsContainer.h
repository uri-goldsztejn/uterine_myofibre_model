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
#include "Stimulator.h"
#include "Cell.h"
#include <string>
#include <iostream>
#include <gsl_vector.h>
#include <gsl_multiroots.h>


// This class represents the myofibre. It is a container for the cells. The container stores the cells and performs calculations that involve the entire array of cells.
// Only one instance of this class is used in the program (singleton).

class StringOfCellsContainer {
private:

	static size_t numberOfCells_;  // The number of cells in the myofibre
	Stimulator stim_;			   // The stimulator used to stimulate the first cell in the myofibre.
	static vector<Cell> cells_;	   // This vector contains the cells of the myofibre

	// constant mechanical parameters used to calculate the force developed
	struct rparams
	{
		double kx1 = 12.5;
		double kx2 = 8.8;
		double beta = 7.5;
		double lopt = 100.0;
		double famp = 130;
		double vx = 5;
		double fam = 85.5;
		double us = 0.01;
		double ks = 0.2;
		double as = 4.5;
		double ls0 = 30.0;
		double kp = 0.1;
		double ap = 0.1;
		double l0 = 40.0;
		double lc0 = 120;
	} p;
	double resistivityMyo;
	double conductivityGJ; // mS*cm

	// Variables for ODE solver
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;
	int status;
	size_t  iter = 0;
	gsl_multiroot_function f;
	gsl_vector *x;

public:
	// Constructor for the container
	StringOfCellsContainer(Stimulator stim, size_t numberOfCells, vector<double> modulations, double resistivity);

	// Used to update the container (and make the calculations).
	void iterateV2();

	// Read states of the cells in the container
	string getVoltages();
	string getCaConcentrations();
	double readCurrent();
	string getLengths();
	string getStringIonCurrents();

	
	// Methods used to manage the solver used for the system of mechanical equations
	static int equationSystem(const gsl_vector * x, void *params, gsl_vector * f); // This function has to be static to be passed to solver
																				   // That requires that the member variables used in this function should be static too (not elegant)
	int	print_state(size_t iter, gsl_multiroot_fsolver * s);
	void setSolver();
	void solveMechanics();

	// Get mechanical states of the cells
	double getFtotalAtCelli(int i);

	string getCellLengths();
	string getForces();
	string getMyosinForceFraction();

	string getLengthsForce();

	string getMyosinStates();
	
};
