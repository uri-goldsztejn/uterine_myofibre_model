/* "A myofiber model for the study of uterine excitation-contraction dynamics"

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
//#include "stdafx.h" //uncomment this line to run in Microsoft's Visual Studio

#include "StringOfCellsContainer.h"

size_t StringOfCellsContainer::numberOfCells_;
vector<Cell> StringOfCellsContainer::cells_;

// Upon initialization, initialize all the cells and add them to the container.
StringOfCellsContainer::StringOfCellsContainer(Stimulator stim, size_t numberOfCells, vector<double> modulations, double resistivity) : stim_(stim){
	numberOfCells_ = numberOfCells;
	resistivityMyo = resistivity;
	conductivityGJ = 1000.0 * areaCrossSection / resistivityMyo; // mS*cm

	if (modulations.size() < 4)
	{
		cout << "missing modulation values" << endl;
		return; 
	}

	for (size_t i = 0; i < numberOfCells_; i++)
	{
		cells_.push_back(Cell());
	}
	for (size_t i = 0; i < numberOfCells_; i++)
	{
		cells_[i].setModulations(modulations[0], modulations[1], modulations[2], modulations[3], modulations[4], modulations[5], modulations[6]);
	}
};

// This method iterates over the container, makes the calculations that involve parameters from multiple cells in the myofiber and then iterates over each cell.
void StringOfCellsContainer::iterateV2() {
	long double previousVoltage, currentVoltage, nextVoltage, previousPreviousVoltage;

	// The following lines calculate the voltage diffusion across cells based on the cable equation.
	// The depolarization of the first cell is initially caused by the stimulator and then propagates to the next cells in the myofiber.
	vector<Cell>::iterator it = cells_.begin();
	double currentStimulus = stim_.getCurrent();
	previousVoltage = it->getVoltage();
	it->iterate(0, currentStimulus); // no diffused voltage in the first cell


	for (it++; it != (cells_.end() -1); ++it) {

		currentVoltage = (it)->getVoltage();
		nextVoltage = (it + 1)->getVoltage();
		
		long double diffusion_num = conductivityGJ*(nextVoltage - 2.0 * currentVoltage + previousVoltage);
		long double diffusion_denom = chi*it->getLc()*it->getLc()*C*pow(10,-8); //length changes
		
		if (it == cells_.end() - 2) {					//one cell before the last one
			previousPreviousVoltage = previousVoltage;	//save two places befroe to do backward derivative in the last cell
		}
		previousVoltage = currentVoltage;
		it->iterate(diffusion_num / diffusion_denom);

	}
	currentVoltage = (it)->getVoltage();
	
	// Voltage diffusion
	long double diffusion_num = conductivityGJ*(currentVoltage - 2.0 * previousVoltage + previousPreviousVoltage);
	long double diffusion_denom = chi*it->getLc()*it->getLc()*C*pow(10, -8); //length changes

	//Iterate over each cell
	it->iterate(diffusion_num / diffusion_denom);

	//solve lengths - after all ionic currents and voltages were updated
	
	//comment out the next line to maintain constant lengths nad make simulations faster (use for debugging only). Not calling the following function makes the simulations faster by ignoring the mechanical component of the model.
	// The results presented in the paper use the following function.
	solveMechanics(); // this function changes the lengths
}


// Returns a string consisting of the voltages of every cell, separated by a tab
string StringOfCellsContainer::getVoltages() {
	
	string output;
	for (vector<Cell>::iterator it = cells_.begin(); it != cells_.end(); ++it) {
		output += to_string(it->getVoltage()) + "\t";
	}
	output += "\n";
	return output;
}

// Returns the intracellular calcium concentrations.
string StringOfCellsContainer::getCaConcentrations()
{
	string output;

	for (vector<Cell>::iterator it = cells_.begin(); it != cells_.end(); ++it) {
		output += to_string(it->getCai()) + "\t";
	}
	output += "\n";
	return output;
}

// Returns the cellular total lenghts.
string StringOfCellsContainer::getLengths()
{
	string output;

	for (vector<Cell>::iterator it = cells_.begin(); it != cells_.end(); ++it) {
		output += to_string(it->getParameters().getLength()) + "\t";
	}
	output += "\n";
	return output;
}

// Returns the total currents.
string StringOfCellsContainer::getStringIonCurrents()
{
	string output;

	for (vector<Cell>::iterator it = cells_.begin(); it != cells_.end(); ++it) {
		output += it->getIonCurrents() + "\t";
	}
	output += "\n";
	return output;
}


int StringOfCellsContainer::equationSystem(const gsl_vector * x, void * params, gsl_vector * f)
{
	double kx1 = ((struct rparams *) params)->kx1;
	double kx2 = ((struct rparams *) params)->kx2;
	double beta = ((struct rparams *) params)->beta;
	double lopt = ((struct rparams *) params)->lopt;
	double famp = ((struct rparams *) params)->famp;
	double vx = ((struct rparams *) params)->vx;
	double fam = ((struct rparams *) params)->fam;
	double us = ((struct rparams *) params)->us;
	double ks = ((struct rparams *) params)->ks;
	double as = ((struct rparams *) params)->as;
	double ls0 = ((struct rparams *) params)->ls0;
	double kp = ((struct rparams *) params)->kp;
	double ap = ((struct rparams *) params)->ap;
	double l0 = ((struct rparams *) params)->l0;
	double lc0 = ((struct rparams *) params)->lc0;

	vector<double> l; // values are stored as la (for all the cells), then ls and then lx.
	for (size_t i = 0; i < numberOfCells_; i++)
	{
		l.push_back(gsl_vector_get(x, i));

	}
	for (size_t i = numberOfCells_; i < 2 * numberOfCells_; i++)
	{
		l.push_back(gsl_vector_get(x, i));

	}
	for (size_t i = 2 * numberOfCells_; i < 3 * numberOfCells_; i++)
	{
		l.push_back(gsl_vector_get(x, i));

	}
	

	int aIdx = 0, sIdx = numberOfCells_, xIdx = 2 * numberOfCells_;

	vector<double> y10; //the n equations #15 in paper
	for (size_t i = 0; i < numberOfCells_; i++)
	{
		

		y10.push_back((kx1*cells_[i].getAMP() + kx2*cells_[i].getAM())*l[xIdx + i] -
			famp*cells_[i].getAMP()*(vx + (l[aIdx + i] - cells_[i].getLa()) / dt) - fam*cells_[i].getAM()*(l[aIdx + i] - cells_[i].getLa()) / dt
		);
	}

	vector<double> y11; //the n equations #16 in the paper
	for (size_t i = 0; i < numberOfCells_; i++)
	{
		y11.push_back((kx1*cells_[i].getAMP() + kx2*cells_[i].getAM())*l[xIdx + i] * exp(-beta*pow((l[aIdx + i] - lopt) / lopt, 2))
			- us*((l[sIdx + i] - cells_[i].getLs()) / dt) - ks*(exp(as*(l[sIdx + i] - ls0) / ls0) - 1.0)
		);
	}

	vector<double> y12; //the n equations #17 in the paper
	for (size_t i = 0; i < numberOfCells_ - 1; i++)
	{
		y12.push_back((kx1*cells_[i].getAMP() + kx2*cells_[i].getAM())*l[xIdx + i] * exp(-beta*pow((l[aIdx + i] - lopt) / lopt, 2))
			+ kp*(exp(ap*(l[aIdx + i] + l[sIdx + i] + l[xIdx + i] - l0) / l0) - 1.0)

			- (kx1*cells_[i + 1].getAMP() + kx2*cells_[i + 1].getAM())*l[xIdx + i + 1] * exp(-beta*pow((l[aIdx + i + 1] - lopt) / lopt, 2))
			- kp*(exp(ap*(l[aIdx + i + 1] + l[sIdx + i + 1] + l[xIdx + i + 1] - l0) / l0) - 1.0)

		);
	}

	// Equation #18 in the paper.
	double y13 = 0;
	for (size_t i = 0; i < numberOfCells_; i++)
	{
		y13 += l[aIdx + i] + l[sIdx + i] + l[xIdx + i];
	}

	y13 -= numberOfCells_ * lc0;

	for (size_t i = 0; i < numberOfCells_; i++) // in f vector the equations are in the following order: first the #15 i =1...n, then #16, then #17 and finally #18
	{
		gsl_vector_set(f, i, y10[i]);
		gsl_vector_set(f, numberOfCells_ + i, y11[i]);
		if (i < numberOfCells_ - 1) {
			gsl_vector_set(f, 2 * numberOfCells_ + i, y12[i]);
		}
	}

	gsl_vector_set(f, 3 * numberOfCells_ - 1, y13); // this value is written over existing value


	return GSL_SUCCESS;
}

// Uncomment the lines in the following method to print the current state of the solver.
int StringOfCellsContainer::print_state(size_t iter, gsl_multiroot_fsolver * s)
{
	/*printf("iter = %3u x = % .3f % .3f % .3f"
		"f(x) = % .3e % .3e % .3e\n",
		iter,
		gsl_vector_get(s->x, 0),
		gsl_vector_get(s->x, numberOfCells_),
		gsl_vector_get(s->x, 2 * numberOfCells_),
		gsl_vector_get(s->f, 0),
		gsl_vector_get(s->f, numberOfCells_),
		gsl_vector_get(s->f, 2 * numberOfCells_));*/

	return 0;
}

// Initialize the PDE solver from GSL.
void StringOfCellsContainer::setSolver()
{
	iter = 0;
	f = { &equationSystem,3 * numberOfCells_, &p };
	x = gsl_vector_alloc(3 * numberOfCells_);

	for (size_t i = 0; i < numberOfCells_; i++)
	{
		
		gsl_vector_set(x, i, 47.3542); //initial conditions
		gsl_vector_set(x, numberOfCells_ + i, 41.365);
		gsl_vector_set(x, 2 * numberOfCells_ + i, 31.28);
		
	
	}

	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc(T, 3 * numberOfCells_);
	gsl_multiroot_fsolver_set(s, &f, x);
	
	print_state(iter, s);
}

// Use the solver to solve the mechanical equations in each iteration.
void StringOfCellsContainer::solveMechanics()
{
	setSolver();

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);

		print_state(iter, s);

		if (status)   /* check if solver is stuck*/
			break;

		status = gsl_multiroot_test_residual(s->f, 1e-7);
	} while (status == GSL_CONTINUE && iter < 1000);
	//printf("status = %s\n", gsl_strerror(status));

	//save results into cells
	for (size_t i = 0; i < numberOfCells_; i++)
	{
		cells_[i].setLa(gsl_vector_get(s->x, i));
		cells_[i].setLs(gsl_vector_get(s->x, numberOfCells_ + i));
		cells_[i].setLx(gsl_vector_get(s->x, 2 * numberOfCells_ + i));

		cells_[i].setLc(cells_[i].getLa() + cells_[i].getLs() + cells_[i].getLx());
	}

	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
}

// Get the total force developed by a given cell.
double StringOfCellsContainer::getFtotalAtCelli(int i)
{
	return (p.kx1*cells_[i].getAMP() + p.kx2*cells_[i].getAM())*cells_[i].getLx()* exp(-p.beta*pow((cells_[i].getLa() - p.lopt) / p.lopt, 2))
		+ p.kp*(exp(p.ap*(cells_[i].getLc() - p.l0) / p.l0) - 1.0);
}




// Returns the lengths of the cell in the container
string StringOfCellsContainer::getCellLengths()
{
	string lengths;
	for (size_t i = 0; i < numberOfCells_; i++)
	{
		lengths += to_string(cells_[i].getLc()) + ' ';
	}
	lengths += '\n';
	return lengths;
}

// Returns the force developed by the myofiber. Form equation #17, we know that all the cells develop the same total force.
string StringOfCellsContainer::getForces()
{
	string forces = to_string(getFtotalAtCelli(1)) + '\n';

	return forces;
	
}

// Get the fraction of myosin in the force producing states as a string to print in the output files.
string StringOfCellsContainer::getMyosinForceFraction()
{
	string fraction;
	for (size_t i = 0; i < numberOfCells_; i++)
	{
		fraction += to_string(cells_[i].getAM()) + ' ' + to_string(cells_[i].getAMP()) + ' ';
	}
	fraction += '\n';
	return fraction;
}

// Returns the cells' compartments lengths, the total force produced by the string, and the myosin concentrations in the force producing states, as a string to print in the output files.

string StringOfCellsContainer::getLengthsForce()
{
	string lengths;
	for (size_t i = 0; i < numberOfCells_; i++)
	{
		lengths += to_string(cells_[i].getLa()) + ' ' + to_string(cells_[i].getLs()) + ' ' + to_string(cells_[i].getLx()) + ' ' + to_string(cells_[i].getLc()) +' ';// +'\t' + to_string(getFtotalAtCelli(i)) + '\n';
	}
	lengths += '\t' + to_string(getFtotalAtCelli(1));
	lengths += '\t' + to_string(cells_[1].getAM());
	lengths += '\t' + to_string(cells_[1].getAMP()) + '\n';

	return lengths;
}

// Returns the myosin states.
string StringOfCellsContainer::getMyosinStates()
{
	string states;
	for (size_t i = 0; i < numberOfCells_; i++)
	{
		states += to_string(cells_[i].getAM()) + ' ' + to_string(cells_[i].getAMP());
	}

	states += '\n';

	return states;
}

// Returns the current produced by the stimulator.
double StringOfCellsContainer::readCurrent() {
	return stim_.readCurrent();
}