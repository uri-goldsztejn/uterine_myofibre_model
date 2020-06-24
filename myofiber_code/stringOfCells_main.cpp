/* "A myofiber model for the study of uterine excitation-contraction dynamics"

The original article can be found at:

If you found this study useful, please cite our work.

The ionic channel model used here is based on the model developed by Tong et al.
Tong, et al. "A computational model of the ionic currents, Ca2+ dynamics and action potentials underlying contraction of isolated uterine smooth muscle."
PloS one 6.4 (2011): e18685.

Copyright (c) 2020 Uri Goldsztejn <uri.goldsztejn@wustl.edu>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



//#include"stdafx.h" //uncomment this line to run in Microsoft's Visual Studio
#include "Stimulator.h"
#include "Cell.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <istream>
#include "StringOfCellsContainer.h"

using namespace std;

int writeMetadata(ofstream& file, double simulationTime);
bool is_number(const std::string& s);

int main()
{
	string tmp;
	// define simulation variables
	cout << "Enter simulation parameters, default values appear in parentheses" << endl;
	
	// simulation name
	string simName;
	cout << "Enter a name for the simulation:" << endl;
	cin >> simName;
	string simulationName = "./simulations/" + simName;

	// Define the intercellular resistivity
	cout << "Intercellular Resistivity (IR) [Ohm*cm] (150): ";
	cin >> tmp;
	double resistivity = 150.0;
	if (is_number(tmp)) {
		resistivity = stod(tmp);
	}

	// The modulations regulate some ionic currents. They are multiplicative factors.
	//Entering a modulation = 1.0, means that the ionic current is not modulated. a modulation of 0.5, halves the conductivity for the specified ionic channel.
	vector<double> modulations;
	vector<string> modulationMessages;
	modulationMessages.push_back("CaL modulation [~] (1): ");
	modulationMessages.push_back("K modulation [~] (1): ");
	modulationMessages.push_back("K(ca) modulation [~] (1): ");
	modulationMessages.push_back("Na modulation [~] (1): ");
	modulationMessages.push_back("CaT modulation [~] (1): ");
	modulationMessages.push_back("MLCK modulation [~] (1): ");
	modulationMessages.push_back("Cl modulation [~] (1): ");

	for (size_t i = 0; i < modulationMessages.size(); i++)
	{
		cout << modulationMessages[i];
		cin >> tmp;
		modulations.push_back(1);
		if (is_number(tmp)) {
			modulations[i] = stod(tmp);
		}
	}

	int simulationTime = 1000; //ms
	size_t numberCells = 70;
	//Stimulator options
	double stON = 70.0, stOFF = 100.0, stBCL = 5000.0, tBON = 0.0, tBOFF = 100.0, stHold = 0.0, stStep = -5.0;
	cout << "Change more parameters? [Y/N]: ";
	cin >> tmp;
	if (tmp == "Y" || tmp =="y")
	{
		// simulation time
		cout << "Enter simulation time [ms] (14000): ";
		cin >> tmp;
		if (is_number(tmp)) {
			simulationTime = stoi(tmp);
		}

		// number of cells
		cout << "Enter number of cells (70): ";
		cin >> tmp;
		if (is_number(tmp)) {
			numberCells = stoi(tmp);
		}

		
		// stimulator ON
		cout << "Enter time to turn on the stimulus [ms] (0): ";
		cin >> tmp;
		if (is_number(tmp)) {
			stON = stoi(tmp);
		}

		// stimulator OFF
		cout << "Enter time to turn off the stimulus [ms] (30): ";
		cin >> tmp;
		if (is_number(tmp)) {
			stOFF = stoi(tmp);
			tBOFF = stOFF; //Comment out this line if uncommenting the following lines
		}

		// Uncomment the following lines to add periodic stimulations to the stimulator

		/*
		// stimulator BCL
		cout << "Enter bsic cycle length, for periodic stimulations [ms] (0): ";
		cin >> tmp;
		if (is_number(tmp)) {
			stBCL = stoi(tmp);
		}

		// stimulator BON
		cout << "Enter time to turn on the stimulator during a periodic stimulation, for periodic stimulations [ms] (0): ";
		cin >> tmp;
		if (is_number(tmp)) {
			tBON = stoi(tmp);
		}

		// stimulator BOFF
		cout << "Enter time to turn off the stimulator during a periodic stimulation, for periodic stimulations [ms] (0): ";
		cin >> tmp;
		if (is_number(tmp)) {
			tBOFF = stoi(tmp);
		}

		// stimulator st_hold
		cout << "Enter hold level of the stimulator [pA/pF] (0): ";
		cin >> tmp;
		if (is_number(tmp)) {
			stHold = stoi(tmp);
		}*/

		// stimulator st_step
		cout << "Enter step level of the stimulator [pA/pF] (-5): ";
		cin >> tmp;
		if (is_number(tmp)) {
			stStep = stoi(tmp);
		}

	}

	Stimulator stimulator(stON, stOFF, stBCL, tBON, tBOFF, stHold, stStep);
	Cell cell;
	
	StringOfCellsContainer stringOfCells(stimulator, numberCells, modulations, resistivity); 
	// Create files to store data
	string voltageFilename = simulationName + "_voltage.txt";
	string calciumFilename = simulationName + "_calcium.txt";
	string metadataFilename = simulationName + "_metadata.txt";
	string lengthsFilename = simulationName + "_lengths.txt";
	string forceFilename = simulationName + "_force.txt";
	string phsosphorilationFilename = simulationName + "_phsosphorilation.txt";

	
	ofstream  voltageData, forceData, metaData, calciumData, mechanicalData, lengthsData, phosphoData;																							  

	voltageData.open(voltageFilename.c_str());
	calciumData.open(calciumFilename.c_str());
	metaData.open(metadataFilename.c_str());
	writeMetadata(metaData, simulationTime);
	lengthsData.open(lengthsFilename.c_str());
	forceData.open(forceFilename.c_str());
	phosphoData.open(phsosphorilationFilename.c_str());

	// The output of this program can be too big. To downsample the output uncomment the following 2 lines and the commented lines in the for loop below.
	//int numberOfLines = 1000;
	//int downsampleRate = round(simulationTime / dt / numberOfLines);

	metaData << "Resistivity = " << resistivity << " [Ohm*cm]" << endl;
	metaData << "CaL modulation = " << modulations[0] << endl;
	metaData << "k modulation = " << modulations[1] << endl;
	metaData << "k(Ca) modulation = " << modulations[2] << endl;
	metaData << "Na modulation = " << modulations[3] << endl;
	metaData << "CaT modulation = " << modulations[4] << endl;
	metaData << "MLCK modulation = " << modulations[5] << endl;
	metaData << "Cl modulation = " << modulations[6] << endl;
	

	metaData << "\n\n\nSimulation based on the model developed in \"A myofibre model for the study of uterine excitation-contraction dynamics\"" << endl;
	metaData << "By Uri Goldsztejn and Arye Nehorai" << endl;
	//metaData << "This simulation was performed using software version " << version << endl;
	metaData << "This software is protected by the GNU General Public License" << endl;

	// write data to output files
	for (int i = 0; i < round(simulationTime / dt); i++)
	{

		//	if (i % downsampleRate == 0) // write to files every "downsampleRate" iterations
		//	{
		voltageData << stringOfCells.getVoltages();
		calciumData << stringOfCells.getCaConcentrations();
		lengthsData << stringOfCells.getCellLengths();
		forceData << stringOfCells.getForces();
		phosphoData << stringOfCells.getMyosinForceFraction();

		//	}


		stringOfCells.iterateV2();
		// Print to screen to keep track of simulation progress
		if (i % int(1000) == 0) {
		
			cout << "iteration: " << i << " = " << i*dt << " ms of simulation" << endl;
		}

		
	}

	voltageData.close();
	calciumData.close();

	lengthsData.close();
	forceData.close();
	phosphoData.close();
	metaData.close();

	return 0;
}

//write timestamp, code version, constants, notes, etc.
int writeMetadata(ofstream& file, double simulationTime) {

	//intro
	file << "#### Uterine myofiber simulation ####\n" << endl;

	cout << "write notes below (they will be saved in the metadata file):" << endl;
	string notes;
	cin >> notes;
	file << "Notes from user:" << endl;
	file << notes << "\n" <<endl;
	

	file << "Parameters used in this simulation:" << endl;
	file << "total simulated time = " << simulationTime << " [ms]" << endl;
	file << "dt = " << dt << " [ms]" << endl;

	return 0;

}

// check if the input is a number, sanity check for input.
bool is_number(const std::string& s)
{

	std::string::const_iterator it = s.begin();
	while (it != s.end() && (isdigit(*it) || *it == '.')) ++it;
	return !s.empty() && it == s.end();
}