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

#include "Cell.h"

int Cell::existingCells = 0;

void Cell::updateAllGates() {
	for (vector<Current *>::iterator it = currents.begin(); it != currents.end(); ++it) {
		(*it)->updateGates(parameters);
	}
}

void Cell::updateMyosinState(double Cai)
{
	

	double K1 = pow(Cai, parameters.nm) / (pow(Cai, parameters.nm) + pow(parameters.CaMLCK, parameters.nm));
	double K6 = K1;
	double dM_dt = -K1*parameters.M_ + parameters.K2*parameters.MP_ + parameters.K7*parameters.AM_;
	double dMp_dt = parameters.K4*parameters.AMP_ + K1*parameters.M_ - (parameters.K2 + parameters.K3)*parameters.MP_;
	double dAmp_dt = parameters.K3*parameters.MP_ + K6*parameters.AM_ - (parameters.K4 + parameters.K5)*parameters.AMP_;
	double dAM_dt = parameters.K5*parameters.AMP_ - (parameters.K7 + K6)*parameters.AM_;

	parameters.M_ += dt * dM_dt;
	parameters.MP_ += dt * dMp_dt;
	parameters.AMP_ += dt * dAmp_dt;
	parameters.AM_ += dt * dAM_dt;

}

Cell::Cell(const Cell & oldCell):cellId(oldCell.cellId)  

{
	currents = { &ina, &iCaL, &iCaT, &ib, &iK1, &iK2, &bKa, &bKab,
		&iKa, &ih, &iC1, &iNSCC, &iNaK, &iNaCaX, &iPMCA };

	this->parameters.id = oldCell.parameters.id;
}

double Cell::getTotalCurrent()
{
	
	double tot = 0;
	for (vector<Current *>::iterator it = currents.begin(); it != currents.end(); ++it) {
		tot += (*it)->getCurrent(parameters);
		
	}
	return tot;
}


double Cell::getCai()
{
	return parameters.getCai();
}

void Cell::iterate(long double diffusedV, double iStim)
{
	double i = getTotalCurrent();
	double icaTotal = getCaFlux();
	double inaTotal = getNaFlux();
	double ikTotal = getKFlux();
	double iclTotal = getClFlux();

	updateAllGates();
	updateMyosinState(parameters.cai);
	parameters.updateIonicConcentrations(icaTotal, inaTotal, ikTotal, iclTotal);

	parameters.setVoltage(parameters.getVoltage() + dt*(diffusedV - (iStim + i)/Cm));

}

string Cell::getIonCurrents() {

	string tmp;
	tmp += "Ina = " + to_string(getNaFlux()) + '\t';
	tmp += "ICa = " + to_string(getCaFlux()) + '\t';
	tmp += "IK = " + to_string(getKFlux()) + '\t';
	tmp += "ICl = " + to_string(getClFlux()) + '\t';
	return tmp;
}

void Cell::setModulations(double CaL_mod, double k_mod, double kCa_mod, double Na_mod, double CaT_mod, double MLCK_mod, double cl_mod)
{
	iCaL.setModulation(CaL_mod);
	iCaT.setModulation(CaT_mod);

	iKa.setModulation(k_mod);
	iK1.setModulation(k_mod);
	iK2.setModulation(k_mod);

	bKa.setModulation(kCa_mod);
	bKab.setModulation(kCa_mod);

	ina.setModulation(Na_mod);

	iC1.setModulation(cl_mod);

	parameters.CaMLCK *= MLCK_mod;
}

double Cell::getCaFlux()
{
	return 	(((AV*Cm*buff) / (zca*frdy))*(iCaL.getCurrent(parameters) + iCaT.getCurrent(parameters) + iNSCC.insca(parameters))
		-  iNaCaX.jnaca(parameters) + iPMCA.getCurrent(parameters));
}

// If we prefer to have the intracellular Na+, K+, and Cl- concentrations updated, uncomment the code in the following 3 methods.
double Cell::getNaFlux() {
	return 0;
	/*double a = ina.getCurrent(parameters);
	double b = iNSCC.insna(parameters);
	double c = isoc.getIsocna(parameters);
	double d = 3.0 * iNaK.getCurrent(parameters);
	double e = 3.0* 2.0* iNaCaX.getCurrent(parameters);
	double f = 0.35*ih.getCurrent(parameters);
	double g = 0.5 * iNaKCl.getCurrent(parameters);

	return ina.getCurrent(parameters) + iNSCC.insna(parameters) + isoc.getIsocna(parameters)
		+ 3.0 * iNaK.getCurrent(parameters) + 3.0* 2.0* iNaCaX.getCurrent(parameters) + 0.35*ih.getCurrent(parameters)
		- 0.5 * iNaKCl.getCurrent(parameters);*/

}

double Cell::getKFlux() {

	return 0;

	/*double a = iK1.getCurrent(parameters);
	double b = iK2.getCurrent(parameters);
	double c = iKa.getCurrent(parameters);
	double d = iKCa.getCurrent(parameters);
	double e = iNSCC.insk(parameters);
	double f = iKleak.getCurrent(parameters);
	double g = iNaKCl.getCurrent(parameters);
	double h = iNaK.getCurrent(parameters);



	return iK1.getCurrent(parameters) + iK2.getCurrent(parameters) + iKa.getCurrent(parameters) + iKCa.getCurrent(parameters)
		- 2.0*iNSCC.insk(parameters) + iKleak.getCurrent(parameters) + iNaKCl.getCurrent(parameters)
		- 2.0*iNaK.getCurrent(parameters);*/
}

double Cell::getClFlux() {
	return 0;
	/*return iCl.getCurrent(parameters) + iNaKCl.getCurrent(parameters);*/
}

// Get and set methods for the internal parameters passed to the cells' container.
double Cell::getVoltage()
{
	return parameters.getVoltage();
}

void Cell::setVoltage(double v)
{
	this->parameters.setVoltage(v);
}



CellParameters Cell::getParameters()
{
	return parameters;
}



double Cell::getAM()
{
	return parameters.AM_;
}

double Cell::getAMP()
{
	return parameters.AMP_;
}

double Cell::getLa()
{
	return parameters.la_;
}

double Cell::getLx()
{
	return parameters.lx_;
}

double Cell::getLs()
{
	return parameters.ls_;
}

double Cell::getLc()
{
	return parameters.lc_;
}

void Cell::setAM(double AM)
{
	parameters.AM_ = AM;
}

void Cell::setAMP(double AMP)
{
	parameters.AMP_ = AMP;
}

void Cell::setLa(double la)
{
	parameters.la_ = la;
}

void Cell::setLx(double lx)
{
	parameters.lx_ = lx;
}

void Cell::setLs(double ls)
{
	parameters.ls_ = ls;
}

void Cell::setLc(double lc)
{
	parameters.lc_ = lc;
}
