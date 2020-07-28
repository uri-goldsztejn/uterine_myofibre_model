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
#include <math.h>

// Constants
#define zca	2.0 		/* Valency for Ca2+ */
#define zna	1.0 		/* Valency for Na+ */
#define zk	1.0			/* Valency for K+ */
#define zcl -1			/* Valency for Cl- */
#define R 	8314.0      	/* Universal Gas Constant (J/kmol*K) */
#define frdy	96485.0   	/* Faraday's Constant (C/mol) */
#define temp	308.0   	/* Temperature (K) */
#define Cm	1.0		/* Specific membrane capacitance (uF/cm2) */
#define PI 3.141592

/* ionic concentrations
modified from Okabe et al 1999 Eur. J Pharmacology 376:101-108*/
const double ko = 6.0;		/* Extracellular K+ (mM) */
const double cao = 2.5;		/* Extracellular Ca2+ (mM) */
const double nao = 130.0;	/* Extracellular Na+ (mM) */
const double clo = 130.0;	/* Extracellular Cl- (mM) */
const double mgo = 0.5;		/* Extracellular Mg2+ (mM) */

const double buff = 0.015;	/* Proportion of free Calcium ions from calcium flux */
const double AV = 4.0;		/* Area to Volume ratio (cm^-1) */

//Calculation constants
const double dt = 0.02;		/* Time step (ms) */

// Potassium conductances
	/* Maximal conductance for total voltage-gated K currents (IK1, IK2, IKa) (nS/pF)  */
const double gkca = 0.8;		/* Maximal conductance for total calcium-dependent K currents (BKa, BKab) (nS/pF)  */


// Constants for AP propagation in the myofibre
const double radius = 7e-4; //cm Testrow et al.
const double length = 120e-4; //cm Testrow et al.
const double areaCrossSection = PI*radius*radius;
const double resistanceGJ = 4e5; //Ohm Shaw and Rudy
const double chi = 7.422e-5 / 2.65e-8; // cm^-1 from Testrow et al.
const double C = Cm*7.422e-5; // uF

//code version
const double version = 1.0;