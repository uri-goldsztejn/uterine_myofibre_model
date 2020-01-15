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


#include "Stimulator.h"

// Returns the current and advances the clock
double Stimulator::getCurrent()
{
	if ((t >= stON_ && t <= stOFF_) && (fmod(t,stBCL_) >= tBON_ && fmod(t,stBCL_) <= tBOFF_))
	{
		t += dt;
		return stStep_;
	}
 	t += dt;
	return stHold_;
}

// Returns the current without advancing the clock
double Stimulator::readCurrent()
{
	if ((t >= stON_ && t <= stOFF_) && (fmod(t, stBCL_) >= tBON_ && fmod(t, stBCL_) <= tBOFF_))
	{
		
		return stStep_;
	}

	return stHold_;
}
