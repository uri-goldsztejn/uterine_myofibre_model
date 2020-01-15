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

#pragma once
#include "Constants.h"
#include "Parameters.h"
// Base class for all ionic currents.
// This abstract class is used to iterate over all currents and allow for future refactorings with better ionic current models.

class Current {
protected:
	
public:

	virtual double getCurrent(CellParameters &parameters) = 0;
	virtual void updateGates(CellParameters &parameters) = 0;

};