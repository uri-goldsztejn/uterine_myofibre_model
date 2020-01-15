#pragma once
#include"Constants.h"
#include <math.h>

// This implementation is adapted form Tong et al.

class Stimulator {
private:
	double stON_, stOFF_, stBCL_, tBON_, tBOFF_, stHold_, stStep_;
	double t = 0;
// ST_ON	1000.0		/* Time to turn on stimulus (ms) */
// ST_OFF	3000.0		/* Time to turn off stimulus (ms) */
// ST_BCL	0.0		/* Basic cycle length, require only for periodic stimulation (ms)*/
// ST_BON	0.0		/* Time to turn on a stimulus during a periodic stimulation (ms) */
// ST_BOFF	0.0		/* Time to turn off a stimulus during a periodic stimulation (ms) */
// ST_HOLD	0.0		/* Holding level of the stimulus (pA/pF) */
// ST_STEP	-0.5		/* Stepping level of the stimulus (pA/pF) */

public:
	Stimulator(double stON, double stOFF, double stBCL, double tBON, double tBOFF, double stHold, double stStep):
		stON_(stON), stOFF_(stOFF), stBCL_(stBCL), tBON_(tBON), tBOFF_(tBOFF), stHold_(stHold), stStep_(stStep){};
	
	double getCurrent();
	double readCurrent();
};