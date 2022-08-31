#ifndef HELPER_FOR_PROGRAMMING_H
#define HELPER_FOR_PROGRAMMING_H

#include "includes.h"
#include "Tools.C"
#include "Arbors.C"


/*
	Function used only for devloppement
*/



void optimalData(std::vector<Arbor*> arbors);
/*
	Set all correct ArborId
*/

void FeedSecondParticle(std::vector<Arbor*> arbors, double range, double timing);
/*
	Put every arbors left in the second bigger arborId
*/


#endif