#ifndef TOOLS_H
#define TOOLS_H

#include "includes.h"
#include "Arbors.C"
#include "HelperForProgramming.h"
#include "ToolsMath.h"


Cluster* GetSeedCluster(Arbor* a);
Cluster* GetTopCluster(Arbor* a);


//I don't know why but the next line have to bet quite

// double distanceSeed(Arbor* a, Arbor* b, double coefx = 1., double coefy = 1., double coefz = 1.);
/*
	return the number of arbor having a number of hit in the interval
*/


int* arborsMaxSizeTop3(std::vector<Arbor*> arbors);
/*
	return the arborId of the 3 bigger arbor
*/


double NumberOfArborInSize(std::vector<Arbor*> arbors, int sizeMin, int sizeMax);
/*
	return the number of arbor having a number of hit in the interval
*/


double distance(std::vector<Arbor*> arbors, int size = 40, bool max = true);
/*
	compute distance between all arbor and return the longest
*/


void IdentifieChargedParticle(std::vector<Arbor*> arbors);
/*
	search if therre is a trace in the first layers on the detector to identifie all possible charged particle
*/


int FindSeedArbor(std::vector<Arbor*> arbors, int id, std::vector<int> toDelete);
/*
	return the id of the seed cell
*/


void SetArborId(std::vector<Arbor*> arbors);
/*
	Set ArborId of every arbor with the own Id
*/

std::vector<Arbor*> DeleteDuplicateData(std::vector<Arbor*> arbors);
/*
	TEMPORARY
	in some situation,some hits are duplicate; (1 every 50 event approx)
	this is an error, I have to find where to fix it.

	but this function just search duplicate hits and delete the duplication
*/

#endif