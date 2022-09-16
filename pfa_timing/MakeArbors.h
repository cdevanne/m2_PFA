#include "includes.h"
#include "Tools.C"
#include "ToolsMath.C"
#include "HelperForProgramming.C"


#include "Arbors.C"


void CreateConnectionsBetweenBranchesAndtrace(std::vector<Arbor*> branches, std::vector<Cluster*> trace);
/*
	use the trace crated and the arbors produce from buildconnections 
	to connect them if the end of a trace is close to the start of an arbor
	and set ArborId if they can fusion
*/

std::vector<Arbor*> MakeFirstArbors(std::vector<Arbor*> branches, std::vector<Cluster*> trace);
/*
	fusion trace and arbors after CreateConnectionsBetweenBranchesAndtrace
*/

std::vector<Arbor*> CleanArbor(std::vector<Arbor*> arbors);
/*
	Fusion all arbors which must fusion acording to their arborId
	Use at the end of almost every next function
*/

std::vector<Arbor*> CommonSeedArbor(std::vector<Arbor*> arbors, double range, double coefZ, double layerMax);
/*
	search if some arbors start from the same area and fusion them
*/

std::vector<Arbor*> BuildArborSeedMissing(std::vector<Arbor*> arbors, bool methodCloseLayers, int maxLayerForFusion, double range, double minSize, double maxSize, double timing, double timingMin);
/*
	search if two arbors came from the same area and are close enought 
	(in case of a neutral particule sending two particule in very different direction for example)
	rarely usefull
*/
std::vector<Arbor*> FusionArbor(std::vector<Arbor*> arbors, double maxDistance, int minSize, int maxSize, double timing, double timingMin);
/*
	use the trajectory of the arbor and search if it can come from another one, 
	care about timing, distance and energy
*/

std::vector<Arbor*> FusionArborFromSeed(std::vector<Arbor*> arbors, double maxDistance, int minSize, int maxSize, double timing, double timingMin);
/*
	almost the same function as FusionArbor but search from the basis of the arbor instead of the mean position
*/

std::vector<Arbor*> ClusteringSmallArbor(std::vector<Arbor*> arbors, double sizeMaxFood, double sizeMinEater, double maxRange, double timing, double timingMin);
/*
	search isolated hit and small arbors close to bigger arbors, and fusions them
*/

std::vector<Arbor*> FinalRecombination(std::vector<Arbor*> arbors, double energyMax, double range, bool protectCharged);
/*
	Last step, Use mainsly a distance parameter, try to fusion arbors and check if the energy
	reconstitution is fine
*/

std::vector<Arbor*> MakeArbors(std::vector<Arbor*> branches, std::vector<Cluster*> trace);
/*
	Use a combinaison of all the previous function to try to build particle showers
*/