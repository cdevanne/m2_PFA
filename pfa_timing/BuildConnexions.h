#include "includes.h"
#include "Tools.C"
#include "ToolsMath.C"

#include "Arbors.C"



std::vector<Cluster*> SelectClustersForConnections(Detector* detector);
/*
	Take all hits that are not clusterised yet
*/


void CreateFirstConnections(std::vector<Cluster*> clusters, double range, double timingBymillimeterMax, double timingBymillimeterMin);
/*
	create connections between cluster if they are in range and if they are in good timing correspondance
*/

void CleanConnection(std::vector<Cluster*> clusters, double maxSpace);
/*
	Clear connection to let only 1 connection backward
*/

void CreateSecondaryConnection(std::vector<Cluster*>clusters, double range, double maxTiming, double minTiming);
/*
	Using the existing connection, search to create new connection in the same direction
*/

std::vector<Arbor*> LinkCellsWithConnections(std::vector<Cluster*> clusters);
/*
	return arbors with all clusters connected
*/

std::vector<Arbor*> MakeBranches(Detector* detector);
/*
	Use all previous functions to create arbors that start from "shower"
*/