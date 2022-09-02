#include "includes.h"
#include "Tools.C"
#include "ToolsMath.C"

#include "Arbors.C"



std::vector<Cluster*>  MakeClustersByLayer(Detector* detector);
/*
	Put in the same cluster all neighbors hits in the same layer (next to another)
*/

std::vector<Cluster*> SelectClustersForPrimaryTrace(std::vector<Cluster*> clustersByLayer);
/*
	return all the cluster with less than 5 hits in
*/

void SetDefaultClusterId(std::vector<Cluster*> potentialUsefullClusters);
/*
	Set every cluster as own clusterId
*/

void SearchClusterOfTheSameTrace(std::vector<Cluster*> potentialUsefullClusters);
/*
	Search a cluster above or below every cluster to create connection, Set clusterId of the clusterBelow
*/

void SetSeedClusterId(std::vector<Cluster*> potentialUsefullClusters);
/*
	Set for all related cluster the same clusterId
*/

std::vector<Cluster*> FusionClusters(std::vector<Cluster*> potentialUsefullClusters);
/*
	Fusion all cluster with ther same clusterId in a single cluster
*/

void DefineUsedCells(std::vector<Cluster*> primaryTrace, Detector* detector);
/*
	Set clustertised = true to all cluster fusionned, and set clustrised = false to all others
*/

std::vector<Cluster*> BuildPrimaryTrace(Detector* detector);
/*
	Use all previous function to create "trace" which is a line
*/