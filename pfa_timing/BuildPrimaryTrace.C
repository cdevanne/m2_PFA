#ifndef BUILD_PRIMARY_TRACE_C
#define BUILD_PRIMARY_TRACE_C


#include "BuildPrimaryTrace.h"

std::vector<Cluster*>  MakeClustersByLayer(Detector* detector)
{
	Cell* cellTested;
	std::vector<Cluster*> clustersByLayer;
	bool clusterised, hited;

	for (int i = 0; i < detector->GetHits().size(); ++i)
	{
		cellTested = detector->GetCell(detector->GetHits()[i]->GetId()[0], detector->GetHits()[i]->GetId()[1], detector->GetHits()[i]->GetId()[2]);
		clusterised = cellTested->GetClusterised();
		hited = cellTested->GetHited();

		if (clusterised == false && hited == true)
		{
			Cluster* cluster = new Cluster();
			cluster->Clustering(cellTested, detector);

			int size = cluster->GetSize();
			int clusterId = clustersByLayer.size();
			cluster->SetClusterId(clusterId);
			cluster->SetPosition();
			clustersByLayer.push_back(cluster);				
		}
	}

	return clustersByLayer;
}


std::vector<Cluster*> SelectClustersForPrimaryTrace(std::vector<Cluster*> clustersByLayer)
{
	std::vector<Cluster*> potentialUsefullClusters; 

	for (int i = 0; i < clustersByLayer.size(); ++i)
	{
		int sizeCluster = clustersByLayer[i]->GetSize();

		if (sizeCluster < 5)
		{
			potentialUsefullClusters.push_back(clustersByLayer[i]);
		}
	}
	clustersByLayer.clear();

	return potentialUsefullClusters;
}

void SetDefaultClusterId(std::vector<Cluster*> potentialUsefullClusters)
{
	int numberOfClusters = potentialUsefullClusters.size();

	for (int i = 0; i < numberOfClusters; ++i)
	{
		potentialUsefullClusters[i]->SetClusterId(i);
	}
}


void SearchClusterOfTheSameTrace(std::vector<Cluster*> potentialUsefullClusters)
{

	int numberOfClusters = potentialUsefullClusters.size();

	for (int i = 0; i < numberOfClusters; ++i)
	{
		int layerI = potentialUsefullClusters[i]->GetMaxLayer();
		for (int j = 0; j < numberOfClusters; ++j)
		{
			if (i != j)
			{
				int layerJ = potentialUsefullClusters[j]->GetMinLayer();
				int diffLayer = layerJ - layerI;
				
				if (diffLayer < 4 && diffLayer > 0)	// diffLayer = 2 allow a hole in hits
				{
					double distance = potentialUsefullClusters[i]->Distance(potentialUsefullClusters[j], 1., 1., 0.);
					if (distance < 15.)
					{
						potentialUsefullClusters[j]->SetClusterId(i);
					}
				}
			}
		}
	}
}

void SetSeedClusterId(std::vector<Cluster*> potentialUsefullClusters)
{

	int numberOfClusters = potentialUsefullClusters.size();

	for (int i = 0; i < numberOfClusters; ++i)
	{
		int j = potentialUsefullClusters[i]->GetClusterId();


		if (j != i)
		{
			int jj =  potentialUsefullClusters[j]->GetClusterId();
			while(jj != j)
			{
				j = jj;
				jj =  potentialUsefullClusters[j]->GetClusterId();
			}
			potentialUsefullClusters[i]->SetClusterId(j);
		}
	}
}


std::vector<Cluster*> FusionClusters(std::vector<Cluster*> potentialUsefullClusters)
{
	std::vector<Cluster*> primaryTrace; 

	int numberOfClusters = potentialUsefullClusters.size();

	SetDefaultClusterId(potentialUsefullClusters);
	SearchClusterOfTheSameTrace(potentialUsefullClusters);
	SetSeedClusterId(potentialUsefullClusters);
	

	for (int i = 0; i < numberOfClusters; ++i)
	{
		if (potentialUsefullClusters[i]->GetClusterId() == i)
		{
			int id = primaryTrace.size();

			primaryTrace.push_back(potentialUsefullClusters[i]);

			int lenght = 1;

			for (int j = 0; j < numberOfClusters; ++j)
			{
				if (i == j) continue;

				if (potentialUsefullClusters[j]->GetClusterId() == i)
				{
					primaryTrace[id]->AddHits(potentialUsefullClusters[j]);
					++lenght;
				}
			}

			if (lenght < 3)
			{
				primaryTrace.erase (primaryTrace.begin() + id);
			}
		}
	}
	potentialUsefullClusters.clear();

	return primaryTrace;
}

void DefineUsedCells(std::vector<Cluster*> primaryTrace, Detector* detector)
{
	std::vector<Cell*> hits = detector->GetHits(); 

	for (int i = 0; i < hits.size(); ++i)
	{
		hits[i]->SetClusterised(false);
	}

	for (int i = 0; i < primaryTrace.size(); ++i)
	{
		for (int j = 0; j < primaryTrace[i]->GetSize(); ++j)
		{
			Cell* cell = primaryTrace[i]->GetCell(j);
			cell->SetClusterised(true);
		}
	}
}



std::vector<Cluster*> BuildPrimaryTrace(Detector* detector)
{
	std::vector<Cluster*> primaryTrace; 
	
	primaryTrace = MakeClustersByLayer(detector);
	primaryTrace = SelectClustersForPrimaryTrace(primaryTrace);
	primaryTrace = FusionClusters(primaryTrace);

	DefineUsedCells(primaryTrace, detector);

	return primaryTrace;
}



#endif