#ifndef CLUSTERS_H
#define CLUSTERS_H

#include "includes.h"

#include "Detector.C"



class Cluster
{
	//groupe de cellule voisinje sur une meme couche
public:
	Cluster();
	~Cluster();


// Setter

	void SetClusterId(int id);
	/*
		Set the clusterId
	*/

	void SetArborId(int Id);
	/*
		Set the arborId
	*/

	void SetPosition(double dx = 0., double dy = 0., double dz = 0.);
	/*
		Set the position
	*/

	void SetDirection();
	/*
		Set the direction of the cluster using TPrincipal
	*/

	void SetTiming(string opt = "mean");
	/*
		Set the timing
	*/

	void AddHits(Cell* cellToAdd);
	/*
		add a hit in the cluster
	*/

	void AddHits(Cluster* clusterToAdd);
	/*
		add all hits from another cluster in this cluster
	*/

	void AddConnectionBackward(int clusterId);
	/*
		add a connection backward
	*/

	void AddConnectionForward(int clusterId);
	/*
		add a connection forward
	*/

	void ClearConnectionBackward();
	/*
		clear the vector connectionForward
	*/

	void ClearConnectionForward();
	/*
		clear the vector connectionBackward
	*/

	void SetSelectedBackwardConnection(int i);
	/*
		set the selected connection
	*/

// Getter
	
	int GetSize();
	/*
		Return the number of cell in the cluster
	*/

	Cell* GetCell(int Id);
	/*
		return a specific cell
	*/

	int GetArborId();
	/*
		return the arborId
	*/

	int GetClusterId();
	/*
		return the clusterId
	*/

	int GetMinLayer();
	/*
		return the layer of the lower cell
	*/

	int GetMaxLayer();
	/*
		return the layer of the top cell
	*/

	double* GetDirection();
	/*
		return the vector of the direction of the cluster
	*/

	double* GetPosition(bool refreshPosition = true);
	/*
		return the position
	*/

	int GetMainParticleId();
	/*
		return the main particle
	*/

	double GetTiming();
	/*
		return the timing
	*/

	std::vector<int> GetConnectionBackward();
	/*
		return the vector connectionForward
	*/

	std::vector<int> GetConnectionForward();
	/*
		return the vector connectionBackward
	*/

	int GetSelectedBackwardConnection();
	/*
		return the selectodConnectionBackward
	*/

//Methods

	void Clustering(Cell* cellToClusterise, Detector* detector);
	/*
		Set cluisterised at a cell,create a cluster then use SearchNeighbors on it
	*/

	void SearchNeighbors(Cell* centralCell, Detector* detector);	
	/*
		Search close hit in the same layer to add in ther same cluster
	*/

	double Distance(Cluster* otherCluster, double xCoef = 1., double yCoef = 1., double zCoef = 1.);
	/*
		return the distance between 2 clusters, default coef are cartesian coordinate system
	*/

	double Distance (Cell* cell, double xCoef = 1., double yCoef = 1., double zCoef = 1.);
	/*
		return the distance a cell and a clusters,  default coef are cartesian coordinate system
	*/
	



private:
	std::vector<Cell*> hitsInThisCluster;
	double timing;
	double direction[3];
	double meanPosition[3];
	int clusterId;
	int arborId;
	std::vector<int> connectionBackward;
	std::vector<int> connectionForward;
	int selectedBackwardConnection;	
};



#endif