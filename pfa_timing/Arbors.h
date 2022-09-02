#ifndef ARBORS_H
#define ARBORS_H

#include "includes.h"

#include "Cluster.C"


class Arbor
{
	//groupe de cluster relier entres eux (dans la pratique les objet clusters et arbors sont presque les même je pourrais modifier ca pour ecrire plus simplmeent)
public:
	Arbor();
	~Arbor();


// Setter

	void AddCluster(Cluster* cluster);
	/*
		add a cluster to clusterInThisArbor and add every cell of the cluster to hitsInThisArbor
	*/

	void SetArborId(int id);
	/*
		Set ArborId
	*/

	void SetPosition(double dx = 0., double dy = 0., double dz = 0.);
	/*
		Set the Position of the arbor, coef a for a projection along the z axis
	*/

	void SetDirection();
	/*
		Set the direction of ther arbor using TPrincipal
	*/

	void SetCharged(bool charged);
	/*
		set Charged
	*/

	void SetSeed(bool so = true);
	/*
		Set seed
	*/

	void SetIsSaved(bool value);
	/*
		Set isSaved
	*/

	void SetTiming(string opt = "mean");
	/*
		SetTiming 
		"mean"	 set the mean timing of all hits
		"min"	 set the minimum timing in all hits
		"max" 	 set the maximum timing in all hits
	*/

// Getter

	int GetArborId();
	/*
		return ArborId
	*/

	int GetMinLayer();
	/*
		return the layer of the lower cell
	*/

	int GetMaxLayer();
	/*
		return the max layer of the cell
	*/

	int GetThreshold(int thr);
	/*
		return the number of hits with threshold = thr in this arbors
	*/

	int GetSize();
	/*
		return the number of hits in this arbors
	*/

	int GetNumberOfClusters();
	/*
		return the number of cluster in this arbors
	*/

	int GetMainParticleId();
	/*
		return the main particle in this arbors
	*/


	Cell* GetCell(int i);
	/*
		return hitsInThisArbor[i];
	*/

	Cluster* GetCluster(int i);
	/*
		return clusterInThisArbor[i];
	*/


	double* GetPosition(bool refreshPosition = true);
	/*
		return the mean position of the arbor
	*/

	double* GetDirection();
	/*
		return the direction
	*/

	double GetTiming();
	/*
		return the timing of the arbors (as it's meaning timing, it's a bit useless)
	*/

	double GetEnergy();
	/*
		return the enrgy of, the arbors
	*/

	double GetSeedTime();
	/*
		return timing of the seed cell
	*/

	bool IsCharged();
	/*
		return isCharged
	*/

	bool GetSeed();
	/*
		return seed
	*/

	bool GetIsSaved();
	/*
		Return isSaved
	*/


// Methods

	double Distance(Arbor* otherCluster, double xCoef = 1., double yCoef = 1., double zCoef = 1.); 
	/*
		return the distance between the arbor and a clusters, default coef are cartesian coordinate system
	*/
											
	double Distance(Cell* cell,  double xCoef = 1., double yCoef = 1., double zCoef = 1.);
	/*
		return the distance between the arbor and a cell, default coef are cartesian coordinate system
	*/
	
private:
	double meanPosition[3];
	double direction[3];
	double time;

	std::vector<Cell*> hitsInThisArbor;
	std::vector<Cluster*> clusterInThisArbor;
	
	int arborId;
	int arborIdtoFusion;


	bool isSaved;
	bool chargedPrimaryTrace;	
	bool seed;
};




#endif