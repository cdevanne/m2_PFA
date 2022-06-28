#include <TH2.h>
#include <TStyle.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <TPrincipal.h>
#include <TRandom.h>





class Cell
{
	//cellule du detecteur avec position, timing, hit etc ...
public:
	Cell();
	Cell(int xID, int yID, int zID);
	Cell(int xID, int yID, int zID, double x, double y, double z, int particle_ID, int thr, double timing);
	~Cell();

	void SetHited(bool isHited);
	void SetClusterised(bool isClusterised);
	void SetArborId(int id);

	int* GetID();
	double* GetPosition();
	bool GetHited();
	bool GetClusterised();
	int GetParticleID();
	int GetThreshold();
	double GetTime();
	int GetArborId();

private:
	int posID[3];
	double pos[3];
	int particleID;
	int threshold;
	double time;
	int arborId;

	bool hited;
	bool clusterised;
};


class Map
{
	//l'enssemble des cellules
public:
	Map(int xmax, int ymax, int zmax);
	~Map();

	void AddCell(Cell* newCell,int xID, int yID, int zID,  bool hited);

	int GetSize(int i);
	Cell* GetCell(int xID, int yID, int zID);
	std::vector<Cell*> GetHits(); 

private:
	int detectorSize[3];
	Cell* cells[96][96][48];
	std::vector<Cell*> cellsHited;
};

class Cluster
{
	//groupe de cellule voisinje sur une meme couche
public:
	Cluster();
	~Cluster();

	void AddHits(Cell* cellToAdd);
	void AddHits(Cluster* clusterToAdd);
	void AddConnectionBackward(int clusterId);
	void AddConnectionForward(int clusterId);
	void ClearConnectionBackward();
	void ClearConnectionForward();
	void Clustering(Cell* cellToClusterise, Map* map);	//pour creer les groupe de hits voisin (cluster)
	void SearchNeighbors(Cell* centralCell, Map* map);	//pour creer les groupe de hits voisin (cluster)

	double Range(bool forLayer = true);
	double Distance(Cluster* otherCluster, double xCoef = 1., double yCoef = 1., double zCoef = 1.);
	double DistanceXY(Cluster* otherCluster);

	void SetClusterId(int id);
	void SetArborID(int ID);
	void SetBranchID(int ID);
	void SetPosition(double dx = 0., double dy = 0., double deltaZ = 0.);
	void SetDirection();
	void SetSelectedBackwardConnection(int i);
	void SetTime(string opt = "mean");
	
	int GetSize();
	Cell* GetCell(int ID);
	int GetArborID();
	int GetBranchID();
	int GetClusterId();
	int GetMinLayer();
	int GetMaxLayer();
	double* GetDirection();
	double* GetPosition(bool refreshPosition = true);
	int GetMainParticleID();
	double GetTime();
	std::vector<int> GetConnectionBackward();
	std::vector<int> GetConnectionForward();
	int GetSelectedBackwardConnection();
	




private:
	std::vector<Cell*> hitsInThisCluster;
	double time;
	double direction[3];
	double meanPosition[3];
	int clusterId;
	int arborID;
	int branchID;
	std::vector<int> connectionBackward;
	std::vector<int> connectionForward;
	int selectedBackwardConnection;
};



class Arbor
{
	//groupe de cluster relier entres eux (dans la pratique les objet clusters et arbors sont presque les même je pourrais modifier ca pour ecrire plus simplmeent)
public:
	Arbor();
	~Arbor();



	void AddCluster(Cluster* cluster);

	double Distance(Arbor* otherCluster, double xCoef = 1., double yCoef = 1., double zCoef = 1.); 	//distance carthésienne ou non
	double DistanceXY(Arbor* otherCluster);															//distance selon le plan XY
	double Distance(Cell* cell,  double xCoef = 1., double yCoef = 1., double zCoef = 1.);


	void SetArborId(int id);
	void SetPosition(double dx = 0., double dy = 0., double dz = 0.);
	void SetDirection();			//utilise la classe TPrincipal
	void SetCharged(bool charged);
	void SetSeed(bool so = true);
	void SetIsSaved(bool value);

	int GetSize();
	int GetNumberOfClusters();
	int GetMainParticleID();
	void SetTime(string opt = "mean");
	double GetTime();
	Cell* GetCell(int i);
	Cluster* GetCluster(int i);
	double* GetPosition(bool refreshPosition = true);
	double* GetDirection();
	bool IsCharged();
	int GetArborId();
	int GetMinLayer();
	int GetMaxLayer();
	int GetThreshold(int thr);
	bool GetIsSaved();
	double GetEnergy();
	bool GetSeed();
	
private:
	double meanPosition[3];
	double direction[3];
	bool chargedPrimaryTrace;
	std::vector<Cell*> hitsInThisArbor;
	std::vector<Cluster*> clusterInThisArbor;
	int arborID;
	int arborIDtoFusion;
	double time;
	bool isSaved;
	bool seed;
};