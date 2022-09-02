#ifndef CLUSTERS_C
#define CLUSTERS_C

#include "Cluster.h"

Cluster::Cluster() 
{
	hitsInThisCluster.clear();
	arborId = -1;
}

Cluster::~Cluster() {}


void Cluster::SetArborId(int Id)
{
	this->arborId = Id;
}

void Cluster::SetClusterId(int id)
{
	clusterId = id;
}

void Cluster::AddHits(Cell* cellToAdd)
{
	this->hitsInThisCluster.push_back(cellToAdd);
}

void Cluster::AddHits(Cluster* clusterToAdd)
{
	int size = clusterToAdd->GetSize();
	for (int i = 0; i < size; ++i)
	{
		Cell* cell = clusterToAdd->GetCell(i);
		this->AddHits(cell);
	}
}

void Cluster::SetPosition(double dx = 0., double dy = 0., double deltaZ = 0.)
{
	int size = this->GetSize();

	for (int j = 0; j < 3; ++j) meanPosition[j] = 0.;

	for (int i = 0; i < size; ++i)	
	{
		for (int j = 0; j < 3; ++j)
		{
			meanPosition[j] += hitsInThisCluster[i]->GetPosition()[j];
		}
	}

	for (int j = 0; j < 3; ++j) meanPosition[j] *= 1./(double)size;

	//options are for projection
	meanPosition[0] += dx*deltaZ;
	meanPosition[1] += dy*deltaZ;
	meanPosition[2] += deltaZ;
}


void Cluster::SetDirection()
{
	int thr1 = 1;	//seuil 0.1pC 	//tester rigoureusement quel valeurs sont à privilégier
	int thr2 = 2;	//seuil 5pC
	int thr3 = 4;	//seuil 15pC
	int n;

	this->SetPosition();
	TPrincipal* principal = new TPrincipal(3, "ND");
	principal->Clear();

	for (int i = 0; i < this->GetSize(); ++i) 
	{	
		if (this->hitsInThisCluster[i]->GetThreshold() == 3)
		{
			n = thr3;
		}else if (this->hitsInThisCluster[i]->GetThreshold() == 2)
		{
			n = thr2;
		}else
		{
			n = thr1;
		}

		double position[3];

		for (int j = 0; j < 3; ++j)
		{
			position[j] = this->hitsInThisCluster[i]->GetPosition()[j];
		}

		for (int i = 0; i < n; ++i)
		{
			principal->AddRow(position);
		}
		
	}

	principal->MakePrincipals();

	double x[3];
	double p[3] = {principal->GetEigenVectors()[0][0][0], principal->GetEigenVectors()[0][1][0], principal->GetEigenVectors()[0][2][0]};
	double x0[3];
	double p0[3] = {0, 0, 0};

	principal->P2X(p, x, 1);
	principal->P2X(p0, x0, 1);

	for (int i = 0; i < 3; ++i)	x[i] += -x0[i];


		// normalisation for z = 1
	if (abs(x[2]) > 0.1)
	{
		for (int i = 0; i < 3; ++i)	this->direction[i] = x[i]/x[2];	
	}else
	{
		for (int i = 0; i < 3; ++i)	this->direction[i] = 0.;	
	}

	delete principal;
		
}


void Cluster::SetTiming(string opt = "mean")
{
	double t, tTest;
	int size = this->GetSize();

	if (opt == "mean")
	{
		t = 0;
		for (int i = 0; i < size; ++i)
		{
			t += hitsInThisCluster[i]->GetTiming();
		}
		t *= 1./(double)size;
	}

	if (opt == "max")
	{
		t = hitsInThisCluster[0]->GetTiming();
		for (int i = 1; i < size; ++i)
		{
			tTest = hitsInThisCluster[i]->GetTiming();
			if (tTest > t)
			{
				t = tTest;
			}
		}
	}

	if (opt == "min")
	{
		t = hitsInThisCluster[0]->GetTiming();
		for (int i = 1; i < size; ++i)
		{
			tTest = hitsInThisCluster[i]->GetTiming();
			if (tTest < t)
			{
				t = tTest;
			}
		}
	}

	if (opt != "mean" && opt != "max" && opt != "min")
	{
		cout << "ERROR in Cluster::GetTiming, option is not \"max\", \"min\" or \"mean\" !" << endl;
		t = 0.;
	}
	timing = t;
}


void Cluster::AddConnectionBackward(int clusterId)
{
	connectionBackward.push_back(clusterId);
}

void Cluster::AddConnectionForward(int clusterId)
{
	connectionForward.push_back(clusterId);
}


void Cluster::ClearConnectionBackward()
{
	connectionBackward.clear();
}

void Cluster::ClearConnectionForward()
{
	connectionForward.clear();
}

void Cluster::SetSelectedBackwardConnection(int i)
{
	selectedBackwardConnection = i;
}














int Cluster::GetSize()
{
	return this->hitsInThisCluster.size();
}

Cell* Cluster::GetCell(int Id)
{
	return this->hitsInThisCluster[Id];
}

int Cluster::GetClusterId()
{
	return this->clusterId;
}

int Cluster::GetArborId()
{
	return this->arborId;
}


double* Cluster::GetPosition(bool refreshPosition = true)
{
	if (refreshPosition == true) this->SetPosition();
	return this->meanPosition;
}

double* Cluster::GetDirection()
{
	this->SetDirection();
	return this->direction;
}

int Cluster::GetMaxLayer()
{
	int memo = hitsInThisCluster[0]->GetId()[2];
	for (int i = 1; i < hitsInThisCluster.size(); ++i)
	{
		if (hitsInThisCluster[i]->GetId()[2] > memo)
		{
			memo = hitsInThisCluster[i]->GetId()[2];
		}
	}
	return memo;
}

int Cluster::GetMinLayer()
{
	int memo = hitsInThisCluster[0]->GetId()[2];
	for (int i = 1; i < hitsInThisCluster.size(); ++i)
	{
		if (hitsInThisCluster[i]->GetId()[2] < memo)
		{
			memo = hitsInThisCluster[i]->GetId()[2];
		}
	}
	return memo;
}


int Cluster::GetSelectedBackwardConnection()
{
	return selectedBackwardConnection;
}

double Cluster::GetTiming()
{
	return timing;
}

std::vector<int> Cluster::GetConnectionBackward()
{
	return connectionBackward;
}

std::vector<int> Cluster::GetConnectionForward()
{
	return connectionForward;
}





double Cluster::Distance(Cluster* otherCluster, double xCoef = 1., double yCoef = 1., double zCoef = 1.)
{
	return sqrt(pow(xCoef *(this->GetPosition(false)[0] - otherCluster->GetPosition(false)[0]), 2.)
		  	+ pow(yCoef *(this->GetPosition(false)[1] - otherCluster->GetPosition(false)[1]), 2.)
			+ pow(zCoef *(this->GetPosition(false)[2] - otherCluster->GetPosition(false)[2]), 2.));
}

double Cluster::Distance(Cell* cell, double xCoef = 1., double yCoef = 1., double zCoef = 1.)
{
	return sqrt(pow(xCoef *(this->GetPosition(false)[0] - cell->GetPosition()[0]), 2.)
		  	+ pow(yCoef *(this->GetPosition(false)[1] - cell->GetPosition()[1]), 2.)
			+ pow(zCoef *(this->GetPosition(false)[2] - cell->GetPosition()[2]), 2.));
}


void Cluster::Clustering(Cell* cellToClusterise, Detector* detector)
{
	cellToClusterise->SetClusterised();
	this->AddHits(cellToClusterise);
	this->SearchNeighbors(cellToClusterise, detector);
}

void Cluster::SearchNeighbors(Cell* centralCell, Detector* detector)
{
	int Id[3] = {centralCell->GetId()[0], centralCell->GetId()[1], centralCell->GetId()[2]};
	int xmax = detector->GetSize(0);
	int ymax = detector->GetSize(1);
	Cell* cellToTest;
	double dt;
	double dtMax = 0.2; // 0.02 ns

	for (int i = -1; i <=1; ++i)
	{for (int j = -1; j <=1; ++j)
	{
		if (Id[0]+i > -1 && Id[0]+i <xmax && Id[1]+j > -1 && Id[1]+j <ymax)
		{
			cellToTest = detector->GetCell(Id[0]+i, Id[1]+j, Id[2]);
			dt = abs(centralCell->GetTiming() - cellToTest->GetTiming());
			if (cellToTest->GetHited() == true && cellToTest->GetClusterised() == false && dt < dtMax)
			{
				this->Clustering(cellToTest, detector);
			}
		}
	}}
}







#endif