#include "Objects.h"



double dtCoefHere = 0.02;	// big don't take time in consideretion


//*********************************************************
///////////////////////////////////////////////////////////
//*********************************************************

Cell::Cell() {}

Cell::Cell(int xID, int yID, int zID)
{
	posID[0] = xID;
	posID[1] = yID;
	posID[2] = zID;

	hited = false;
	clusterised = false;
}

Cell::Cell(int xID, int yID, int zID, double x, double y, double z, int particle_ID, int thr, double timing)
{
	posID[0] = xID;
	posID[1] = yID;
	posID[2] = zID;

	pos[0] = x;
	pos[1] = y;
	pos[2] = z;

	hited = true;
	clusterised = false;
	particleID = particle_ID;
	threshold = thr;

	TRandom *R = new TRandom(); 
	double mean = timing;
	double sigma = 0.05;	//0.05ns

	double gaussTime = R->Gaus(mean, sigma);

	time = timing;
}

Cell::~Cell() {}

//*********************************************************

void Cell::SetHited(bool isHited = true)
{
	hited = isHited;
}

void Cell::SetClusterised(bool isClusterised = true)
{
	clusterised = isClusterised;
}

void Cell::SetArborId(int id)
{
	arborId = id;
}

//*********************************************************

bool Cell::GetHited()
{
	return hited;
}


bool Cell::GetClusterised()
{
	return clusterised;
}

int* Cell::GetID()
{
	return posID;
}

double* Cell::GetPosition()
{
	return pos;
}

int Cell::GetParticleID()
{
	return particleID;
}

int Cell::GetThreshold()
{
	return threshold;
}

double Cell::GetTime()
{
	return time;
}

int Cell::GetArborId()
{
	return arborId;
}

//*********************************************************
///////////////////////////////////////////////////////////
//*********************************************************

Map::Map(int xmax, int ymax, int zmax)
{
	detectorSize[0] = xmax;
	detectorSize[1] = ymax;
	detectorSize[2] = zmax;
	
	cellsHited.clear();
}

Map::~Map() {}

//*********************************************************

void Map::AddCell(Cell* newCell, int xID, int yID, int zID, bool hited = false)
{
	cells[xID][yID][zID] = newCell;
	if (hited == true) cellsHited.push_back(newCell);
}

//*********************************************************

Cell* Map::GetCell(int xID, int yID, int zID)
{
	return cells[xID][yID][zID];
}

std::vector<Cell*> Map::GetHits()
{
	return cellsHited;
}

int Map::GetSize(int i)
{
	return detectorSize[i];
}

//*********************************************************
///////////////////////////////////////////////////////////
//*********************************************************

Cluster::Cluster() 
{
	hitsInThisCluster.clear();
	arborID = -1;
	branchID = -1;
}

Cluster::~Cluster() {}

void Cluster::SetArborID(int ID)
{
	this->arborID = ID;
}

void Cluster::SetBranchID(int ID)
{
	this->branchID = ID;
}

void Cluster::AddHits(Cell* cellToAdd)
{
	this->hitsInThisCluster.push_back(cellToAdd);
}

void Cluster::AddHits(Cluster* clusterToAdd)
{
	for (int i = 0; i < clusterToAdd->GetSize(); ++i)
	{
		this->hitsInThisCluster.push_back(clusterToAdd->GetCell(i));
	}
}




double Cluster::Distance(Cluster* otherCluster, double xCoef = 1., double yCoef = 1., double zCoef = 1.)
{
	return sqrt(xCoef * pow(this->GetPosition(false)[0] - otherCluster->GetPosition(false)[0], 2.)
		  	+ yCoef * pow(this->GetPosition(false)[1] - otherCluster->GetPosition(false)[1], 2.)
			+ zCoef * pow(this->GetPosition(false)[2] - otherCluster->GetPosition(false)[2], 2.));
}

double Cluster::DistanceXY(Cluster* otherCluster)
{
	return sqrt(pow(this->GetPosition(false)[0] - otherCluster->GetPosition(false)[0], 2.)
		  	+ pow(this->GetPosition(false)[1] - otherCluster->GetPosition(false)[1], 2.));
}

double Cluster::Range(bool forLayer = true)
{

	double coef;
	int min = hitsInThisCluster[0]->GetID()[2];
	int max = hitsInThisCluster[0]->GetID()[2];
	int const s = 48;	//number of layer


	for (int i = 0; i < this->hitsInThisCluster.size(); i++)
	{
		int layer = hitsInThisCluster[i]->GetID()[2];
		if(min > layer) min = layer;
		if(min < layer) max = layer;
	}

	int AllLayer[s];

	for (int i = 0; i < s; ++i) AllLayer[i] = 0;

	for (int i = 0; i < this->hitsInThisCluster.size(); ++i)
	{
		AllLayer[hitsInThisCluster[i]->GetID()[2]] += 1;
	}

	max = AllLayer[0];

	for (int i = 0; i < s; ++i)
	{
		if (max < AllLayer[i])
		{
			max = AllLayer[i];
		}
	}

	if (forLayer == true)
	{
		coef = (1.+(double)this->hitsInThisCluster.back()->GetID()[2]/48.);
	}else coef = 1.;

	return 11.+sqrt(max)*11.*coef;

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

void Cluster::SetClusterId(int id)
{
	clusterId = id;
}

void Cluster::SetDirection()
{
	int thr1 = 1;	//seuil 0.1pC 	//tester rigoureusement quel valeurs sont à privilégier
	int thr2 = 4;	//seuil 5pC
	int thr3 = 20;	//seuil 15pC
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

void Cluster::SetTime(string opt = "mean")
{
	double t, tTest;
	int size = this->GetSize();

	if (opt == "mean")
	{
		t = 0;
		for (int i = 0; i < size; ++i)
		{
			t += hitsInThisCluster[i]->GetTime();
		}
		t *= 1./(double)size;
	}

	if (opt == "max")
	{
		t = hitsInThisCluster[0]->GetTime();
		for (int i = 1; i < size; ++i)
		{
			tTest = hitsInThisCluster[i]->GetTime();
			if (tTest > t)
			{
				t = tTest;
			}
		}
	}

	if (opt == "min")
	{
		t = hitsInThisCluster[0]->GetTime();
		for (int i = 1; i < size; ++i)
		{
			tTest = hitsInThisCluster[i]->GetTime();
			if (tTest < t)
			{
				t = tTest;
			}
		}
	}

	if (opt != "mean" && opt != "max" && opt != "min")
	{
		cout << "ERROR in Cluster::GetTime, option is not \"max\", \"min\" or \"mean\" !" << endl;
		t = 0.;
	}
	time = t;
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




void Cluster::Clustering(Cell* cellToClusterise, Map* map)
{
	cellToClusterise->SetClusterised();
	this->AddHits(cellToClusterise);
	this->SearchNeighbors(cellToClusterise, map);
}

void Cluster::SearchNeighbors(Cell* centralCell, Map* map)
{
	int ID[3] = {centralCell->GetID()[0], centralCell->GetID()[1], centralCell->GetID()[2]};
	int xmax = map->GetSize(0);
	int ymax = map->GetSize(1);
	Cell* cellToTest;
	double dt;
	double dtMax = dtCoefHere; // 0.02 ns

	for (int i = -1; i <=1; ++i)
	{for (int j = -1; j <=1; ++j)
	{
		if (ID[0]+i > -1 && ID[0]+i <xmax && ID[1]+j > -1 && ID[1]+j <ymax)
		{
			cellToTest = map->GetCell(ID[0]+i, ID[1]+j, ID[2]);
			dt = abs(centralCell->GetTime() - cellToTest->GetTime());
			if (cellToTest->GetHited() == true && cellToTest->GetClusterised() == false && dt < dtMax)
			{
				this->Clustering(cellToTest, map);
			}
		}
	}}
}


//*********************************************************


int Cluster::GetSize()
{
	return this->hitsInThisCluster.size();
}

Cell* Cluster::GetCell(int ID)
{
	return this->hitsInThisCluster[ID];
}

int Cluster::GetClusterId()
{
	return this->clusterId;
}

int Cluster::GetArborID()
{
	return this->arborID;
}

int Cluster::GetBranchID()
{
	return this->branchID;
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
	int memo = hitsInThisCluster[0]->GetID()[2];
	for (int i = 1; i < hitsInThisCluster.size(); ++i)
	{
		if (hitsInThisCluster[i]->GetID()[2] > memo)
		{
			memo = hitsInThisCluster[i]->GetID()[2];
		}
	}
	return memo;
}

int Cluster::GetMinLayer()
{
	int memo = hitsInThisCluster[0]->GetID()[2];
	for (int i = 1; i < hitsInThisCluster.size(); ++i)
	{
		if (hitsInThisCluster[i]->GetID()[2] < memo)
		{
			memo = hitsInThisCluster[i]->GetID()[2];
		}
	}
	return memo;
}

int Cluster::GetSelectedBackwardConnection()
{
	return selectedBackwardConnection;
}

double Cluster::GetTime()
{
	return time;
}

std::vector<int> Cluster::GetConnectionBackward()
{
	return connectionBackward;
}

std::vector<int> Cluster::GetConnectionForward()
{
	return connectionForward;
}



//*********************************************************
///////////////////////////////////////////////////////////
//*********************************************************

Arbor::Arbor() {};

Arbor::~Arbor() {}


void Arbor::AddCluster(Cluster* cluster)
{
	clusterInThisArbor.push_back(cluster);
	for (int i = 0; i < cluster->GetSize(); ++i)
	{
		hitsInThisArbor.push_back(cluster->GetCell(i));
	}
}


double Arbor::Distance(Arbor* otherCluster, double xCoef = 1., double yCoef = 1., double zCoef = 1.)
{
	return sqrt(xCoef*pow(this->GetPosition(false)[0] - otherCluster->GetPosition(false)[0], 2.)
		  	+ yCoef*pow(this->GetPosition(false)[1] - otherCluster->GetPosition(false)[1], 2.)
			+ zCoef*pow(this->GetPosition(false)[2] - otherCluster->GetPosition(false)[2], 2.));
}

double Arbor::DistanceXY(Arbor* otherCluster)
{
	return sqrt(pow(this->GetPosition(false)[0] - otherCluster->GetPosition(false)[0], 2.)
		  	+ pow(this->GetPosition(false)[1] - otherCluster->GetPosition(false)[1], 2.));
}

double Arbor::Distance(Cell* cell,  double xCoef = 1., double yCoef = 1., double zCoef = 1.)
{
	return sqrt(xCoef * pow(this->GetPosition(false)[0] - cell->GetPosition()[0], 2.)
		  	+ yCoef * pow(this->GetPosition(false)[1] - cell->GetPosition()[1], 2.)
			+ zCoef * pow(this->GetPosition(false)[2] - cell->GetPosition()[2], 2.));
}



void Arbor::SetIsSaved(bool value)
{
	isSaved = value;
}


void Arbor::SetTime(string opt = "mean")
{
	double t, tTest;
	int size = this->GetSize();

	if (opt == "mean")
	{
		t = 0;
		for (int i = 0; i < size; ++i)
		{
			t += hitsInThisArbor[i]->GetTime();
		}
		t *= 1./(double)size;
	}

	if (opt == "max")
	{
		t = hitsInThisArbor[0]->GetTime();
		for (int i = 1; i < size; ++i)
		{
			tTest = hitsInThisArbor[i]->GetTime();
			if (tTest > t)
			{
				t = tTest;
			}
		}
	}

	if (opt == "min")
	{
		t = hitsInThisArbor[0]->GetTime();
		for (int i = 1; i < size; ++i)
		{
			tTest = hitsInThisArbor[i]->GetTime();
			if (tTest < t)
			{
				t = tTest;
			}
		}
	}

	if (opt != "mean" && opt != "max" && opt != "min")
	{
		cout << "ERROR in Cluster::GetTime, option is not \"max\", \"min\" or \"mean\" !" << endl;
		t = 0.;
	}
	time = t;
}

void Arbor::SetDirection()
{
	int thr1 = 1;	//seuil 0.1pC 	//tester quel valeurs sont à privilégier
	int thr2 = 4;	//seuil 5pC
	int thr3 = 40;	//seuil 15pC
	int n;

	this->SetPosition();
	TPrincipal* principal = new TPrincipal(3, "ND");
	principal->Clear();

	for (int i = 0; i < this->GetSize(); ++i) 
	{	
		if (this->hitsInThisArbor[i]->GetThreshold() == 3)
		{
			n = thr3;
		}else if (this->hitsInThisArbor[i]->GetThreshold() == 2)
		{
			n = thr2;
		}else
		{
			n = thr1;
		}

		double position[3];

		for (int j = 0; j < 3; ++j)
		{
			position[j] = this->hitsInThisArbor[i]->GetPosition()[j];
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
	if (abs(x[2]) > 0.01)
	{
		for (int i = 0; i < 3; ++i)	this->direction[i] = x[i]/x[2];	
	}else
	{
		for (int i = 0; i < 3; ++i)	this->direction[i] = 0.;	
	}

	delete principal;
		
}


void Arbor::SetArborId(int id)
{
	arborID = id;
}

void Arbor::SetCharged(bool charged)
{
	chargedPrimaryTrace = charged;
}


void Arbor::SetPosition(double dx = 0., double dy = 0., double dz = 0.)
{
	int size = this->GetSize();

	for (int j = 0; j < 3; ++j) meanPosition[j] = 0.;

	for (int i = 0; i < size; ++i)	
	{
		for (int j = 0; j < 3; ++j)
		{
			meanPosition[j] += hitsInThisArbor[i]->GetPosition()[j];
		}
	}

	for (int j = 0; j < 3; ++j) meanPosition[j] *= 1./(double)size;

	//options are for projection
	meanPosition[0] += dx*dz;
	meanPosition[1] += dy*dz;
	meanPosition[2] += dz;
}


///////////////////////////


double* Arbor::GetPosition(bool refreshPosition = true)
{
	if (refreshPosition == true) this->SetPosition();
	return this->meanPosition;
}

bool Arbor::GetIsSaved()
{
	return isSaved;
}


double Arbor::GetTime()
{
	return time;
}

int Arbor::GetSize()
{
	return hitsInThisArbor.size();
}

void Arbor::SetSeed(bool so = true)
{
	seed = so;
}

bool Arbor::GetSeed()
{
	return seed;
}

int Arbor::GetNumberOfClusters()
{
	return clusterInThisArbor.size();
}

Cell* Arbor::GetCell(int i)
{
	return hitsInThisArbor[i];
}

Cluster* Arbor::GetCluster(int i)
{
	return clusterInThisArbor[i];
}

int Arbor::GetMinLayer()
{
	int memo = hitsInThisArbor[0]->GetID()[2];
	for (int i = 1; i < hitsInThisArbor.size(); ++i)
	{
		int id = hitsInThisArbor[i]->GetID()[2];
		if (id < memo)
		{
			memo = id;
		}
	}
	return memo;
}

int Arbor::GetMaxLayer()
{
	int memo = hitsInThisArbor[0]->GetID()[2];
	for (int i = 1; i < hitsInThisArbor.size(); ++i)
	{
		if (hitsInThisArbor[i]->GetID()[2] > memo)
		{
			memo = hitsInThisArbor[i]->GetID()[2];
		}
	}
	return memo;
}


double* Arbor::GetDirection()
{
	return direction;
}

int Arbor::GetArborId()
{
	return arborID;
}

bool Arbor::IsCharged()
{
	return chargedPrimaryTrace;
}



int Arbor::GetMainParticleID()
{
	double pID = 0;
	double undefined = 0;
	for (int i = 0; i < this->GetSize(); ++i)
	{
		if (hitsInThisArbor[i]->GetParticleID() != 0) pID += (double)hitsInThisArbor[i]->GetParticleID();
		else ++undefined;
	}
	pID *= 1./((double)this->GetSize() - undefined);
	if (pID > 1.5 ) return 2;
		else  return 1;
}


int Arbor::GetThreshold(int thr)
{
	int threshold = 0;

	for (int i = 0; i < this->GetSize(); ++i)
	{
		int thrHit = hitsInThisArbor[i]->GetThreshold();
		if (thrHit == thr)
		{
			++threshold;
		}
	}
	return threshold;
}


double Arbor::GetEnergy()
{
	double energy = 0;
	int size = this->GetSize();

	//Coef
	double alpha[3] = {0.0387388, 3.38241e-05, -1.58769e-08};
	double beta[3] = {0.0757203, -0.42138e-05, -1.50741e-10};
	double gamma[3] = {0.107314, 1.28111e-05, 3.23606e-08};

	for (int i = 0; i < size; ++i)
	{
		int thr = this->hitsInThisArbor[i]->GetThreshold();

		if (thr == 1)
		{
			energy += (alpha[0] + alpha[1]*(double)size + alpha[2]*pow((double)size, 2));
		}else if (thr == 2)
		{
			energy += (beta[0] + beta[1]*(double)size + beta[2]*pow((double)size, 2));
		}else if (thr == 3)
		{
			energy += (gamma[0] + gamma[1]*(double)size + gamma[2]*pow((double)size, 2));
		}else
		{
			cerr << "Error: Arbor::EnergyReconstruction(), threshold must be 1, 2 or 3" << endl;
		}
	}

	return energy;
}