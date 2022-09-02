#ifndef CELLS_C
#define CELLS_C

#include "Cells.h"

Cell::Cell() {}

Cell::Cell(int xId, int yId, int zId)
{
	positionId[0] = xId;
	positionId[1] = yId;
	positionId[2] = zId;

	hited = false;
	clusterised = false;
}

Cell::Cell(int xId, int yId, int zId, double x, double y, double z, int particleId, int thr, double time)
{
	positionId[0] = xId;
	positionId[1] = yId;
	positionId[2] = zId;

	position[0] = x;
	position[1] = y;
	position[2] = z;

	hited = true;
	clusterised = false;
	particle = particleId;
	threshold = thr;

	TRandom *R = new TRandom(); 
	double mean = time;
	double sigma = .05;	// +- 100ps

	double gaussTime = R->Gaus(mean, sigma);

	timing = gaussTime;
	timing = time;
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

int* Cell::GetId()
{
	return positionId;
}

double* Cell::GetPosition()
{
	return position;
}

int Cell::GetParticleId()
{
	return particle;
}

int Cell::GetThreshold()
{
	return threshold;
}

double Cell::GetTiming()
{
	return timing;
}

int Cell::GetArborId()
{
	return arborId;
}

double Cell::Distance(Cell* cell, double xCoef = 1., double yCoef = 1., double zCoef = 1.)
{
	return sqrt(
				pow(xCoef * (this->GetPosition()[0] - cell->GetPosition()[0]) ,2)
				+ pow(yCoef * (this->GetPosition()[1] - cell->GetPosition()[1]) ,2)
				+ pow(zCoef * (this->GetPosition()[2] - cell->GetPosition()[2]) ,2)
				);
}


#endif
