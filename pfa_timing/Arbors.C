#ifndef ARBORS_C
#define ARBORS_C

#include "Arbors.h"



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


void Arbor::SetTiming(string opt = "mean")
{
	double t, tTest;
	int size = this->GetSize();

	if (opt == "mean")
	{
		t = 0;
		for (int i = 0; i < size; ++i)
		{
			t += hitsInThisArbor[i]->GetTiming();
		}
		t *= 1./(double)size;
	}

	if (opt == "max")
	{
		t = hitsInThisArbor[0]->GetTiming();
		for (int i = 1; i < size; ++i)
		{
			tTest = hitsInThisArbor[i]->GetTiming();
			if (tTest > t)
			{
				t = tTest;
			}
		}
	}

	if (opt == "min")
	{
		t = hitsInThisArbor[0]->GetTiming();
		for (int i = 1; i < size; ++i)
		{
			tTest = hitsInThisArbor[i]->GetTiming();
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
	time = t;
}

void Arbor::SetDirection()
{
	int thr1 = 1;	//seuil 0.1pC 	//tester quel valeurs sont à privilégier
	int thr2 = 1;	//seuil 5pC
	int thr3 = 1;	//seuil 15pC
	int n;

	int thr[3] = {0};

	for (int i = 0; i < this->GetSize(); ++i) 
	{
		++thr[this->hitsInThisArbor[i]->GetThreshold() - 1];
	}

	if (thr[1] > 5)
	{
		thr2 = 4;
		thr3 = 8;
	}

	if (thr[2] > 5)
	{
		thr2 = 4;
		thr3 = 40;
	}

	this->SetPosition();
	TPrincipal* principal = new TPrincipal(3, "ND");
	principal->Clear();

	for (int i = 0; i < this->GetSize(); ++i) 
	{	
		if (this->hitsInThisArbor[i]->GetThreshold() <2.5 || this->hitsInThisArbor[i]->GetId()[2] > this->GetMinLayer() + 6) continue;
		

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
	arborId = id;
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


double Arbor::GetTiming()
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
	int memo = hitsInThisArbor[0]->GetId()[2];
	for (int i = 1; i < hitsInThisArbor.size(); ++i)
	{
		int id = hitsInThisArbor[i]->GetId()[2];
		if (id < memo)
		{
			memo = id;
		}
	}
	return memo;
}

int Arbor::GetMaxLayer()
{
	int memo = hitsInThisArbor[0]->GetId()[2];
	for (int i = 1; i < hitsInThisArbor.size(); ++i)
	{
		if (hitsInThisArbor[i]->GetId()[2] > memo)
		{
			memo = hitsInThisArbor[i]->GetId()[2];
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
	return arborId;
}

bool Arbor::IsCharged()
{
	return chargedPrimaryTrace;
}

double Arbor::GetSeedTime()
{
	int minLayer = this->GetMinLayer();
	for (int i = 0; i < this->GetNumberOfClusters(); ++i)
	{
		Cluster* cluster = this->GetCluster(i);
		if (cluster->GetMinLayer() == minLayer)
		{
			return cluster->GetTiming();
		}
	}
	return 0.;
}



int Arbor::GetMainParticleId()
{
	double pId = 0;
	double undefined = 0;
	for (int i = 0; i < this->GetSize(); ++i)
	{
		if (hitsInThisArbor[i]->GetParticleId() != 0) pId += (double)hitsInThisArbor[i]->GetParticleId();
		else ++undefined;
	}
	pId *= 1./((double)this->GetSize() - undefined);
	if (pId > 1.5 ) return 2;
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


#endif