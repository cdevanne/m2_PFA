#ifndef TOOLS_C
#define TOOLS_C

#include "Tools.h"






Cluster* GetSeedCluster(Arbor* a)
{
	int layerI = a->GetMinLayer();
	Cluster* c1;

	for (int i = 0; i < a->GetNumberOfClusters(); ++i)
	{
		c1 = a->GetCluster(i);
		if (c1->GetMinLayer() == layerI) break;
	}
	return c1;
}

Cluster* GetTopCluster(Arbor* a)
{
	int layerI = a->GetMaxLayer();
	Cluster* c1;

	for (int i = 0; i < a->GetNumberOfClusters(); ++i)
	{
		c1 = a->GetCluster(i);
		if (c1->GetMinLayer() == layerI) break;
	}
	return c1;
}


double distanceSeed(Arbor* a, Arbor* b, double coefx = 1., double coefy = 1., double coefz = 1.)
{
	int layerI = a->GetMinLayer();
	int layerJ = b->GetMinLayer();
	double distance;
	Cluster* c1;
	Cluster* c2;

	for (int i = 0; i < a->GetNumberOfClusters(); ++i)
	{
		c1 = a->GetCluster(i);
		if (c1->GetMinLayer() == layerI) break;
		
	}
	for (int i = 0; i < b->GetNumberOfClusters(); ++i)
	{
		c2 = b->GetCluster(i);
		if (c2->GetMinLayer() == layerJ) break;
	}

	return c1->Distance(c2, coefx, coefy, coefz);
}


double TotalEnergy(Arbor* a, Arbor* b)
{
	Arbor* fusion = new Arbor();

	for (int i = 0; i < a->GetNumberOfClusters(); ++i) fusion->AddCluster(a->GetCluster(i));
	for (int i = 0; i < b->GetNumberOfClusters(); ++i) fusion->AddCluster(b->GetCluster(i));
	
	return fusion->GetEnergy();
}


int* arborsMaxSizeTop3(std::vector<Arbor*> arbors)
{

	int max = 0;
	int max2 = max;
	int max3 = max2;
	int id1 = -1;
	int id2 = -1;
	int id3 = -1;
	
	for (int i = 0; i < arbors.size(); ++i)	
	{
		arbors[i]->SetArborId(i);
		int size = arbors[i]->GetSize();
		if (size > max)	
		{
			max3 = max2;
			id3 = id2;
			max2 = max;
			id2 = id1;
			max = size;
			id1 = i;

		}else if (size > max2)
		{
			max3 = max2;
			id3 = id2;
			max2 = size;
			id2 = i;
		}
		else if (size > max3)
		{

			max3 = size;
			id3 = i;
		}
	}

	int id[3];

	id[0] = id1;
	id[1] = id2;
	id[2] = id3;
	return id;
}


double NumberOfArborInSize(std::vector<Arbor*> arbors, int sizeMin, int sizeMax)
{
	int count = 0;
	for (int i = 0; i < arbors.size(); ++i)
	{
		if (arbors[i]->GetSize() > sizeMin && arbors[i]->GetSize() < sizeMax)
		{
			++count;
		}
	}

	return count;
}


double distance(std::vector<Arbor*> arbors, int size = 40, bool max = true)
{
	double dist = -1;
	if (max == false) dist = 1e20;

	for (int i = 0; i < arbors.size(); ++i)
	{
		if (arbors[i]->GetSize() > size)
		{
			for (int j = i + 1 ; j < arbors.size(); ++j)
			{
				if (arbors[j]->GetSize() > size)
				{
					double dist2 = distanceSeed(arbors[j], arbors[i], 1., 1., 0.05);
					if (dist2 > dist && max == true) dist = dist2;
					if (dist2 < dist && max == false) dist = dist2;
				}
				
			}
		}
	}	
	return dist;
}



void IdentifieChargedParticle(std::vector<Arbor*> arbors, int numberOfChargedParticles)
{
	for (int i = 0; i < arbors.size(); ++i)
	{
		int countHit = 0;
		int size = arbors[i]->GetSize();

		for (int j = 0; j < size; ++j)
		{
				
				Cell* cell = arbors[i]->GetCell(j);
				int layer = cell->GetId()[2];

				if (layer < 4)
				{
					++countHit;
				}

		}
	
		if (countHit >= 3)
		{
			arbors[i]->SetCharged(true);
		}else
		{
			arbors[i]->SetCharged(false);
		}
	}
	int id[numberOfChargedParticles];

}

int FindSeedArbor(std::vector<Arbor*> arbors, int id, std::vector<int> toDelete)
{
	//les fonction qui vont suivre vont attitrer des arbres parent aux arbres parents, pour faciliter les association cette fonction serta trouver l'arbre parent d'un arbre (au cas ou un arbres parents soit lui même le fils d'un arbre tiers)
	if (id < 0 || id >= arbors.size())
	{
		cerr << "ERROR: FindSeedArbor has a wrong value id = " << id << endl;
		return 0;
	}
	int arborId = id;
	int previousArborId = -1;

	std::vector<int> memory;
	std::vector<int> loopMemory;
	bool loop = false;

	memory.push_back(arborId);

	while ((arborId != previousArborId && loop == false))
	{
		previousArborId = arborId;
		arborId = arbors[previousArborId]->GetArborId();
		memory.push_back(arborId);


		for (int i = memory.size() - 1; i > 0; --i)
		{
			for (int j = i-1; j >= 0; --j)
			{
				if (memory[i] == memory[j])
				{
					int difference = i-j;
					if (j - difference > 0)
					{
						if(memory[j - difference] == memory[j])
						{
							loop = true;
							for (int k = 0; k < difference; ++k)
							{
								loopMemory.push_back(memory[j+k]);
							}
						}
					}
				}
			}			
		}
	}

	bool helper = false;
	int j = 0;
	if (loop == true)
	{
		sort(loopMemory.begin(), loopMemory.end());
		while (helper == false && j < loopMemory.size())
		{
			helper = true;
			arborId = loopMemory[j];
			for (int i = 0; i < toDelete.size(); ++i)
			{
				if (arborId == toDelete[i])
				{
					helper = false;
				}
			}
			++j;
		}

	}

	if (j == toDelete.size())
	{
		arborId = -1;
	}

	for (int i = 0; i < toDelete.size(); ++i)
	{
		if (arborId == toDelete[i])
		{
			arborId = -1;
			break;
		}
	}


	return arborId;
}


void SetArborId(std::vector<Arbor*> arbors)
{
	for (int i = 0; i < arbors.size(); ++i)
	{
		int size = arbors[i]->GetSize();
		for (int j = 0; j < size; ++j)
		{
			Cell* cell = arbors[i]->GetCell(j);
			if (size > 9)
			{
				cell->SetArborId(j);
			}else
			{
				cell->SetArborId(-1);
			}
		}
	}
}


#endif