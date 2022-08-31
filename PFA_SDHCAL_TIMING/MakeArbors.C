#ifndef MAKE_ARBORS_C
#define MAKE_ARBORS_C


#include "MakeArbors.h"

void CreateConnectionsBetweenBranchesAndtrace(std::vector<Arbor*> branches, std::vector<Cluster*> trace)
{
	for (int i = 0; i < branches.size(); ++i) branches[i]->SetArborId(-1);

	//search at the top of trace

	for (int i = 0; i < trace.size(); ++i)
	{
		int layer = trace[i]->GetMaxLayer();
		int size = trace[i]->GetSize();
		Cell* cellTrace;

		for (int j = 0; j < size; ++j)
		{
			cellTrace = trace[i]->GetCell(j);
			if (cellTrace->GetId()[2] == layer) break;
		}

		for (int j = 0; j < branches.size(); ++j)
		{
			int layerArbor = branches[j]->GetMinLayer();

			if (layer > layerArbor || layerArbor - layer > 7) continue;

			size = branches[j]->GetSize();
			Cell* cellArbor;

			for (int k = 0; k < size; ++k)
			{
				cellArbor = branches[j]->GetCell(k);
				if (cellArbor->GetId()[2] == layerArbor)
				{
					double distance = cellTrace->Distance(cellArbor, 1., 1., 0.2);

					if (distance < 40.)
					{
						branches[j]->SetArborId(i);
					}
				}
			}
		}
	}

	//search at the bot of trace
	double distanceMin = 40.;

	for (int i = 0; i < trace.size(); ++i)
	{
		int layer = trace[i]->GetMinLayer();
		int size = trace[i]->GetSize();
		Cell* cellTrace;

		for (int j = 0; j < size; ++j)
		{
			cellTrace = trace[i]->GetCell(j);
			if (cellTrace->GetId()[2] == layer) break;
		}

		for (int j = 0; j < branches.size(); ++j)
		{
			int layerArbor = branches[j]->GetMinLayer();

			if (layer < layerArbor || - layerArbor + layer > 2) continue;

			size = branches[j]->GetSize();
			Cell* cellArbor;

			for (int k = 0; k < size; ++k)
			{
				cellArbor = branches[j]->GetCell(k);
				if (cellArbor->GetId()[2] == layerArbor)
				{
					double distance = cellTrace->Distance(cellArbor, 1., 1., 0.);

					if (distance < distanceMin)
					{
						branches[j]->SetArborId(i);
						distance = distanceMin;
					}
				}
			}
		}
	}

	return branches;
}

std::vector<Arbor*> MakeFirstArbors(std::vector<Arbor*> branches, std::vector<Cluster*> trace)
{
	std::vector<Arbor*> arbors;

	for (int i = 0; i < trace.size(); ++i)
	{
			int arborId = arbors.size();
			trace[i]->SetArborId(arborId);
			Arbor* newArbor = new Arbor();
			newArbor->AddCluster(trace[i]);
			newArbor->SetArborId(arborId);
			arbors.push_back(newArbor);
	}

	for (int i = 0; i < branches.size(); ++i)
	{
		int id = branches[i]->GetArborId();
		if (id == -1)
		{
			arbors.push_back(branches[i]);
		}else
		{
			for (int j = 0; j < branches[i]->GetNumberOfClusters(); ++j)
			{
				Cluster* cluster = branches[i]->GetCluster(j);
				arbors[id]->AddCluster(cluster);
			}
		}
	}
	return arbors;
}




std::vector<Arbor*> CleanArbor(std::vector<Arbor*> arbors)
{
	//une fois les arbres parent initiaux trouvercette fonction se charge de combiner tout le monde
	std::vector<int> toDelete;

	for (int i = 0; i < arbors.size() ; ++i)
	{
		int arborId = arbors[i]->GetArborId();

		if (arborId != i)
		{
			if (arborId < 0 /*|| (arbors[arborId]->GetSeed() == true && arbors[i]->GetSeed() == true)*/)
			{
				arbors[i]->SetArborId(i);
			}else
			{
				toDelete.push_back(i);
			}
			
		}			
	}
	
	sort(toDelete.begin(), toDelete.end(), greater<int>());

	for (int i = 0; i < arbors.size(); ++i)
	{
		arbors[i]->SetIsSaved(false);
	}

	for (int i = 0; i < toDelete.size(); ++i)
	{
		int targetId = FindSeedArbor(arbors, toDelete[i], toDelete);
		if (targetId >= 0)
		{
			int size = arbors[toDelete[i]]->GetNumberOfClusters();

			for (int j = 0; j < size; ++j)
			{
				Cluster* cluster = arbors[toDelete[i]]->GetCluster(j);
				arbors[targetId]->AddCluster(cluster);
				
			}
			arbors[toDelete[i]]->SetIsSaved(true);
		}
	}

	for (int i = 0; i < toDelete.size(); ++i)
	{
		if (arbors[toDelete[i]]->GetIsSaved() == true)
		{
			arbors.erase (arbors.begin() + toDelete[i]);
		}
	}

	return arbors;
}


std::vector<Arbor*> CommonSeedArbor(std::vector<Arbor*> arbors, double range, double coefZ, double layerMax)
{
	double distanceMin = range;

	for (int i = 0; i < arbors.size(); ++i)
	{
		arbors[i]->SetArborId(i);
	}


	for (int i = 0; i < arbors.size(); ++i)
	{
		Cluster* c1 = GetSeedCluster(arbors[i]);


		for (int j = i; j < arbors.size(); ++j)
		{
			if (arbors[j]->GetSize() < 10) continue;
			Cluster* c2 = GetSeedCluster(arbors[j]);

			if (abs(c1->GetMinLayer() - c2->GetMinLayer()) > layerMax) continue;

			if (c1->GetMinLayer() < 2 && c2->GetMinLayer() < 2) 
			{
				range = 25.;
				coefZ = 0.;

			}
		

			double distanceSeed = c1->Distance(c2, 1., 1., coefZ);

			if (distanceSeed < range)
			{
				arbors[i]->SetArborId(j);
				distanceMin = distanceSeed;
			}		
		}
	}
	arbors = CleanArbor(arbors);

	return arbors;
}



std::vector<Arbor*> BuildArborSeedMissing(std::vector<Arbor*> arbors, bool methodCloseLayers, int maxLayerForFusion, double range, double minSize, double maxSize, double timing, double timingMin)
{
	//deux fonction utiles pour les arbres neutre notement, depend de la valeur du booleen "methodCloseLayers":
	
	//		1)si un arbres a sa racine sur une couche, cherche dans un petit perimetre si un autres arbres a également sa racine ici pour les associer
	
	//		2)a partie de l'orientation de deux arbres, regarde si ils sembles pouvoir venir d'une même région en amont

	std::vector<Cluster*> potentialPrimaryTrace; 

	for (int i = 0; i < arbors.size(); ++i)
	{
		arbors[i]->SetDirection();
		arbors[i]->SetArborId(i);
		arbors[i]->SetTiming();
	}

	for (int i = arbors.size() - 1; i > 0 ; --i)
	{

		int layerI = arbors[i]->GetMinLayer();

		//if (arbors[i]->GetMaxLayer() - layerI < 4) continue;
		if (arbors[i]->GetSize() > maxSize || arbors[i]->GetSize() < minSize) continue;

		double* directionI = arbors[i]->GetDirection();
		double* positionI = arbors[i]->GetPosition();

		double t1 = 0; 
		int helper = 0;
		for (int k = 0; k < arbors[i]->GetNumberOfClusters(); ++k)
		{

			Cluster* clusterI = arbors[i]->GetCluster(k);
			if (clusterI->GetMinLayer() < layerI + 2)
			{
				t1 += clusterI->GetTiming();
				++helper;
			}
		}

		t1 *= 1/(double)helper;
		

		Cluster* clusterI;
		for (int k = 0; k < arbors[i]->GetNumberOfClusters(); ++k)
		{
			clusterI = arbors[i]->GetCluster(k);
			if (clusterI->GetMinLayer() == layerI)
			{
					
				
				for (int j = arbors.size() - 1; j >= 0; --j)
				{
					if (i == j) continue;

					int layerJ = arbors[j]->GetMinLayer();
					int diffLayer = layerI - layerJ; 
					double distanceMax = range;
					double* directionJ = arbors[j]->GetDirection();
					double* positionJ = arbors[j]->GetPosition();

					if (diffLayer < 0) continue;
					

					double t2; 

					helper = 0;
					Cluster* clusterJ;
					for (int k = 0; k < arbors[j]->GetNumberOfClusters(); ++k)
					{
						clusterJ = arbors[j]->GetCluster(k);
						if (clusterJ->GetMinLayer() == layerJ)
						{
							t2 += clusterJ->GetTiming();
							++helper;
						}
					}



					t2 *= 1/(double)helper;

					double dtMax = timing * distanceSeed(arbors[j], arbors[i], 1., 1., 1.);
					double dtMin = timingMin * distanceSeed(arbors[j], arbors[i], 0.1, 0.1, 1.);

					double dt = t1 - t2;

					if (methodCloseLayers == false)
					{
						
						double distanceMin = MinimalDistanceLines(directionI, directionJ, positionI, positionJ);
						double error = 50.;
						double diffZ = MinimalDistanceLinePoint(directionI, directionJ, positionI, positionJ, distanceMin, error);

						

						if (distanceMin < distanceMax && dt < dtMax && dt > dtMin && positionJ[2] + diffZ > -error && diffZ < 50.)
						{
							distanceMax = distanceMin;		
							arbors[i]->SetArborId(j);
							
						}
					
					}else
					{
						Cluster* clusterJ;
						for (int k = 0; k < arbors[j]->GetNumberOfClusters(); ++k)
						{
							clusterJ = arbors[j]->GetCluster(k);
							if (clusterJ->GetMinLayer() < layerJ+3)
							{
								double distance = clusterI->Distance(clusterJ, 1., 1., 0.);

								if (diffLayer <= maxLayerForFusion && distance < distanceMax && dt < dtMax && dt > dtMin )
								{
									distanceMax = distance;
									arbors[i]->SetArborId(j);
								}
							}
						}
					}
				}	
			}
		}
	}

	arbors = CleanArbor(arbors);
	return arbors;
}





std::vector<Arbor*> FusionArbor(std::vector<Arbor*> arbors, double maxDistance, int minSize, int maxSize, double timing, double timingMin)
{
	for (int i = 0; i < arbors.size(); ++i)
	{
		arbors[i]->SetArborId(i);
	}

	IdentifieChargedParticle(arbors);

	for (int i = 0; i < arbors.size(); ++i)
	{
		int size = arbors[i]->GetSize();

		if (size < maxSize && size > minSize) 
		{
			double minDistanceBetweenArbors = 1000000.;
			int arborId = -1;

			arbors[i]->SetPosition();
			double  positionArbor[3];
			double directionArbor[3];

			for (int j = 0; j < 3; ++j)
			{
				directionArbor[j] = arbors[i]->GetDirection()[j];
			}

			for (int j = 0; j < arbors.size(); ++j)
			{	
				double distance;
				double realDistance;



					if (i != j)	
					{

						if (arbors[i]->IsCharged() == true)
						{
							double energy = arbors[i]->GetEnergy() + arbors[j]->GetEnergy();
							if (energy > 30.+sqrt(30.)*0.5 )	
							{
//!!!!!!!!!  30. IS BECAUSE CHARGED PARTICLE IS 30 GEV !!!
								continue;
							}
						}
						

						
						int seedCorrelation = 0;
						int sizeJ =  arbors[j]->GetSize();

						for (int k = 0 ; k < sizeJ; ++k)
						{
							Cell* cell = arbors[j]->GetCell(k);
							double* positionHit = cell->GetPosition();
							double dz = positionArbor[2] -  positionHit[2];

							realDistance = arbors[i]->Distance(cell);


							if (dz >= 0.)
							{
								arbors[i]->SetPosition(directionArbor[0], directionArbor[1], -dz);
								distance = arbors[i]->Distance(cell, 1., 1., 0.1);
										
								if (distance < maxDistance)
								{
									++seedCorrelation;
								}

							}
						}

						int seedCorrelationMax;

						
						if (size < 20)
						{
							seedCorrelationMax = 1;
						}else if (size < 40)
						{
							seedCorrelationMax = 2;
						}else if (size < 60)
						{
							seedCorrelationMax = 4;
						}else if (size < 100)
						{
							seedCorrelationMax = 4;
						}else 
						{
							seedCorrelationMax = 4;
						}
					

						if (seedCorrelation > seedCorrelationMax)
						{
							double distanceBetweenArbors = distance;
							if (distanceBetweenArbors < minDistanceBetweenArbors)
							{
								distanceBetweenArbors = minDistanceBetweenArbors;
								arborId = j;
							}
						}
					}
			}

			if (arborId != -1)
			{
				arbors[i]->SetArborId(arborId);
			}else
			{
			arbors[i]->SetArborId(i);
			}

		}
	}


	arbors = CleanArbor(arbors);



	return arbors;
}



std::vector<Arbor*> FusionArborFromSeed(std::vector<Arbor*> arbors, double maxDistance, int minSize, int maxSize, double timing, double timingMin)
{
	for (int i = 0; i < arbors.size(); ++i)
	{
		arbors[i]->SetArborId(i);
	}
	IdentifieChargedParticle(arbors);

	for (int i = 0; i < arbors.size(); ++i)
	{
		int size = arbors[i]->GetSize();

		Cluster* cluster = GetSeedCluster(arbors[i]);
		double t1 = cluster->GetTiming();

		if (size < maxSize && size > minSize) 
		{
			double minDistanceBetweenArbors = 1000000.;
			int arborId = -1;

			arbors[i]->SetPosition();
			double  positionArbor[3];
			double directionArbor[3];

			for (int j = 0; j < 3; ++j)
			{
				positionArbor[j] = cluster->GetPosition()[j];
				directionArbor[j] = arbors[i]->GetDirection()[j];
			}

			for (int j = 0; j < arbors.size(); ++j)
			{	
				int iId = arbors[i]->GetMinLayer();
				int jId = arbors[j]->GetMinLayer();
				if (iId < jId) continue;

				double distance;
				double realDistance;

					if (i != j)	
					{

						if (arbors[i]->IsCharged() == true)
						{
							double energy = arbors[i]->GetEnergy() + arbors[j]->GetEnergy();
							if (energy > 30.+sqrt(30.)*0.5 )	
							{
//!!!!!!!!!  30. IS BECAUSE CHARGED PARTICLE IS 30 GEV !!!
								continue;
							}
						}


						int seedCorrelation = 0;
						int sizeJ =  arbors[j]->GetSize();

						for (int k = 0 ; k < sizeJ; ++k)
						{
							Cell* cell = arbors[j]->GetCell(k);
							double* positionHit = cell->GetPosition();
							double dz = positionArbor[2] -  positionHit[2];

							realDistance = cluster->Distance(cell);


							if (dz >= 0.)
							{
								cluster->SetPosition(directionArbor[0], directionArbor[1], -dz);
								distance = cluster->Distance(cell, 1., 1., 0.1);

								double dt = t1 - cell->GetTiming();
								double dtMax = timing * realDistance;
								double dtMin = timingMin * realDistance;
										
								if (distance < maxDistance && dt < dtMax && dt > dtMin)
								{
									++seedCorrelation;
								}

							}
						}

						int seedCorrelationMax;

						
						if (size < 20)
						{
							seedCorrelationMax = 1;
						}else if (size < 40)
						{
							seedCorrelationMax = 2;
						}else if (size < 60)
						{
							seedCorrelationMax = 4;
						}else if (size < 100)
						{
							seedCorrelationMax = 4;
						}else 
						{
							seedCorrelationMax = 4;
						}
					

						if (seedCorrelation > seedCorrelationMax)
						{
							double distanceBetweenArbors = distance;
							if (distanceBetweenArbors < minDistanceBetweenArbors)
							{
								distanceBetweenArbors = minDistanceBetweenArbors;
								arborId = j;
							}
						}
					}
			}

			if (arborId != -1)
			{
				arbors[i]->SetArborId(arborId);
			}else
			{
			arbors[i]->SetArborId(i);
			}

		}
	}


	arbors = CleanArbor(arbors);



	return arbors;
}



std::vector<Arbor*> ClusteringSmallArbor(std::vector<Arbor*> arbors, double sizeMaxFood, double sizeMinEater, double maxRange, double timing, double timingMin)
{
	//processus vers la fin, tente d'associer les hit isolé et les petits arbres au particules reconstruite, prend uniquement la distance en compte
	int nbArbor = arbors.size();

	for (int i = 0; i < nbArbor; ++i)
	{
		arbors[i]->SetDirection();
		arbors[i]->SetArborId(i);
		arbors[i]->SetTiming();
	}
	

	for (int i = 0; i < nbArbor; ++i)
	{
		double minParam = 10000000.;
		int jj = -1;


		int sizeArbor1 = arbors[i]->GetSize();
		if (sizeArbor1 <= sizeMaxFood)
		{
			for (int h = 0; h < sizeArbor1; ++h)
			{
				
				Cell* cell1 = arbors[i]->GetCell(h);

				int layer = cell1->GetId()[2];

				for (int j = 0; j < nbArbor; ++j)
				{
					int sizeArbor2 = arbors[j]->GetSize();
					if (i != j && sizeArbor2 > sizeMinEater)
					{

						for (int k = 0; k < sizeArbor2; ++k)
						{
							Cell* cellTested = arbors[j]->GetCell(k);

							if (layer <= cellTested->GetId()[2]) continue;
					
							double distance = cell1->Distance(cellTested, 1., 1., 0.33);
							double range = maxRange;
							double dtMax = (timing)*(cell1->Distance(cellTested, 1., 1., 1.)+11.);
							double dtMin = (timingMin)*(cell1->Distance(cellTested, 0., 0., 1.)-11.);
							double dt = -cellTested->GetTiming() + cell1->GetTiming();

							if (distance < range && dt < dtMax && dt > dtMin)
							{ 
								double param = pow(distance/maxRange, 2);

								if (param < minParam)
								{
									jj = j;
									minParam = param;	
								}
							}
						}
					}
				}
			}
		}


		if (jj != -1)
		{
			arbors[i]->SetArborId(jj);
		}
	}
	arbors = CleanArbor(arbors);

	return arbors;
}



std::vector<Arbor*> RepairTimingFragmentation(std::vector<Arbor*> arbors, double dMax)
{
	for (int i = 0; i < arbors.size(); ++i)
	{
		arbors[i]->SetArborId(i);
	}


	int nbOfArbor = arbors.size();
	for (int i = 0; i < nbOfArbor; ++i)
	{

		Cluster* c1 = GetSeedCluster(arbors[i]);

		for (int j = 0; j < nbOfArbor; ++j)
		{
			if (i == j) continue;

			for (int k = 0; k < arbors[j]->GetSize(); ++k)
			{
				Cell* c2 = arbors[j]->GetCell(k);

				double distance = c1->Distance(c2);

				if (distance < dMax)
				{
					arbors[i]->SetArborId(j);		
				}
			}
		}
	}
	
	arbors = CleanArbor(arbors);

	return arbors;
}






std::vector<Arbor*> FinalRecombination(std::vector<Arbor*> arbors, double energyMax, double range) 
{
	IdentifieChargedParticle(arbors);
	if (arbors.size() == 1) return arbors;

	int id;


	int oldSize = arbors.size();
	int size = 0;


	


	while (size != oldSize)
	{
		SetArborId(arbors);
		id =  arborsMaxSizeTop3(arbors)[0];
		oldSize = size;
		for (int i = 0; i < arbors.size(); ++i)
		{	
			if(i == id) continue;
	
			double energy = arbors[i]->GetEnergy() + arbors[id]->GetEnergy();
			double dz = arbors[i]->GetPosition()[2] - arbors[id]->GetPosition()[2];

			if (dz > 0 && arbors[id]->Distance(arbors[i], 1., 1., 0.025) < range && energy < energyMax+sqrt(energyMax)*0.5 && arbors[i]->IsCharged() == false)
			{
				arbors[i]->SetArborId(id);
			}
		}
		arbors = CleanArbor(arbors);
		size = arbors.size();
	}
	SetArborId(arbors);

	if (arbors.size() < 3) return arbors;

	id =  arborsMaxSizeTop3(arbors)[1];

	if (arbors[id]->GetSize() < 30) return arbors;

	oldSize = arbors.size();
	size = 0;

		for (int i = 0; i < arbors.size(); ++i)
		{	
			if(i == id) continue;
	
			double energy = arbors[i]->GetEnergy() + arbors[id]->GetEnergy();
			double dz = arbors[i]->GetPosition()[2] - arbors[id]->GetPosition()[2];

			if (dz > 0 && arbors[id]->Distance(arbors[i], 1., 1., 0.025) < range && energy < energyMax+sqrt(energyMax)*0.5 && arbors[i]->IsCharged() == false)
			{
				arbors[i]->SetArborId(id);
			}
		}
	arbors = CleanArbor(arbors);
	SetArborId(arbors);
	
	if (arbors.size() < 4) return arbors;

	id =  arborsMaxSizeTop3(arbors)[2];

	if (arbors[id]->GetSize() < 30) return arbors;

	oldSize = arbors.size();
	size = 0;

		for (int i = 0; i < arbors.size(); ++i)
		{	
			if(i == id) continue;
	
			double energy = arbors[i]->GetEnergy() + arbors[id]->GetEnergy();
			double dz = arbors[i]->GetPosition()[2] - arbors[id]->GetPosition()[2];

			if (dz > 0 && arbors[id]->Distance(arbors[i], 1., 1., 0.025) < range && energy < energyMax+sqrt(energyMax)*0.5 && arbors[i]->IsCharged() == false)
			{
				arbors[i]->SetArborId(id);
			}
		}
		arbors = CleanArbor(arbors);


	return arbors;
}




std::vector<Arbor*> MakeArbors(std::vector<Arbor*> branches, std::vector<Cluster*> trace)
{
	std::vector<Arbor*> arbors;
	int size;
	int oldSize;

	double timingMin = .0032;
	double timingMax = .008;


	CreateConnectionsBetweenBranchesAndtrace(branches, trace);
	arbors = MakeFirstArbors(branches, trace);


	oldSize = arbors.size();
	size = 0;
	while (size != oldSize)
	{
		oldSize = size;
		// arbors = RepairTimingFragmentation(arbors, 30.);
		size = arbors.size();
	}

	arbors = CommonSeedArbor(arbors, 35., 1, 80);

	arbors = BuildArborSeedMissing(arbors, true, 3, 15., 10, 10000, timingMax, timingMin);
	arbors = BuildArborSeedMissing(arbors, false, 3, 15., 10, 10000, timingMax, timingMin);





	oldSize = arbors.size();
	size = 0;
	while (size != oldSize)
	{
		oldSize = size;
		arbors = FusionArbor(arbors, 25., 5, 1000, timingMax, timingMin);
		size = arbors.size();
	}


	oldSize = arbors.size();
	size = 0;
	while (size != oldSize)
	{
		oldSize = size;
		arbors = ClusteringSmallArbor(arbors, 8, 25, 65., timingMax, timingMin);
		size = arbors.size();
	}


	oldSize = arbors.size();
	size = 0;
	while (size != oldSize)
	{
		oldSize = size;
		arbors = FusionArbor(arbors, 20., 20, 1000, timingMax, timingMin);
		size = arbors.size();
	}

	oldSize = arbors.size();
	size = 0;
	while (size != oldSize)
	{
		oldSize = size;
		arbors = FusionArborFromSeed(arbors, 25., 20, 1000, timingMax, timingMin);
		size = arbors.size();
	}

	
	oldSize = arbors.size();
	size = 0;
	while (size != oldSize)
	{
		oldSize = size;
		arbors = ClusteringSmallArbor(arbors,30, 30, 300., 10, -10);
		size = arbors.size();
	}

	arbors = FinalRecombination(arbors, 30., 50.);
	arbors = FinalRecombination(arbors, 30., 100.);
	arbors = FinalRecombination(arbors, 30., 150.);

	// optimalData(arbors);
	// arbors = CleanArbor(arbors);
	return arbors;

	






	












	return arbors;
}


#endif