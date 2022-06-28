#include "HitGroupingByLayer.h"
#include "Tools.C"
#include "parameter.C"

double xydistance = 10.5;
double zdistance = 26.5;
int seuil = 10;

double coef1 = .75;	//range layer by layer
double coef2 = .5;	//range cluster to cluster
double dtCoef = 1000000;	// big don't take time in consideretion
double weightForward = 1.;
double weightBackward = 1.;




std::vector<Cluster*>  MakeClustersByLayer(Map* map)
{
	//cherche les hits voisins sur une couche
	Cell* cellTested;
	std::vector<Cluster*> clustersByLayer;
	bool clusterised, hited;

	for (int i = 0; i < map->GetHits().size(); ++i)
	{
		cellTested = map->GetCell(map->GetHits()[i]->GetID()[0], map->GetHits()[i]->GetID()[1], map->GetHits()[i]->GetID()[2]);
		clusterised = cellTested->GetClusterised();
		hited = cellTested->GetHited();

		if (clusterised == false && hited == true)
		{
			Cluster* cluster = new Cluster();
			cluster->Clustering(cellTested, map);
			int clusterId = clustersByLayer.size();
			cluster->SetClusterId(clusterId);
			cluster->SetPosition();
			clustersByLayer.push_back(cluster);
		}
	}

	return clustersByLayer;
}


std::vector<Cluster*> SelectClustersForPrimaryTrace(std::vector<Cluster*> clustersByLayer, double timeMin, double timeMax)
{
	//selectionne des données selon un critère a définir (timing)
	std::vector<Cluster*> potentialPrimaryTrace; 

	for (int i = 0; i < clustersByLayer.size(); ++i)
	{
		int layerId = clustersByLayer[i]->GetMinLayer();
		int sizeCluster = clustersByLayer[i]->GetSize();
		double time = clustersByLayer[i]->GetTime();

		if (time < timeMax && time > timeMin)
		{
			potentialPrimaryTrace.push_back(clustersByLayer[i]);
		}
	}

	return potentialPrimaryTrace;
}

void CreateConnectionsForPrimaryTrace(std::vector<Cluster*>potentialPrimaryTrace, double coefSpace, double coefTime)
{
	//création de connexion entre les cluster(groupe de hit sur une meme couche) selon un distance max coefSpace et une distance temporel max coefTime
	for (int i = 0; i < potentialPrimaryTrace.size(); ++i)
	{
		potentialPrimaryTrace[i]->SetTime();
		potentialPrimaryTrace[i]->SetPosition();
		potentialPrimaryTrace[i]->ClearConnectionBackward();
		potentialPrimaryTrace[i]->ClearConnectionForward();
	}



	for (int i = 0; i < potentialPrimaryTrace.size(); ++i)
	{

		for (int j = 0; j < potentialPrimaryTrace.size(); ++j)
		{
			int layerI = potentialPrimaryTrace[i]->GetMinLayer();
			int layerJ = potentialPrimaryTrace[j]->GetMinLayer();
			int diffLayer = layerJ - layerI;

			if (layerI < layerJ && diffLayer <= 3)
			{
				double distance = potentialPrimaryTrace[i]->Distance(potentialPrimaryTrace[j]);
				double dtMax = coefTime * potentialPrimaryTrace[i]->Distance(potentialPrimaryTrace[j],1, 1, 0.5);
				double dt = abs(potentialPrimaryTrace[j]->GetTime() - potentialPrimaryTrace[i]->GetTime());

				double range = coefSpace;

				

				if (distance < range && dt < dtMax)
				{ 
					potentialPrimaryTrace[j]->AddConnectionBackward(i);
					potentialPrimaryTrace[i]->AddConnectionForward(j);
				}
			}
		}
	}
}


void CleanConnection(std::vector<Cluster*>potentialPrimaryTrace, double maxSpace, double maxTime)
{
	// garde la connexion vers l'arriere la plus probable
	for (int i = 0; i < potentialPrimaryTrace.size(); ++i)
	{
		double orientation[3] = {0, 0, 0};
		std::vector<int> bck = potentialPrimaryTrace[i]->GetConnectionBackward();

		

		if (bck.size() > 1)
		{
			std::vector<int> fwd = potentialPrimaryTrace[i]->GetConnectionForward();

			double *position = potentialPrimaryTrace[i]->GetPosition(); 

			for (int j = 0; j < bck.size(); ++j)
			{
				double *bckPosition = potentialPrimaryTrace[bck[j]]->GetPosition();
				double bckOrientation[3];

				for (int k = 0; k < 3; ++k)
				{
					bckOrientation[k] = bckPosition[k] - position[k];
				}

				double norm = Norm(bckOrientation);

				for (int k = 0; k < 3; ++k)
				{
					orientation[k] += 2*bckOrientation[k]/norm;
				}
			}

			for (int j = 0; j < fwd.size(); ++j)
			{
				double *fwdPosition = potentialPrimaryTrace[fwd[j]]->GetPosition();
				double fwdOrientation[3];

				for (int k = 0; k < 3; ++k)
				{
					fwdOrientation[k] = position[k] - fwdPosition[k] ;
				}

				double norm = Norm(fwdOrientation);

				for (int k = 0; k < 3; ++k)
				{
					orientation[k] += 3*fwdOrientation[k]/norm;
				}
			}

			double param = 0.;
			int bckId = -1;

			for (int j = 0; j < bck.size(); ++j)
			{
				double *bckPosition = potentialPrimaryTrace[bck[j]]->GetPosition();
				double bckOrientation[3];



				for (int k = 0; k < 3; ++k)
				{
					bckOrientation[k] =  position[k] - bckPosition[k];
				}

				double paramAngle = acos(CosAngleBetweenVectors(orientation, bckOrientation))/M_PI;
				double paramDistance = (potentialPrimaryTrace[bck[j]]->Distance(potentialPrimaryTrace[i]))/maxSpace;
				double paramTiming = abs(potentialPrimaryTrace[i]->GetTime() - potentialPrimaryTrace[bck[j]]->GetTime())/maxTime;

				double newParam = pow(paramDistance, 5);// * pow(paramTiming, 4);

				if (newParam > param)
				{
					param = newParam;
					bckId = bck[j];

				}
			}	

			if (bckId != -1)
			{
				potentialPrimaryTrace[i]->ClearConnectionBackward();		
				potentialPrimaryTrace[i]->AddConnectionBackward(bckId);	
			}
			
		}

	}

	for (int i = 0; i <  potentialPrimaryTrace.size(); ++i)
	{
		potentialPrimaryTrace[i]->ClearConnectionForward();
	}

	for (int i = 0; i <  potentialPrimaryTrace.size(); ++i)
	{
		if (potentialPrimaryTrace[i]->GetConnectionBackward().size() > 0)
		{
			int backwardId = potentialPrimaryTrace[i]->GetConnectionBackward()[0];
			potentialPrimaryTrace[backwardId]->AddConnectionForward(backwardId);
		}
	}
}




void CreateSecondaryConnection(std::vector<Cluster*>potentialPrimaryTrace, double maxSpace, double maxTime)
{
	//utilise les premieres connexions crée pour regarder si d'autre connexion peuvent être creer dans la direction de la premiere, le paterne a des similitudes avec les fonction précédente mais elle est indépendante (je doit corriger ca quand tout marchera)

	std::vector<int> connections;
	for (int i = 0; i < potentialPrimaryTrace.size(); ++i)
	{
		if (potentialPrimaryTrace[i]->GetConnectionBackward().size() > 0) 
		{
			int id = potentialPrimaryTrace[i]->GetConnectionBackward()[0];
			connections.push_back(id);
		}else
		{
			connections.push_back(-1);
		}
	}

	//clear les vecteurs


	for (int i = 0; i < potentialPrimaryTrace.size(); ++i)
	{
		potentialPrimaryTrace[i]->ClearConnectionBackward();
		potentialPrimaryTrace[i]->ClearConnectionForward();
		potentialPrimaryTrace[i]->SetSelectedBackwardConnection(connections[i]);
	}

	//creer les connection


	double maxAngle = 0.7; //rad

	for (int i = 0; i < potentialPrimaryTrace.size(); ++i)
	{
	int layerI = potentialPrimaryTrace[i]->GetMinLayer();

		for (int j = 0; j < potentialPrimaryTrace.size(); ++j)
		{
			
			int layerJ = potentialPrimaryTrace[j]->GetMinLayer();
			int diffLayer = layerJ - layerI;
			int id  = potentialPrimaryTrace[i]->GetSelectedBackwardConnection();

			if (layerI < layerJ && diffLayer <= 5)
			{
				double distance = potentialPrimaryTrace[i]->Distance(potentialPrimaryTrace[j]);
				double range = maxSpace;

				double dtMax = maxTime * potentialPrimaryTrace[i]->Distance(potentialPrimaryTrace[j]);
				double dt = abs(potentialPrimaryTrace[j]->GetTime() - potentialPrimaryTrace[i]->GetTime());
				
				double fwdOrientation[3];
				double orientation[3];		

				if (id == -1) 
				{
					orientation[0] = 0;
					orientation[0] = 0;
					orientation[0] = 1;
				}

				
				for (int k = 0; k < 3; ++k)
				{
					
					if (id != -1) orientation[k] = potentialPrimaryTrace[i]->GetPosition()[k] - potentialPrimaryTrace[id]->GetPosition()[k];
					fwdOrientation[k] = potentialPrimaryTrace[j]->GetPosition()[k] - potentialPrimaryTrace[i]->GetPosition()[k];
				}


				double angle = acos(CosAngleBetweenVectors(orientation, fwdOrientation));
				
				if (distance < range && dt < dtMax && angle < maxAngle)
				{ 
					potentialPrimaryTrace[j]->AddConnectionBackward(i);
					potentialPrimaryTrace[i]->AddConnectionForward(j);
				}

				
			}
		}
	}

	// Regarde les vecteurt derrier et devant ceux de lka connextion etudier et cherche dans cette zone
	//de layer en layer pour ne garder que els bonnes connection vers l'avant


	for (int layer = 47; layer > 0; --layer)
	{
		
	

		for (int i = 0; i < potentialPrimaryTrace.size(); ++i)
		{
			

			int connexionUsed = potentialPrimaryTrace[i]->GetSelectedBackwardConnection();
			int layerI = potentialPrimaryTrace[i]->GetMinLayer();

			

			if (connexionUsed != -1 && layerI == layer)
			{
				double orientation[3] = {0, 0, 0};
				double *position = potentialPrimaryTrace[i]->GetPosition(); 

				std::vector<int> bck = potentialPrimaryTrace[connexionUsed]->GetConnectionBackward();
				std::vector<int> fwd = potentialPrimaryTrace[i]->GetConnectionForward();

					

					for (int j = 0; j < bck.size(); ++j)
					{
						double *bckPosition = potentialPrimaryTrace[bck[j]]->GetPosition();
						double bckOrientation[3];

						for (int k = 0; k < 3; ++k)
						{
							bckOrientation[k] = bckPosition[k] - position[k];
						}

						double norm = Norm(bckOrientation);

						for (int k = 0; k < 3; ++k)
						{
							orientation[k] += .2*bckOrientation[k]/norm;
						}
					}

					for (int j = 0; j < fwd.size(); ++j)
					{
						double *fwdPosition = potentialPrimaryTrace[fwd[j]]->GetPosition();
						double fwdOrientation[3];

						for (int k = 0; k < 3; ++k)
						{
							fwdOrientation[k] = position[k] - fwdPosition[k] ;
						}

						double norm = Norm(fwdOrientation);

						for (int k = 0; k < 3; ++k)
						{
							orientation[k] += 5*fwdOrientation[k]/norm;
						}
					}

				if (orientation[2] >= 0)
				{
					double param = 0.;
					int bckId = -1;

					for (int j = 0; j < bck.size(); ++j)
					{
						double *bckPosition = potentialPrimaryTrace[bck[j]]->GetPosition();
						double bckOrientation[3];

						double paramAngle = acos(CosAngleBetweenVectors(orientation, bckOrientation))/M_PI;
						double paramDistance = (potentialPrimaryTrace[bck[j]]->Distance(potentialPrimaryTrace[i]))/maxSpace;
						double paramTiming = abs(potentialPrimaryTrace[i]->GetTime() - potentialPrimaryTrace[bck[j]]->GetTime())/maxTime;

						double newParam = pow(paramAngle, 1) * pow(paramDistance, 5) /** pow(paramTiming, 10)*/;


						if (maxAngle < maxAngle && param < newParam)
						{
							param = newParam;
							bckId = bck[j];

						}
					}	

					if (bckId != -1)
					{

						for (int j = 0; j < bck.size(); ++j)
						{
							if (bck[j] != bckId)
							{
								std::vector<int> tool = potentialPrimaryTrace[bck[j]]->GetConnectionForward();
								potentialPrimaryTrace[bck[j]]->ClearConnectionForward();
								for (int k = 0; k < tool.size(); ++k)
								{
									if (bck[j] != tool[k])
									{
										potentialPrimaryTrace[bck[j]]->AddConnectionForward(bck[j]);
									}
								}
							}
							
						}

						potentialPrimaryTrace[connexionUsed]->ClearConnectionBackward();		
						potentialPrimaryTrace[connexionUsed]->AddConnectionBackward(bckId);	
					}
				}
			}
		}
	}
	for (int i = 0; i < potentialPrimaryTrace.size(); ++i)
	{
		if (potentialPrimaryTrace[i]->GetConnectionBackward().size() == 0 && connections[i] != -1)
		{
			potentialPrimaryTrace[i]->AddConnectionBackward(connections[i]);	
		}
	}


}







std::vector<Arbor*> BuildPotentialPrimaryTrace(std::vector<Cluster*> potentialPrimaryTrace, std::vector<Cluster*> clustersByLayer)
{
	//regroupe les cluster lier par des connexions entres eux au seins d'une même arbre
	std::vector<Arbor*> primaryTrace;

	for (int i = 0; i < potentialPrimaryTrace.size(); ++i)
	{
		int nbBackwardConnections = potentialPrimaryTrace[i]->GetConnectionBackward().size();

		if (nbBackwardConnections == 0)
		{
			int arborId = primaryTrace.size();
			potentialPrimaryTrace[i]->SetArborID(arborId);
			Arbor* newArbor = new Arbor();
			newArbor->AddCluster(potentialPrimaryTrace[i]);
			newArbor->SetArborId(arborId);
			primaryTrace.push_back(newArbor);
		}
	}
	for (int i = 0; i < potentialPrimaryTrace.size(); ++i)
	{
		int nbBackwardConnections = potentialPrimaryTrace[i]->GetConnectionBackward().size();

		if (nbBackwardConnections > 0)
		{
			int clusterBackwardId = potentialPrimaryTrace[i]->GetConnectionBackward()[0];
			int arborId = clustersByLayer[clusterBackwardId]->GetArborID();

			while (arborId == -1)
			{
				clusterBackwardId = clustersByLayer[clusterBackwardId]->GetConnectionBackward()[0];
				arborId = clustersByLayer[clusterBackwardId]->GetArborID();
			}

			if (arborId >= 0 && arborId < primaryTrace.size())
			{
				primaryTrace[arborId]->AddCluster(potentialPrimaryTrace[i]);
			}
		}
	}

	return primaryTrace;
}

std::vector<Arbor*>  CleanPotentialPrimaryTrace(std::vector<Arbor*> primaryTrace)
{
	//OBSELETE: la fonction précédente etait créait des datas en double, je les effacais ici (je la garde sous la main au cas ou)
	int nbOfTrace = primaryTrace.size();

	for (int i = nbOfTrace-1; i >= 0; --i)
	{
		int numberOfLayer = 0;
		std::vector<int> layerMemory;
		int numberOfCluster = primaryTrace[i]->GetNumberOfClusters();

		for (int j = 0; j < numberOfCluster; ++j)
		{
			Cluster* cluster = primaryTrace[i]->GetCluster(j);
			int layerId = cluster->GetMinLayer();
			bool layerAlreadyExists = false;

			for (int k = 0; k < layerMemory.size(); ++k)
			{
				if (layerId == layerMemory[k])
				{
					layerAlreadyExists = true;
				}
			}
			
			if (layerAlreadyExists == false)
			{
				++numberOfLayer;
				layerMemory.push_back(layerId);
			}
		}
		if (numberOfLayer <= -1)
		{
			primaryTrace.erase (primaryTrace.begin()+i);
		}
	}
	return primaryTrace;
}


void IdentifieChargedParticle(std::vector<Arbor*> arbors)
{
	//essaye d'identifier les particules gerbe issus des particules charger (je l'utilise pas)
	for (int i = 0; i < arbors.size(); ++i)
	{
		int countLayerHited = 0;
		int size = arbors[i]->GetNumberOfClusters();


		for (int j = 0; j < size; ++j)
		{
			Cluster* cluster = arbors[i]->GetCluster(j);
			int nbOfConnexionsForward = cluster->GetConnectionForward().size();
			int layerId = cluster->GetMinLayer();
			int sizeCluster = cluster->GetSize();

			if (layerId <= 4)
			{
				++countLayerHited;
			}
		}
	
		if ((countLayerHited >= 3 && size > 20) || (size > 300))
		{
			arbors[i]->SetCharged(true);
		}else
		{
			arbors[i]->SetCharged(false);
		}
	}
}


std::vector<Arbor*> BuildPrimaryTrace(std::vector<Cluster*> clustersByLayer, double coefSpace, double coefTime, double timeMin, double timeMax, bool reset)
{
	//premiere fonction mère, elle utilise une partie des fonction ci dessus pour creer les premiers arbres
	std::vector<Cluster*> potentialPrimaryTrace; 
	std::vector<Arbor*> primaryTrace;

	for (int i = 0; i < clustersByLayer.size(); ++i) clustersByLayer[i]->SetTime();

	potentialPrimaryTrace = SelectClustersForPrimaryTrace(clustersByLayer, timeMin, timeMax);
	CreateConnectionsForPrimaryTrace(potentialPrimaryTrace, coefSpace, coefTime);
	if (reset == false) CreateSecondaryConnection(potentialPrimaryTrace, 1.5 * coefSpace, coefTime);
	CleanConnection(potentialPrimaryTrace, coefSpace, coefTime);

	

	primaryTrace = BuildPotentialPrimaryTrace(potentialPrimaryTrace, potentialPrimaryTrace);


	primaryTrace = CleanPotentialPrimaryTrace(primaryTrace);
	
	



	return primaryTrace;
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


std::vector<Arbor*> BuildArborSeedMissing(std::vector<Arbor*> arbors, bool methodCloseLayers, int maxLayerForFusion, double range, double minSize, double maxSize)
{
	//deux fonction utiles pour les arbres neutre notement, depend de la valeur du booleen "methodCloseLayers":
	
	//		1)si un arbres a sa racine sur une couche, cherche dans un petit perimetre si un autres arbres a également sa racine ici pour les associer
	
	//		2)a partie de l'orientation de deux arbres, regarde si ils sembles pouvoir venir d'une même région en amont

	std::vector<Cluster*> potentialPrimaryTrace; 

	for (int i = 0; i < arbors.size(); ++i)
	{
		arbors[i]->SetDirection();
		arbors[i]->SetArborId(i);
		arbors[i]->SetTime();
	}

	for (int i = arbors.size() - 1; i > 0 ; --i)
	{

		int layerI = arbors[i]->GetMinLayer();

		if (arbors[i]->GetMaxLayer() - layerI < 4) continue;

		double* directionI = arbors[i]->GetDirection();
		double* positionI = arbors[i]->GetPosition();

		if (arbors[i]->GetSize() > maxSize || arbors[i]->GetSize() < minSize)
		{
			continue;
		}

		Cluster* clusterI;
		for (int k = 0; k < arbors[i]->GetNumberOfClusters(); ++k)
		{
			clusterI = arbors[i]->GetCluster(k);
			if (clusterI->GetMinLayer() == layerI)
			{
				break;
			}	
		}

		for (int j = i-1; j >= 0; --j)
		{
			int layerJ = arbors[j]->GetMinLayer();
			int diffLayer = abs(layerI - layerJ);
			double distanceMax = range;
			double* directionJ = arbors[j]->GetDirection();
			double* positionJ = arbors[j]->GetPosition();
			double dtMax = .01*(double)(diffLayer + 1)*zdistance;;
			double dt = abs(arbors[j]->GetTime() - arbors[i]->GetTime());

			if (methodCloseLayers == false)
			{
				
				double distanceMin = MinimalDistanceLines(directionI, directionJ, positionI, positionJ);
				double error = 50.;
				double diffZ = MinimalDistanceLinePoint(directionI, directionJ, positionI, positionJ, distanceMin, error);

				//double cosTheta = CosAngleBetweenVectors(directionI, directionJ);

				if (distanceMin < distanceMax && dt < dtMax && positionJ[2] - diffZ >= -error)
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
					if (clusterJ->GetMinLayer() == layerJ)
					{
						double distance = clusterI->DistanceXY(clusterJ);

						if (diffLayer < maxLayerForFusion && distance < distanceMax)
						{
							distanceMax = distance;
							arbors[i]->SetArborId(j);
						}
					}	
				}
			}
		}
	}

	arbors = CleanArbor(arbors);
	return arbors;
}



std::vector<Arbor*> FusionArbor(std::vector<Arbor*> primaryTrace, bool FusionNeutralOnly, double distanceConditon, int minSize, int maxSize, int thr = 1)
{
	// a partir de la direction de l'arbre, regarde si il peut venir d'un arbes derrière lui
	for (int i = 0; i < primaryTrace.size(); ++i)
	{
		primaryTrace[i]->SetArborId(i);
	}

	int maxDistance =  distanceConditon*sqrt(2.)*xydistance;

	for (int i = 0; i < primaryTrace.size(); ++i)
	{
		int size = primaryTrace[i]->GetSize();
		int threshold = primaryTrace[i]->GetThreshold(thr);

		if (size < maxSize && size > minSize && threshold >= 0) 
		{
			double minDistanceBetweenArbors = 1000000.;
			int arborId = -1;

			primaryTrace[i]->SetPosition();
			double  positionArbor[3];
			double directionArbor[3];

			for (int j = 0; j < 3; ++j)
			{
				positionArbor[j] = primaryTrace[i]->GetPosition()[j];
				directionArbor[j] = primaryTrace[i]->GetDirection()[j];
			}

			for (int j = 0; j < primaryTrace.size(); ++j)
			{	
				int iId = primaryTrace[i]->GetMinLayer();
				int jId = primaryTrace[j]->GetMinLayer();
				if (iId < jId)
				{
					continue;
				}

				double distanceAAA = primaryTrace[i]->Distance(primaryTrace[j]);
				/*if ((primaryTrace[i]->IsCharged() == true && FusionNeutralOnly == true ))
				{
					continue;
				}*/
					if (i != j)	
					{
						int seedCorrelation = 0;
						int sizeJ =  primaryTrace[j]->GetSize();

						for (int k = 0 ; k < sizeJ; ++k)
						{
							Cell* cell = primaryTrace[j]->GetCell(k);
							double* positionHit = cell->GetPosition();
							//double dz = primaryTrace[i]->Distance(primaryTrace[j]);
							double dz = positionArbor[2] -  positionHit[2];
							if (dz > 0.)
							{
								primaryTrace[i]->SetPosition(directionArbor[0], directionArbor[1], dz);
								double distance = primaryTrace[i]->Distance(cell);

								double dtMax = dtCoef*dz + 787878;
								double dt = abs(primaryTrace[j]->GetTime() - primaryTrace[i]->GetTime());
									
								if (distance < maxDistance && dt < dtMax)
								{
									++seedCorrelation;
								}
							}
						}

						int seedCorrelationMax;

						/*
						if (size < 20)
						{
							seedCorrelationMax = 2;
						}else if (size < 40)
						{
							seedCorrelationMax = 3;
						}else if (size < 60)
						{
							seedCorrelationMax = 6;
						}else if (size < 100)
						{
							seedCorrelationMax = 10;
						}else 
						{
							seedCorrelationMax = 20;
						}
						*/
						seedCorrelationMax = 1;

						if (seedCorrelation > seedCorrelationMax)
						{
							double distanceBetweenArbors = primaryTrace[i]->Distance(primaryTrace[j]);
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
				primaryTrace[i]->SetArborId(arborId);
			}else
			{
			primaryTrace[i]->SetArborId(i);
			}

		}
	}


	primaryTrace = CleanArbor(primaryTrace);



	return primaryTrace;
}


std::vector<Arbor*> FusionBigArbor(std::vector<Arbor*> arbors, double sizeMin, double distanceMaxXY, double coef)
{
	//	essaye d'utiliser l'enrientation et les angles entres les arbres pour faires des association (souvent inutiles mais pas toujours, souvent surpasser par "TryFusionArbor")
	int nbArbor = arbors.size();

	for (int i = 0; i < nbArbor; ++i)
	{
		arbors[i]->SetDirection();
		arbors[i]->SetArborId(i);
		arbors[i]->SetTime();
	}

	for (int i = 0; i < nbArbor; ++i)
	{
		int sizeI = arbors[i]->GetSize();

		if (sizeI > sizeMin)
		{
			double* directionI = arbors[i]->GetDirection();
			double* positionI = arbors[i]->GetPosition();
			int layerMinI = arbors[i]->GetMinLayer();
			int layerMaxI = arbors[i]->GetMaxLayer();

			for (int j = i + 1; j < nbArbor; ++j)
			{
				int sizeJ = arbors[j]->GetSize();

				if (sizeJ > sizeMin)
				{
					double* directionJ = arbors[j]->GetDirection();
					double* positionJ = arbors[j]->GetPosition();
					double ez[3] = {0., 0 , 1};
					int layerMinJ = arbors[j]->GetMinLayer();
					int layerMaxJ = arbors[j]->GetMaxLayer();

					Arbor tested;

					for (int k = 0; k < arbors[j]->GetNumberOfClusters(); ++k)
					{
						tested.AddCluster(arbors[j]->GetCluster(k));
					}
					for (int k = 0; k < arbors[i]->GetNumberOfClusters(); ++k)
					{
						tested.AddCluster(arbors[i]->GetCluster(k));
					}
					tested.SetDirection();
					double* directionIJ = tested.GetDirection();

					double cosIZ = CosAngleBetweenVectors(directionI, ez);
					double cosJZ = CosAngleBetweenVectors(directionJ, ez);
					double cosIJZ = CosAngleBetweenVectors(directionIJ, ez);
					double cosMean = (cosIZ + cosJZ)/2.;
					double dt = abs(arbors[j]->GetTime() - arbors[i]->GetTime());
					double dtMax = dtCoef*(arbors[i]->Distance(arbors[j]));

					double dist = DistanceXY(positionI, positionJ);

					if (cosIJZ > cosMean && cosIJZ > coef && dist < distanceMaxXY)
					{
						arbors[i]->SetArborId(j);
					}
				}
			}
		}
	}

	arbors = CleanArbor(arbors);

	return arbors;
}

std::vector<Arbor*> ClusteringSmallArbor(std::vector<Arbor*> arbors, double sizeMax, double maxRange, double coefTime)
{
	//processus vers la fin, tente d'associer les hit isolé et les petits arbres au particules reconstruite, prend uniquement la distance en compte
	int nbArbor = arbors.size();

	for (int i = 0; i < nbArbor; ++i)
	{
		arbors[i]->SetDirection();
		arbors[i]->SetArborId(i);
		arbors[i]->SetTime();
	}
	

	for (int i = 0; i < nbArbor; ++i)
	{
		double minParam = 10000000.;
		int jj = -1;

		if (arbors[i]->GetSize() <= sizeMax)
		{

			for (int j = 0; j < nbArbor; ++j)
			{
				int sizeArbor = arbors[j]->GetSize();
				if (i != j && sizeArbor > sizeMax )
				{

					for (int k = 0; k < sizeArbor; ++k)
					{
						Cell* cellTested = arbors[j]->GetCell(k);
											
						double distance = arbors[i]->Distance(cellTested, 1., 1., 0.5);
						double range = maxRange;
						double dtMax = coefTime*distance;
						double dt = abs(cellTested->GetTime() - arbors[i]->GetTime());

						if (distance < range && dt < dtMax)
						{ 
							double param = pow(distance/maxRange, 2) * pow(dt, 1);
							if (sizeArbor < sizeMax)	param = pow(param, 2); 

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


		if (jj != -1)
		{
			arbors[i]->SetArborId(jj);
		}
	}
	arbors = CleanArbor(arbors);

	return arbors;
}

std::vector<Arbor*> TryFusionArbor(std::vector<Arbor*> arbors, double energyParticle, double range)
{
	//Utilise l'angle ainsi que l'énergie pour tenter de recombiner des arbres entres eux
	int nbArbor = arbors.size();

	for (int i = 0; i < nbArbor; ++i)
	{
		arbors[i]->SetDirection();
		arbors[i]->SetArborId(i);
		arbors[i]->SetTime();
	}

	// AJOUTER UNE FONCTION QUI IDENTIFIE LA PARTICULES CHARGER (en fonction de la position initial etc...)
		int size = 0;
		int id = 0;

		for (int i = 0; i < arbors.size(); ++i)
		{
			if (arbors[i]->GetMinLayer() < 10 && arbors[i]->GetSize() > size)
			{
				size = arbors[i]->GetSize();
				id = i;
			}
		}
	double energyMax = energyParticle + sqrt(energyParticle)*0.5;

	if (energyParticle > energyMax) return arbors;
	


	for (int i = 0; i < nbArbor; ++i)
	{
		if (i != id) continue;


		double* directionI = arbors[i]->GetDirection();
		double* positionI = arbors[i]->GetPosition();
		int layerMinI = arbors[i]->GetMinLayer();
		int layerMaxI = arbors[i]->GetMaxLayer();


		for (int j = 0; j < nbArbor; ++j)
		{
			if (j == i) continue;

			double* directionJ = arbors[j]->GetDirection();
			double* positionJ = arbors[j]->GetPosition();
			double ez[3] = {0., 0 , 1};
			int layerMinJ = arbors[j]->GetMinLayer();
			int layerMaxJ = arbors[j]->GetMaxLayer();

			Arbor tested;

			for (int k = 0; k < arbors[j]->GetNumberOfClusters(); ++k)
			{
				tested.AddCluster(arbors[j]->GetCluster(k));
			}
			for (int k = 0; k < arbors[i]->GetNumberOfClusters(); ++k)
			{
				tested.AddCluster(arbors[i]->GetCluster(k));
			}
			tested.SetDirection();
			double* directionIJ = tested.GetDirection();

			double cosIZ = CosAngleBetweenVectors(directionI, ez);
			double cosJZ = CosAngleBetweenVectors(directionJ, ez);
			double cosIJZ = CosAngleBetweenVectors(directionIJ, ez);

			double energy = tested.GetEnergy();
			

			double dist = DistanceXY(positionI, positionJ);

			if (dist < range && 1.7*cosIJZ < cosIZ + cosJZ && energy < energyMax)
			{
				arbors[i]->SetArborId(j);
			}
		}
	}
	arbors = CleanArbor(arbors);

	return arbors;
}




std::vector<Arbor*> EnergyReconstruction(std::vector<Arbor*> arbors, double energyFromChargedParticle, double coefRangeMax, double size, double id)
{
	//OBSELETE: ne prend que l'énergie en compte, est surpasser par "TryFusionArbor"
	double energy;
	double energyLimit = energyFromChargedParticle + sqrt(energyFromChargedParticle)*0.5;
	double distanceMax = 1e20;
	std::vector<double> potentialParticles;
	
	int ii = -1;
	int jj = -1;

	for (int i = 0; i < arbors.size(); ++i)
	{
		arbors[i]->SetDirection();
		arbors[i]->SetArborId(i);
		if (arbors[i]->GetSize() > size)
		{
			potentialParticles.push_back(i);
		}
	}
	
	for (int i = 0; i < arbors.size(); ++i)
	{

		if (i != id) continue;

			double* positionI = arbors[i]->GetPosition();
			double Ei = arbors[i]->GetEnergy();

			if (arbors[i]->GetSize() < size ) continue;


			int minLayerI = arbors[i]->GetMinLayer();

			Cluster* clusterI;
			for (int j = 0; j < arbors[i]->GetNumberOfClusters(); ++j)
			{
				if(arbors[i]->GetCluster(j)->GetMinLayer() == minLayerI) clusterI = arbors[i]->GetCluster(j);
			}

			for (int j = 0; j < arbors.size(); ++j)
			{
				if (i == j || arbors[j]->GetSize() < size ) continue;

				int minLayerJ = arbors[j]->GetMinLayer();


				if (minLayerJ < minLayerI) continue;


				Cluster* clusterJ;
				for (int k = 0; k < arbors[j]->GetNumberOfClusters(); ++k)
				{
					
					if(arbors[j]->GetCluster(k)->GetMinLayer() == minLayerI) clusterI = arbors[j]->GetCluster(k);
				}


				double distance = clusterI->Distance(clusterJ, 1., 1., 0.2);

				double range = coefRangeMax;

				energy = arbors[j]->GetEnergy() + Ei;




				if (distance < range && distance < distanceMax && energy < energyLimit)
				{
					jj = j;
					ii = i;
					distanceMax = distance;
				}
			}
	}

	if (jj != -1 && ii != -1)
	{
		arbors[jj]->SetArborId(ii);
	}

	arbors = CleanArbor(arbors);


	return arbors;
}

int* arborsMaxSizeTop3(std::vector<Arbor*> arbors)
{

	int max = 0;
	int max2 = max;
	int id1 = -1;
	int id2 = -1;
	
	for (int i = 0; i < arbors.size(); ++i)	
	{
		arbors[i]->SetArborId(i);
		int size = arbors[i]->GetSize();
		if (size > max)	
		{
			max2 = max;
			id2 = id1;
			max = size;
			id1 = i;

		}else if (size > max2)
		{

			max2 = size;
			id2 = i;
		}
	}

	int id[2];

	id[0] = id1;
	id[1] = id2;
	return id;
}

double distance(std::vector<Arbor*> arbors)
{
	double dist = -1;


	for (int i = 0; i < arbors.size(); ++i)
	{
		if (arbors[i]->GetSize() > 40)
		{
			for (int j = i + 1 ; j < arbors.size(); ++j)
			{
				if (arbors[j]->GetSize() > 40)
				{
					double dist2 = arbors[i]->Distance(arbors[j], 1., 1., 0.05);
					if (dist2 > dist) dist = dist2;
				}
				
			}
		}
	}	
	return dist;
}


std::vector<Arbor*> FinalRecombination(std::vector<Arbor*> arbors, double energyMax) 
{
	if (arbors.size() == 1) return arbors;

	int id1, id2;


	id1 = arborsMaxSizeTop3(arbors)[0];
	for (int i = 0; i < arbors.size(); ++i)
	{	
		double energy = arbors[i]->GetEnergy() + arbors[id1]->GetEnergy();
		if (arbors[id1]->Distance(arbors[i], 1., 1., 0.05) < distance(arbors)/3. && energy < energyMax+sqrt(energyMax)*0.5 )
		{
			arbors[i]->SetArborId(id1);
		}
	}

	arbors = CleanArbor(arbors);
	id1 = arborsMaxSizeTop3(arbors)[0];
	id2 = arborsMaxSizeTop3(arbors)[1];
	for (int i = 0; i < arbors.size(); ++i)
	{	
		double energy = arbors[i]->GetEnergy() + arbors[id2]->GetEnergy();
		if (arbors[id2]->Distance(arbors[i], 1., 1., 0.05) < distance(arbors)/2. && energy < energyMax+sqrt(energyMax)*0.5 && arbors[i]->GetArborId() != id1)
		{
			arbors[i]->SetArborId(id2);
		}
	}

	arbors = CleanArbor(arbors);
	id1 = arborsMaxSizeTop3(arbors)[0];
	for (int i = 0; i < arbors.size(); ++i)
	{	
		double energy = arbors[i]->GetEnergy() + arbors[id1]->GetEnergy();
		if (arbors[id1]->Distance(arbors[i], 1., 1., 0.05) < distance(arbors)/1.5 && energy < energyMax+sqrt(energyMax)*0.5 )
		{
			arbors[i]->SetArborId(id1);
		}
	}

	arbors = CleanArbor(arbors);



	return arbors;
}

std::vector<Arbor*> FeedNeutral(std::vector<Arbor*> arbors) 
{
	if (arbors.size() == 1) return arbors;
	int id2 = arborsMaxSizeTop3(arbors)[1];
	int id1 = arborsMaxSizeTop3(arbors)[0];
	for (int i = 0; i < arbors.size(); ++i)
	{	
		double energy = arbors[i]->GetEnergy() + arbors[id2]->GetEnergy();

		if (i != id1)
		{
			arbors[i]->SetArborId(id2);
		}
	}

	arbors = CleanArbor(arbors);
	return arbors;
}

std::vector<Arbor*> MakeArbors(std::vector<Cluster*> clustersByLayer, double coefS, double coefT, bool reset = false) 
{
	/*
	Fonction Principale, elle combine et utilise l'enssembles des fonction défini afain de chercher a reconstruire le plus fidelement les gerbes
	*/
double energyMax = 30. + sqrt(30.)*0.5;
double energyMin = 30. - sqrt(30.)*0.5;

int id1;
int id2;

std::vector<Arbor*> arbors;

		//35. si mangé a 5 cm
		arbors = BuildPrimaryTrace(clustersByLayer, 45., 1000, 0., 20., reset);

		arbors = BuildArborSeedMissing(arbors, true, 3, 25., 10, 10000);
		arbors = FusionArbor(arbors, false, 2., 10, 1000, 1);
		arbors = ClusteringSmallArbor(arbors, 5, 25., .4);
		arbors = FusionArbor(arbors, false, 2., 10, 1000, 1);
		arbors = BuildArborSeedMissing(arbors, false, 30, 19., 10, 10000);

		for (double i = 0.5; i < 4 ; i = ++i) arbors = ClusteringSmallArbor(arbors, 20, coefS*i, .4);

		arbors = FusionArbor(arbors, false, 4., 10, 1000, 1);

		arbors = TryFusionArbor(arbors, 30, 100.);
		arbors = EnergyReconstruction(arbors, 30., distance(arbors), 10000, arborsMaxSizeTop3(arbors)[0]);
		arbors = TryFusionArbor(arbors, 30, 100);

		arbors = FinalRecombination(arbors, 30.);

		



		id1 = arborsMaxSizeTop3(arbors)[0];//Juste une fonction pour trouver les arbres les plus gros
			



//si il y a trop d'énergie ou pas asser dans la particule 1, alors on recommence avec d'autre parametre plus ou moin strict selon le besoin

		if (arbors[id1]->GetEnergy() > energyMax &&  distance(arbors) < 500.)
		{
			arbors = BuildPrimaryTrace(clustersByLayer, 35., 1000, 0., 20., reset);

			arbors = BuildArborSeedMissing(arbors, true, 3, 25., 10, 10000);
			arbors = FusionArbor(arbors, false, 1., 10, 1000, 1);
			arbors = ClusteringSmallArbor(arbors, 5, 25., .4);
			arbors = FusionArbor(arbors, false, 2., 10, 1000, 1);
			arbors = BuildArborSeedMissing(arbors, false, 30, 19., 10, 10000);

			for (double i = 0.5; i < 3 ; i = ++i) arbors = ClusteringSmallArbor(arbors, 20, coefS*i, .4);

			arbors = FusionArbor(arbors, false, 3., 10, 1000, 1);

			arbors = TryFusionArbor(arbors, 30, 40.);
			arbors = EnergyReconstruction(arbors, 30., distance(arbors), 10000, arborsMaxSizeTop3(arbors)[0]);
			
			arbors = FinalRecombination(arbors, 30.);

		}else if (arbors[id1]->GetEnergy() < energyMin && distance(arbors) > 200.)
		{
			arbors = BuildPrimaryTrace(clustersByLayer, 45., 1000, 0., 20., reset);

			arbors = BuildArborSeedMissing(arbors, true, 3, 40., 10, 10000);
			arbors = FusionArbor(arbors, false, 3., 10, 1000, 1);
			arbors = ClusteringSmallArbor(arbors, 5, 50., 1.);
			arbors = FusionArbor(arbors, false, 5., 10, 1000, 1);
			arbors = BuildArborSeedMissing(arbors, false, 30, 40., 10, 10000);


			for (double i = 0.5; i < 8 ; i = ++i) arbors = ClusteringSmallArbor(arbors, 20, coefS*i, .4);

			arbors = FusionArbor(arbors, false, 5., 10, 1000, 1);
			arbors = EnergyReconstruction(arbors, 30., distance(arbors), 10000, arborsMaxSizeTop3(arbors)[0]);
			arbors = TryFusionArbor(arbors, 30, 300.);

			arbors = FinalRecombination(arbors, 30.);
		}

		arbors = FeedNeutral(arbors);




	







	/*else
	{

		arbors = BuildPrimaryTrace(clustersByLayer, coefS, coefT, 0., 30., reset);

		arbors = BuildArborSeedMissing(arbors, true, 3, coefS-10., 10, 10000);

		arbors = ClusteringSmallArbor(arbors, 5, 25., .4);
	
		arbors = FusionArbor(arbors, false, 1., 10, 1000, 1);
		arbors = FusionArbor(arbors, false, 2., 10, 1000, 1);

		arbors = FusionArbor(arbors, false, 3., 10, 1000, 1);


		for (double i = 0.5; i < 3 ; i = ++i) arbors = ClusteringSmallArbor(arbors, 20, coefS*i, .4);

		arbors = FusionArbor(arbors, false, 4., 10, 1000, 1);

		if (arbors.size() > 2) arbors = BuildArborSeedMissing(arbors, false, 50, 25., 10, 10000);
		if (arbors.size() > 2) arbors = BuildArborSeedMissing(arbors, true, 50, 12., 10, 10000);

		arbors = TryFusionArbor(arbors);

		
	}*/


	
/*		







		arbors = FusionArbor(arbors, false, 1., 10, 1000, 1);
		arbors = FusionArbor(arbors, false, 1., 10, 1000, 1);
		arbors = FusionArbor(arbors, false, 1., 10, 1000, 1);

		arbors = BuildArborSeedMissing(arbors, true, 1, coefS-20., 10, 10000);

		//for (int i = 1; i < 10 ; ++i) arbors = ClusteringSmallArbor(arbors, 5, coefS*i, .5);

		int size = 0;
		for (int i = 0; i < arbors.size(); ++i)
		{
			if (arbors[i]->GetSize() > 10) ++size;
		}	

		double coef = 1;
		while (size > 8) 
		{
			arbors = FusionArbor(arbors, false, coef, 10, 1000, 1);

			coef +=1.;
			size = 0;
			for (int i = 0; i < arbors.size(); ++i)
			{
				if (arbors[i]->GetSize() > 10) ++size;
			}
			if (coef > 10) break;
		}

		arbors = FusionBigArbor(arbors, 10, coefS, 0.97);


		//arbors = BuildArborSeedMissing(arbors, true, 3, 25., 10, 10000);
		//arbors = ClusteringSmallArbor(arbors, 20, coefS*5., .5);

		//arbors = BuildArborSeedMissing(arbors, false, 30, 25., 10, 10000);


		
	}
*/
/*	
	

	std::vector<int> potentialParticles;
	for (int i = 0; i < arbors.size(); ++i)
	{
		if (arbors[i]->GetSize() > 40)
		{
			potentialParticles.push_back(i);
		}
	}

	if(potentialParticles.size() > 3) //suspect qu'il y ai au moin 4 particules dans la même région
	{
		arbors = FusionArbor(arbors, false, 3., 10, 1000, 1);
		arbors = TryFusionArbor(arbors);
	}


*/
	







		//arbors = CleanArbor(arbors);	//souvent inutiles, mais si il restes des hits à la fin qui n'ont pas été associer aucunes particules on les donnes a la plus petites des deux gerbes créer si jamais il y a deux gerbes







//info

/*
		for (int i = 0; i < arbors.size(); i++)
		{
			int size = arbors[i]->GetSize();
			int minLayer = arbors[i]->GetMinLayer();
			int maxLayer = arbors[i]->GetMaxLayer();
			if (size < seuil) continue;

			bool min =false;
			bool max = false;

			for (int j = 0; j < size; ++j)
			{
				Cell* cell = arbors[i]->GetCell(j);
				int* pos = cell->GetID();
				if (pos[2] == minLayer && min == false)
				{
					cout << "min : [" << pos[0] << " ; " << pos[1] << " ; " << pos[2] << "] " << size << endl;
					min = true;
				}
				if (pos[2] == maxLayer && max == false)
				{
					//cout << "max : [" << pos[0] << " ; " << pos[1] << " ; " << pos[2] << "] " << size << endl;
					max = true;
				}
			}
		}

		for (int i = 0; i < arbors.size(); ++i)
		{
			if (arbors[i]->GetSize() < seuil) continue;
			double* dir = arbors[i]->GetDirection();
			double ez[3] = {0, 0, 1};

			cout << arbors[i]->GetSize() << "	theta = " << acos(CosAngleBetweenVectors(dir, ez)) << endl;

			for (int j = i+1; j < arbors.size(); ++j)
			{
				if (arbors[j]->GetSize() < seuil) continue;
				double* dir2 = arbors[j]->GetDirection();
				cout << endl  << i << ", " << j <<"	theta = " << acos(CosAngleBetweenVectors(dir, dir2)) << endl;
			}
		}



*/


	return arbors;
}