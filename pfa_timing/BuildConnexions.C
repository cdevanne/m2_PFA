#ifndef BUILD_CONNEXIONS_C
#define BUILD_CONNEXIONS_C

#include "BuildConnexions.h"




std::vector<Cluster*> SelectClustersForConnections(Detector* detector)
{
	std::vector<Cluster*> clusterForConnexions; 

	int nbOfCells = detector->GetHits().size();

	for (int i = 0; i < nbOfCells; ++i)
	{
		Cell* cell = detector->GetCell(detector->GetHits()[i]->GetId()[0], detector->GetHits()[i]->GetId()[1], detector->GetHits()[i]->GetId()[2]);
		if (cell->GetHited() == true && cell->GetClusterised() == false)
		{
			Cluster* cluster = new Cluster();
			cluster->AddHits(cell);
			clusterForConnexions.push_back(cluster);
		}
	}
	return clusterForConnexions;
}

void CreateFirstConnections(std::vector<Cluster*> clusters, double range, double timingBymillimeterMax, double timingBymillimeterMin)
{
	//création de connexion entre les cluster(groupe de hit sur une meme couche) selon un distance max coefSpace et une distance temporel max coefTiming
	for (int i = 0; i < clusters.size(); ++i)
	{
		clusters[i]->SetTiming();
		clusters[i]->SetPosition();
		clusters[i]->ClearConnectionBackward();
		clusters[i]->ClearConnectionForward();
	}



	for (int i = 0; i < clusters.size(); ++i)
	{

		for (int j = 0; j < clusters.size(); ++j)
		{
			int layerI = clusters[i]->GetMinLayer();
			int layerJ = clusters[j]->GetMinLayer();
			int diffLayer = layerJ - layerI;

			if (layerI < layerJ && diffLayer <= 2)
			{
				double tolerance = 11.;	//mm

				double distance = clusters[i]->Distance(clusters[j],1, 1, 1.);
				double dtMax = timingBymillimeterMax * (clusters[i]->Distance(clusters[j]) + tolerance);
				double dtMin = timingBymillimeterMin * (clusters[i]->Distance(clusters[j]) - tolerance);
				double t1 = clusters[j]->GetTiming();
				double t2 = clusters[i]->GetTiming();
				double dt = t1 - t2;			

				if (distance < range && dt < dtMax && dt > dtMin && dt !=0)
				{ 
					clusters[j]->AddConnectionBackward(i);
					clusters[i]->AddConnectionForward(j);
				}
			}
		}
	}
}


void CleanConnection(std::vector<Cluster*> clusters, double maxSpace)
{
	// garde la connexion vers l'arriere la plus probable
	for (int i = 0; i < clusters.size(); ++i)
	{
		double orientation[3] = {0, 0, 0};
		std::vector<int> bck = clusters[i]->GetConnectionBackward();

		

		if (bck.size() > 1)
		{
			std::vector<int> fwd = clusters[i]->GetConnectionForward();

			double *position = clusters[i]->GetPosition(); 

			for (int j = 0; j < bck.size(); ++j)
			{
				double *bckPosition = clusters[bck[j]]->GetPosition();
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
				double *fwdPosition = clusters[fwd[j]]->GetPosition();
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
				double *bckPosition = clusters[bck[j]]->GetPosition();
				double bckOrientation[3];



				for (int k = 0; k < 3; ++k)
				{
					bckOrientation[k] =  position[k] - bckPosition[k];
				}

				bckOrientation[0] = 0.;
				bckOrientation[1] = 0.;
				bckOrientation[2] = 1.;

				double paramAngle = acos(CosAngleBetweenVectors(orientation, bckOrientation))/M_PI;
				double paramDistance = (clusters[bck[j]]->Distance(clusters[i]))/maxSpace;

				double newParam = pow(paramDistance, 5) * pow(paramAngle, 2);

				if (newParam > param)
				{
					param = newParam;
					bckId = bck[j];

				}
			}	

			if (bckId != -1)
			{
				clusters[i]->ClearConnectionBackward();		
				clusters[i]->AddConnectionBackward(bckId);	
			}
			
		}

	}

	for (int i = 0; i <  clusters.size(); ++i)
	{
		clusters[i]->ClearConnectionForward();
	}

	for (int i = 0; i <  clusters.size(); ++i)
	{
		if (clusters[i]->GetConnectionBackward().size() > 0)
		{
			int backwardId = clusters[i]->GetConnectionBackward()[0];
			clusters[backwardId]->AddConnectionForward(backwardId);
		}
	}
}




void CreateSecondaryConnection(std::vector<Cluster*>clusters, double range, double maxTiming, double minTiming)
{
	//utilise les premieres connexions crée pour regarder si d'autre connexion peuvent être creer dans la direction de la premiere, le paterne a des similitudes avec les fonction précédente mais elle est indépendante (je doit corriger ca quand tout marchera)

	std::vector<int> connections;
	for (int i = 0; i < clusters.size(); ++i)
	{
		if (clusters[i]->GetConnectionBackward().size() > 0) 
		{
			int id = clusters[i]->GetConnectionBackward()[0];
			connections.push_back(id);
		}else
		{
			connections.push_back(-1);
		}
	}

	//clear les vecteurs


	for (int i = 0; i < clusters.size(); ++i)
	{
		clusters[i]->ClearConnectionBackward();
		clusters[i]->ClearConnectionForward();
		clusters[i]->SetSelectedBackwardConnection(connections[i]);
	}

	//creer les connection


	double maxAngle = 1.1; //rad

	for (int i = 0; i < clusters.size(); ++i)
	{
		int layerI = clusters[i]->GetMinLayer();

		bool condition = false;
		int maxDiffLayer = 3;

		// while (condition == false && maxDiffLayer < 4)
		{
			// ++maxDiffLayer;
			for (int j = 0; j < clusters.size(); ++j)
			{
				
				int layerJ = clusters[j]->GetMinLayer();
				int diffLayer = layerJ - layerI;
				int id  = clusters[i]->GetSelectedBackwardConnection();

				



				if (layerI < layerJ && diffLayer <= maxDiffLayer)
				{
					double distance = clusters[i]->Distance(clusters[j]);

					double tolerance = 11.;	//mm

			

					double dtMax = maxTiming * (clusters[i]->Distance(clusters[j]) + tolerance);
					double dtMin = minTiming * (clusters[i]->Distance(clusters[j]) - tolerance);
					double t1 = clusters[j]->GetTiming();
					double t2 = clusters[i]->GetTiming();
					double dt = t1 - t2;
					
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
						
						if (id != -1) orientation[k] = clusters[i]->GetPosition()[k] - clusters[id]->GetPosition()[k];
						fwdOrientation[k] = clusters[j]->GetPosition()[k] - clusters[i]->GetPosition()[k];
					}


					double angle = acos(CosAngleBetweenVectors(orientation, fwdOrientation));
					
					if (distance < range && dt < dtMax && dt > dtMin && angle < maxAngle)
					{ 
						clusters[j]->AddConnectionBackward(i);
						clusters[i]->AddConnectionForward(j);
					}
				}
			}
		}
	}

	// Regarde les vecteurt derrier et devant ceux de la connextion etudier et cherche dans cette zone
	//de layer en layer pour ne garder que els bonnes connection vers l'avant


	for (int layer = 47; layer > 0; --layer)
	{
		
	

		for (int i = 0; i < clusters.size(); ++i)
		{
			

			int connexionUsed = clusters[i]->GetSelectedBackwardConnection();
			int layerI = clusters[i]->GetMinLayer();

			

			if (connexionUsed != -1 && layerI == layer)
			{
				double orientation[3] = {0, 0, 0};
				double *position = clusters[i]->GetPosition(); 

				std::vector<int> bck = clusters[connexionUsed]->GetConnectionBackward();
				std::vector<int> fwd = clusters[i]->GetConnectionForward();

					

					for (int j = 0; j < bck.size(); ++j)
					{
						double *bckPosition = clusters[bck[j]]->GetPosition();
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
						double *fwdPosition = clusters[fwd[j]]->GetPosition();
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
					int diffLayerMax = 1000;

					for (int j = 0; j < bck.size(); ++j)
					{
						int diffLayer = layerI - clusters[bck[j]]->GetMinLayer();

						double *bckPosition = clusters[bck[j]]->GetPosition();
						double bckOrientation[3];

						for (int k = 0; k < 3; ++k)
						{
							bckOrientation[k] =  position[k] - bckPosition[k];
						}

						double paramAngle = acos(CosAngleBetweenVectors(orientation, bckOrientation))/M_PI;
						double paramDistance = (clusters[bck[j]]->Distance(clusters[i]))/range;
						double paramTiming = (abs(clusters[i]->GetTiming() - clusters[bck[j]]->GetTiming())/clusters[bck[j]]->Distance(clusters[i]))/maxTiming;

						double newParam = pow(paramAngle, 1) * pow(paramDistance, 2);


						if (param < newParam || diffLayer < diffLayerMax )
						{
							param = newParam;
							bckId = bck[j];
							diffLayerMax = diffLayer;
						}
					}	

					if (bckId != -1)
					{

						for (int j = 0; j < bck.size(); ++j)
						{
							if (bck[j] != bckId)
							{
								std::vector<int> tool = clusters[bck[j]]->GetConnectionForward();
								clusters[bck[j]]->ClearConnectionForward();
								for (int k = 0; k < tool.size(); ++k)
								{
									if (bck[j] != tool[k])
									{
										clusters[bck[j]]->AddConnectionForward(bck[j]);
									}
								}
							}
							
						}

						clusters[connexionUsed]->ClearConnectionBackward();		
						clusters[connexionUsed]->AddConnectionBackward(bckId);	
					}
				}
			}
		}
	}
	for (int i = 0; i < clusters.size(); ++i)
	{
		if (clusters[i]->GetConnectionBackward().size() == 0 && connections[i] != -1)
		{
			clusters[i]->AddConnectionBackward(connections[i]);	
		}
	}


}



std::vector<Arbor*> LinkCellsWithConnections(std::vector<Cluster*> clusters)
{
	//regroupe les cluster lier par des connexions entres eux au seins d'une même arbre
	std::vector<Arbor*> branches;

	for (int i = 0; i < clusters.size(); ++i)
	{
		int nbBackwardConnections = clusters[i]->GetConnectionBackward().size();

		if (nbBackwardConnections == 0)
		{
			int arborId = branches.size();
			clusters[i]->SetArborId(arborId);
			Arbor* newArbor = new Arbor();
			newArbor->AddCluster(clusters[i]);
			newArbor->SetArborId(arborId);
			branches.push_back(newArbor);
		}
	}
	for (int i = 0; i < clusters.size(); ++i)
	{
		int nbBackwardConnections = clusters[i]->GetConnectionBackward().size();

		if (nbBackwardConnections > 0)
		{
			int clusterBackwardId = clusters[i]->GetConnectionBackward()[0];
			int arborId = clusters[clusterBackwardId]->GetArborId();

			while (arborId == -1)
			{
				clusterBackwardId = clusters[clusterBackwardId]->GetConnectionBackward()[0];
				arborId = clusters[clusterBackwardId]->GetArborId();
			}

			if (arborId >= 0 && arborId < branches.size())
			{
				branches[arborId]->AddCluster(clusters[i]);
			}
		}
	}

	return branches;
}



std::vector<Arbor*> MakeBranches(Detector* detector)
{
	std::vector<Cluster*> clusters;

	double rangeMax = 40.;
	double secondRangeMax = 55.;
	double timingBymillimeterMin = .0032;
	double timingBymillimeterMax = .008;

	clusters = SelectClustersForConnections(detector);

	for (int i = 0; i < clusters.size(); ++i)
	{
		clusters[i]->SetTiming();
	}

	CreateFirstConnections(clusters, rangeMax, timingBymillimeterMax, timingBymillimeterMin);
	CleanConnection(clusters, rangeMax);
	CreateSecondaryConnection(clusters, secondRangeMax, timingBymillimeterMax, timingBymillimeterMin);

	std::vector<Arbor*> branches;

	branches = LinkCellsWithConnections(clusters);



	return branches;
}




#endif