#define SDHCAL_PFA_cxx


#include "SDHCAL_PFA.h"
#include <TH2.h>
#include <TStyle.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <TPrincipal.h>

#include "HitGroupingByLayer.C"

//taille cellule detecteur
int xmax = 96;
int ymax = 96;
int zmax = 48;




//datas
bool print = true;
int countEvent = 0.;
double meanCluster = 0.;
double pur1 = 0.;
double pur2 = 0.;
double meanDataUsed = 0.;
double completion = 0.;
double comp1 = 0;
double comp2 = 0;


void SDHCAL_PFA::Begin(TTree * /*tree*/)
{
    TString option = GetOption();
    cout << "Process begin" << endl;
}

void SDHCAL_PFA::SlaveBegin(TTree * /*tree*/){}







Bool_t SDHCAL_PFA::Process(Long64_t entry)
{
	fReader.SetLocalEntry(entry);
	GetEntry(entry);

	std::vector<Cluster*> clustersByLayer; 
	std::vector<Arbor*> arbors;


	//datas
	int hitsLost = 0;
	int hitFromParticle1 = 0;
	int hitFromParticle2 = 0;
	int countCluster = 0;


//map du detecteur avec chaque cellule toucher ou non

	int IJK[96][96][48]={0};
	Map* map;
	map = new Map(xmax, ymax, zmax);
	for (int hit = 0; hit < K.GetSize(); ++hit)
	{
		if (particle1Tag[hit] == 1) 
		{
			IJK[I[hit]][J[hit]][K[hit]] += 1;
			++hitFromParticle1;
		}
		if (particle2Tag[hit] == 1 ) 
		{
			IJK[I[hit]][J[hit]][K[hit]] += 1;
			++hitFromParticle2;
		}
	}

	for (int k = 0; k < zmax; ++k) { 
	for (int j = 0; j < ymax; ++j) { 
	for (int i = 0; i < xmax; ++i)
	{
		int hitID = 0;
		if (IJK[i][j][k] == 0)
		{
			Cell* newCell = new Cell(i, j, k);
			map->AddCell(newCell, i, j, k);
		}else
		{
			while((I[hitID] != i && hitID < K.GetSize()) || (J[hitID] != j && hitID < K.GetSize()) || (K[hitID] != k && hitID < K.GetSize())) hitID += 1;
			
			int pID = 0;
			if (hitID < K.GetSize())
			{
				if (particle1Tag[hitID] == 1) pID = 1;
					else pID = 2;
				Cell* newCell = new Cell(i, j, k, x[hitID], y[hitID], z[hitID], pID, thr[hitID], time[hitID]);
				map->AddCell(newCell, i, j, k, true);
			}else
			{	//fix somes error of lost data
				Cell* newCell = new Cell(i, j, k);
				map->AddCell(newCell, i, j, k);
			}
		}
	}}}
	






//Fonction principal, se trouve a la toute fin de HitGroupingByLayer.C

	clustersByLayer = MakeClustersByLayer(map);
	arbors = MakeArbors(clustersByLayer, 45., 0.005);

	if (arbors.size() == 0)
	{
		cerr << "NO DATA" << endl;
		return kTRUE;
	}

	
//si trop d'énergie on tente de refaire avec des parametre plus fin (marche pas super)
	double energy;
	int reset = 0;

	for (int i = 0; i < arbors.size(); ++i)
	{
		energy = arbors[i]->GetEnergy();
		if (energy > 30. + 0.5*sqrt(30.)) ++reset;
		break;
		
	}
	if (reset/(double)arbors.size() > 0) 
	{
		//cout << " RESET: energy = " << energy << endl;
		//arbors = MakeArbors(clustersByLayer, 45., 0.0035, true);
	}









//***************************************************************
/////////////////////////////////////////////////////////////////
//***************************************************************


// datas to save and print

	++countEvent;

	int nbOfHits = 0;
	for (int i = 0; i < arbors.size(); ++i) nbOfHits += arbors[i]->GetSize();
	if (print == true)
	{
		cout << "--------------------------------------------" << endl;
		cout << "event: "<< countEvent << " -> " << endl ;
		cout << "Number of hits: " << nbOfHits << "	(" << K.GetSize() << ")" << endl;
		cout << endl << "Algorithm:" << endl;
	}

//initialisation

	int biggestArborID = 0;
	int biggestArborID2 = 0;

	double purity1 = 0.;
	double purity2 = 0.;

	double size1 = 0;
	double size2 = 0;

//identification des 2 particules si possible
	for (int i = 0; i < arbors.size(); ++i)
	{
		int size = arbors[i]->GetSize();
		if (size < 30) continue;

		int particleID = arbors[i]->GetMainParticleID();

		if (particleID != 1 && particleID != 2)
		{
			cerr << "Error: SDHCAL_PFA: Particle ID is not 1 or 2" << endl;;
			continue;
		}


		if (size > size1 && particleID == 1)
		{
			size1 = size;
			biggestArborID = i;
		}

		if (size > size2 && particleID == 2)
		{
			size2 = size;
			biggestArborID2 = i;
		}
		
	}

	if (biggestArborID == biggestArborID2)
	{
		if ( arbors.size() > 1) biggestArborID2 = 1;
	}

	int ID1, ID2;

	if (arbors[biggestArborID]->GetMainParticleID() == 1)
	{
		ID1 = 1;
		ID2 = 2;
	}else 
	{
		ID1 = 2;
		ID2 = 1;
	}

//calcul de la pureté = nombre de hit de la particule / nombre de hit 

	for (int i = 0; i < arbors[biggestArborID]->GetSize(); ++i)
	{
		if (arbors[biggestArborID]->GetCell(i)->GetParticleID() == ID1) ++purity1;
	}
	for (int i = 0; i < arbors[biggestArborID2]->GetSize(); ++i)
	{
		if (arbors[biggestArborID2]->GetCell(i)->GetParticleID() == ID2) ++purity2;
	}



//calcul de l'efficacité = nombre de hit de la particule dans la gerbe / nombre de hit total de la particule

	double completions1;
	double completions2;

	if (hitFromParticle1 != 0)  completions1 = purity1 / (double)hitFromParticle1;
		else  completions1 = 1.;
	if (hitFromParticle2 != 0)  completions2 = purity2 / (double)hitFromParticle2;
	else  completions2 = 1.;

	completion += 100.*(completions1 + completions2)/2.;
	
	purity1 = 100.*purity1 / (double)arbors[biggestArborID]->GetSize();
	purity2 = 100.*purity2 / (double)arbors[biggestArborID2]->GetSize();


	if (arbors.size() == 1)
	{
		completions2 = 0.;
	}


//affichage
	for (int i = 0; i < arbors.size(); i++)
	{
			int th1 = 0;
			int th2 = 0;
			int th3 = 0;
			for (int j = 0; j < arbors[i]->GetSize(); ++j)
			{
				if(arbors[i]->GetCell(j)->GetThreshold() == 1)
				{
					++th1;
				}
				if(arbors[i]->GetCell(j)->GetThreshold() == 2)
				{
					++th2;
				}
				if(arbors[i]->GetCell(j)->GetThreshold() == 3)
				{
					++th3;
				}
			}

		if (arbors[i]->GetSize() > seuil)
		{
			++countCluster;
			if (print == true) 
			{
				double p = 0;
				for (int c = 0; c < arbors[i]->GetSize(); ++c) 
					if (arbors[i]->GetCell(c)->GetParticleID() == 1)
						++p;

				if (i == biggestArborID) cout << "1->";
				if (i == biggestArborID2) cout << "2->";
				cout << "	n°" << countCluster << ":	" << arbors[i]->GetSize() << "	| "	 << (int)(100.*p/(double)arbors[i]->GetSize()) << "%		| thr = (" << th1 << ", " << th2 << ", " << th3 << ")" << endl;
					
				
			}
		}else hitsLost += arbors[i]->GetSize();
	}

	double dataUsed = 100.*(1 - ((hitsLost)/(double)K.GetSize()));
	meanDataUsed += (double)dataUsed;

	if (biggestArborID == biggestArborID2 && biggestArborID2 == 0)
	{
		purity2 = 0.;
		completions2 = 0.;
	}


	//if (K.GetSize() == nbOfHits)
	{
		pur1 += purity1;
		pur2 += purity2;
		comp1 += 100.*completions1;
		comp2 += 100.*completions2;
		meanCluster += (double)countCluster;
	}



	if (print == true)
	{

		cout << "	Number of clusters: " << countCluster << endl << endl;
		cout << "	hits not in cluster: " << hitsLost << endl;
		cout << "Reality:"<< endl;
		cout << "	hits from particle 1: " << hitFromParticle1 << endl;
		cout << "	hits from particle 2: " << hitFromParticle2 << endl << endl;
		cout << "Performance:" << endl;
		cout << "	data used: " <<  dataUsed << "%" << endl << endl;
		cout << "	p1 -purity: " <<  purity1 << "%" << endl <<"	   -completion: " << 100. * completions1 << "%" <<endl;
		cout << "	p2 -purity: " <<  purity2 << "%" << endl <<"	   -completion: " << 100. * completions2 << "%" <<endl << endl;
	}





//outils
	/*
	bool plot = false;

	if (plot == true && arbors.size() > 1 && countEvent == 18)
	{
	   TCanvas *c1 = new TCanvas("c1","Graph2D example",0,0,600,400);
	   Double_t P = 5.;
	   Int_t npx  = arbors[biggestArborID2]->GetSize() ;
	   Int_t npy  = arbors[biggestArborID2]->GetSize() ;
	   Double_t x;
	   Double_t y;
	   Double_t z;
	   Int_t k = 1;
	   TGraph2D *dt = new TGraph2D(npx*npy);
	   dt->SetNpy(41);
	   dt->SetNpx(40);
	   for (Int_t i=0; i<npx; i++) {

	         x = arbors[biggestArborID2]->GetCell(i)->GetPosition()[0];
	         y = arbors[biggestArborID2]->GetCell(i)->GetPosition()[1];
	         z = arbors[biggestArborID2]->GetCell(i)->GetPosition()[2];
	         dt->SetPoint(k,x,y,z);
	         ++k;

	   }
	   gStyle->SetPalette(1);
	   dt->SetMarkerStyle(1);
	   dt->Draw();
	   return c1;

}


  */















	clustersByLayer.clear();
	arbors.clear();
	
    return kTRUE;
}

void SDHCAL_PFA::SlaveTerminate(){}









void SDHCAL_PFA::Terminate()
{
	
	bool printFinal = true;

	cout << "PARAM = " << PARAM;

	if (printFinal == true)
	{

		cout << "--------------------------------------------" << endl;
		cout << "--------------- PERFORMANCE ----------------" << endl;
		cout << "--------------------------------------------" << endl << endl;

		cout << "cluster mean: 	" << meanCluster/(double)countEvent <<endl;
		cout << "data used:	" << meanDataUsed/(double)countEvent << "%" <<endl;
		cout << "purity:		" << "[" << pur1/(double)countEvent << ", " << pur2/(double)countEvent << "] " << endl;
		cout << "completion:	" << "[" << comp1/(double)countEvent << ", " << comp2/(double)countEvent << "] " << endl;

	}
	ofstream myfile;
	myfile.open ("datas.txt", ios::app);
	myfile << comp2/(double)countEvent << "	" << pur2/(double)countEvent << endl;
	myfile.close();

		
	cout << endl << "Process end: ";



















}
