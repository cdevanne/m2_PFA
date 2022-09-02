#define PFA_SDHCAL_TIMING_cxx
// The class definition in PFA_SDHCAL_TIMING.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("PFA_SDHCAL_TIMING.C")
// root> T->Process("PFA_SDHCAL_TIMING.C","some options")
// root> T->Process("PFA_SDHCAL_TIMING.C+")


#include "PFA_SDHCAL_TIMING.h"
#include <TH2.h>
#include <TStyle.h>

#include "BuildPrimaryTrace.C"
#include "BuildConnexions.C"
#include "MakeArbors.C"
#include "Tools.C"
#include "ToolsMath.C"


int seuil = 10;
bool print = false;
bool plot2 = false;
int countEvent = 0.;
double meanCluster = 0.;
double pur1 = 0.;
double pur2 = 0.;
double meanDataUsed = 0.;
double completion = 0.;
double comp1 = 0;
double comp2 = 0;

int event = 0;
int eventDone = 0;


void PFA_SDHCAL_TIMING::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   cout << "Process begin" << endl;
}

void PFA_SDHCAL_TIMING::SlaveBegin(TTree * /*tree*/) //USELESS
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t PFA_SDHCAL_TIMING::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);
   GetEntry(entry);


   ++event;
   ++eventDone;

   cout << "Processing event " << eventDone << endl;

   if (eventDone < 0 || eventDone > 200) return kTRUE;

   

   //datas
   int hitsLost = 0;
   int countCluster = 0;

   //Creating the detector and saving datas


   int IJK[96][96][48]={0};

   Detector* detector;

   detector = new Detector();

   int numberOfHits =  K.GetSize();
   int hitFromParticle1 = 0;
   int hitFromParticle2 = 0;

   for (int hit = 0; hit < numberOfHits; ++hit)
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

   for (int k = 0; k < 48; ++k) { 
      for (int j = 0; j < 96; ++j) { 
         for (int i = 0; i < 96; ++i)
         {
            int hitId = 0;
            if (IJK[i][j][k] == 0)
            {
               Cell* newCell = new Cell(i, j, k);
               detector->AddCell(newCell, i, j, k);
            }else
            {
               while((I[hitId] != i && hitId < numberOfHits) || (J[hitId] != j && hitId < numberOfHits) || (K[hitId] != k && hitId < numberOfHits)) hitId += 1;
               
               int pId = 0;
               if (hitId < numberOfHits && time[hitId] < 12.)
               {
                  if (particle1Tag[hitId] == 1) pId = 1;
                     else pId = 2;
                  
                  if (pId == 1) 
                  {
                     Cell* newCell = new Cell(i, j, k, x[hitId], y[hitId], z[hitId], pId, thr[hitId], time[hitId]);
                     // Cell* newCell = new Cell(i, j, k, x[hitId]-500., y[hitId]+500., z[hitId], pId, thr[hitId], time[hitId]);
                     detector->AddCell(newCell, i, j, k, true);
                  }

                  if (pId == 2) 
                  {
                     Cell* newCell = new Cell(i, j, k, x[hitId], y[hitId], z[hitId], pId, thr[hitId], time[hitId]);
                     detector->AddCell(newCell, i, j, k, true);
                  }

               }else
               {  //fix somes errors
                  Cell* newCell = new Cell(i, j, k);
                  detector->AddCell(newCell, i, j, k);
               }
            }
         }
      }
   }

   //MAIN ALGORITHM



   std::vector<Cluster*> trace; 
   std::vector<Arbor*> branches;
   std::vector<Arbor*> arbors;

   // trace = BuildPrimaryTrace(detector);
   branches = MakeBranches(detector);
   arbors = MakeArbors(branches, trace);

   if (arbors.size() == 0)
   {
      cerr << "NO DATA" << endl;
      return kTRUE;
   }


//***************************************************************
/////////////////////////////////////////////////////////////////
//***************************************************************


// ANALYSIS



// datas to save and print

   ++countEvent;

   int nbOfHits = 0;
   for (int i = 0; i < arbors.size(); ++i) nbOfHits += arbors[i]->GetSize();

   if (print == true)
   {
      cout << "--------------------------------------------" << endl;
      cout << "event: "<< countEvent << " -> " << endl ;
      cout << "Number of hits: " << nbOfHits << "  (" << K.GetSize() << ")" << endl;
      cout << endl << "Algorithm:" << endl;
   }

//initialisation

   int biggestArborId = 0;
   int biggestArborId2 = -1;

   double purity1 = 0.;
   double purity2 = 0.;

   double size1 = 0;
   double size2 = 0;

//identification des 2 particules si possible
   for (int i = 0; i < arbors.size(); ++i)
   {
      int size = arbors[i]->GetSize();
      if (size < 30) continue;

      int particleId = arbors[i]->GetMainParticleId();

      if (particleId != 1 && particleId != 2)
      {
         cerr << "Error: SDHCAL_PFA: Particle Id is not 1 or 2" << endl;;
         continue;
      }


      if ((size > size1 && particleId == 1))
      {
         size1 = size;
         biggestArborId = i;
      }

      if ((size > size2 && particleId == 2))
      {
         size2 = size;
         biggestArborId2 = i;
      }
      
   }
   bool help = false;
   if (biggestArborId2 == -1)
   {
      if ( arbors.size() > 1)
      {
      help = true;
       biggestArborId2 = arborsMaxSizeTop3(arbors)[1];
      }
   }

   if ( arbors.size() ==1 )
   {
      biggestArborId2 = 0;
      biggestArborId = 0;
   }

   int Id1 = 1;
   int Id2 = 2;

//calcul de la pureté = nombre de hit de la particule / nombre de hit 
//calcul de l'efficacité = nombre de hit de la particule dans la gerbe / nombre de hit total de la particule


   double memoPurity2 = 0;
   double memoPurity1 = 0;

   for (int i = 0; i < arbors[biggestArborId]->GetSize(); ++i)
   {
      if (arbors[biggestArborId]->GetCell(i)->GetParticleId() == Id1) ++purity1;
      else ++memoPurity2;
   }
   for (int i = 0; i < arbors[biggestArborId2]->GetSize(); ++i)
   {
      if (arbors[biggestArborId2]->GetCell(i)->GetParticleId() == Id2) ++purity2;
      else ++memoPurity1;
   }




// for analysis, check if the id of the two particle are correctly associated   

   double helper;

   double completions1;
   double completions2;

   double memoCompletions1;
   double memoCompletions2;



   if (hitFromParticle1 != 0)  completions1 = purity1 / (double)hitFromParticle1;
      else  completions1 = 1.;
   if (hitFromParticle2 != 0)  completions2 = purity2 / (double)hitFromParticle2;
      else  completions2 = 1.;


      //helper
   if (hitFromParticle1 != 0)  memoCompletions1 = memoPurity1 / (double)hitFromParticle1;
      else  memoCompletions1 = 1.;
   if (hitFromParticle2 != 0)  memoCompletions2 = memoPurity2 / (double)hitFromParticle2;
      else  memoCompletions2 = 1.;



   
   if ((memoPurity2 > purity2 && completions2 < 0.2 && memoPurity1 > 50. && completions1 < 0.65) || (purity2 < .2 && completions1 < .2))
   {
      helper = biggestArborId;
      biggestArborId = biggestArborId2;
      biggestArborId2 = helper;

      purity1 = memoPurity1;
      purity2 = memoPurity2;     

      completions1 = memoCompletions1;
      completions2 = memoCompletions2;

      // cout << "DONE here" << endl;
   }

   

   

   purity1 = 100.*purity1 / (double)arbors[biggestArborId]->GetSize();
   purity2 = 100.*purity2 / (double)arbors[biggestArborId2]->GetSize();


   
   //cout << purity2 << " "  << 100.-purity1 << " ----- " << completions1 << " " << completions2 << endl;

   

   if (NumberOfArborInSize(arbors, seuil, 10000) == 1 && hitFromParticle2 > 0)
   {
      purity2 = 0;
      completions2 = 0;
   }

   if (completions2 < 0.05 && purity2 > 70.)
   {
      purity2 = 0.;
      // bad datas
   }


   completion += 100.*(completions1 + completions2)/2.;



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

      if (arbors[i]->GetSize() >= seuil)
      {
         ++countCluster;
         if (print == true) 
         {
            double p = 0;
            for (int c = 0; c < arbors[i]->GetSize(); ++c) 
               if (arbors[i]->GetCell(c)->GetParticleId() == 1)
                  ++p;

            if (i == biggestArborId) cout << "1->";
            else if (i == biggestArborId2) cout << "2->";
            else cout << "  ";
            cout << "   n°" << countCluster << ":  " << arbors[i]->GetSize() << "   | "    << (int)(100.*p/(double)arbors[i]->GetSize()) << "%     | thr = (" << th1 << ", " << th2 << ", " << th3 << ")" << endl;
               
            
         }
      }else hitsLost += arbors[i]->GetSize();
   }

   double dataUsed = 100.*(1 - ((hitsLost)/(double)(hitFromParticle1 + hitFromParticle2)));
   meanDataUsed += (double)dataUsed;

   if (biggestArborId == biggestArborId2 && biggestArborId2 == 0)
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

      cout << "   Number of clusters: " << countCluster << endl << endl;
      cout << "   hits not in cluster: " << hitsLost << endl;
      cout << "Reality:"<< endl;
      cout << "   hits from particle 1: " << hitFromParticle1 << endl;
      cout << "   hits from particle 2: " << hitFromParticle2 << endl << endl;
      cout << "Performance:" << endl;
      cout << "   data used: " <<  dataUsed << "%" << endl << endl;
      cout << "   p1 -purity: " <<  purity1 << "%" << endl <<"    -completion: " << 100. * completions1 << "%" <<endl;
      cout << "   p2 -purity: " <<  purity2 << "%" << endl <<"    -completion: " << 100. * completions2 << "%" <<endl << endl;
   }







//outils


   
   bool plot = false;

   if (plot2 == true)
   {
      int sizeMin = 11;
      double timeMin = 200.;
      for (int i = 0; i < time.GetSize(); ++i)
      {
         if (time[i] < timeMin) timeMin = time[i];
      }

      // TCanvas *c1 = new TCanvas("c1","Graph2D example",0,0,2400,800);
      TCanvas *c1 = new TCanvas("c1","Graph2D example",0,0,1000,500);

      c1->Divide(3,1);
      
      for (double tt = timeMin+10.; tt < timeMin + 10.1; tt = tt+.2)
      {


         TNtuple *ntuple = new TNtuple("ntuple","timing","x:y:z:time");
         TNtuple *ntuple2 = new TNtuple("ntuple","timing","x:y:z:particle1Tag");
         TNtuple *ntuple3 = new TNtuple("ntuple","timing","x:y:z:particle");

         for (int i = 0; i < x.GetSize(); ++i)
         {
            if (time[i] < tt && time[i] > 0.)
            {
               ntuple->Fill(x[i], y[i], z[i], time[i]);
               if (particle1Tag[i] == 1)
               {
                  ntuple2->Fill(x[i], y[i], z[i], 1);
               }else
               {
                  ntuple2->Fill(x[i], y[i], z[i], 0);
               }
            }
         }



         int arborId = 1;
         for (int i = 0; i < arbors.size(); ++i)
         {
            
            int size = arbors[i]->GetSize();
            if (size > sizeMin) ++arborId;
            for (int j = 0; j < size; ++j)
            {
               Cell* cell = arbors[i]->GetCell(j);
               double* X = cell->GetPosition();
               if (size > sizeMin)
               {
                  ntuple3->Fill(X[0], X[1], X[2], arborId);
               }else
               {
                  ntuple3->Fill(X[0], X[1], X[2], 0);
               }
            }
         }


         //scale
         ntuple->Fill(0, 0, 0, timeMin);
         ntuple2->Fill(0, 0, 0, 0);
         ntuple3->Fill(0, 0, 0, 0);

         ntuple->Fill(700, 700, 1200, timeMin);
         ntuple2->Fill(700, 700, 1200, 0);
         ntuple3->Fill(700, 700, 1200, 0);




         
         c1->cd(1);

         ntuple->SetMarkerStyle(21);
         ntuple->Draw("z:y:x:time","","COLZ");

         c1->cd(2);

         ntuple2->SetMarkerStyle(21);
         ntuple2->Draw("z:y:x:particle1Tag","","COLZ");

         c1->cd(3);

         ntuple3->SetMarkerStyle(21);
         ntuple3->Draw("z:y:x:particle","","COLZ");
         
         c1->Update();
         std::this_thread::sleep_for(std::chrono::milliseconds(4000));
         // int a;
         // cin >> a ;

      }

   c1->SaveAs("plot/file.root");

   // char str[128];
   // strcpy(str , "plot/event_XXXX.root" );

   int ints[] = {(event/1000)%10, (event/100)%10, (event/10)%10, (event)%10};

   char th, h, te, u;

   for (int i = 0; i < 10; ++i)
   {
      if (ints[0] == i) th = (ints[0]) + '0';
      if (ints[1] == i) h = (ints[1]) + '0';
      if (ints[2] == i) te = (ints[2]) + '0';
      if (ints[3] == i) u = (ints[3]) + '0';
   }

   cout << th << h << te << u << endl;

   char str[] = {'p','l','o','t','/','e','v','e','n','t','_', th, h, te, u, '.','r','o','o','t','\0'};


   gSystem->Rename("plot/file.root", str);


   }



   arbors.clear();
   trace.clear();
   branches.clear();


   return kTRUE;
}

void PFA_SDHCAL_TIMING::SlaveTerminate()  //USELESS
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void PFA_SDHCAL_TIMING::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.



   bool printFinal = true;
   cout << "PARAM = " << 0.;

   if (printFinal == true)
   {

      cout << "--------------------------------------------" << endl;
      cout << "--------------- PERFORMANCE ----------------" << endl;
      cout << "--------------------------------------------" << endl << endl;

      cout << "cluster mean:  " << meanCluster/(double)countEvent <<endl;
      cout << "data used:  " << meanDataUsed/(double)countEvent << "%" <<endl;
      cout << "purity:     " << "[" << pur1/(double)countEvent << ", " << pur2/(double)countEvent << "] " << endl;
      cout << "completion: " << "[" << comp1/(double)countEvent << ", " << comp2/(double)countEvent << "] " << endl;

   }
   ofstream myfile;
   myfile.open ("datas.txt", ios::app);
   myfile << comp2/(double)countEvent << "   " << pur2/(double)countEvent << endl;
   myfile.close();

      
   cout << endl << "Process end: ";

}