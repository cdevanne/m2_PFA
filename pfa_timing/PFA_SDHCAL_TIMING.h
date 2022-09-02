//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Aug 27 22:50:02 2022 by ROOT version 6.24/06
// from TTree tree/tree
// found on file: rootFiles/100events20cm30pi-10kaon0.root
//////////////////////////////////////////////////////////

#ifndef PFA_SDHCAL_TIMING_h
#define PFA_SDHCAL_TIMING_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>



class PFA_SDHCAL_TIMING : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> eventNumber = {fReader, "eventNumber"};
   TTreeReaderValue<Int_t> particle1PDG = {fReader, "particle1PDG"};
   TTreeReaderValue<Int_t> particle2PDG = {fReader, "particle2PDG"};
   TTreeReaderArray<int> I = {fReader, "I"};
   TTreeReaderArray<int> J = {fReader, "J"};
   TTreeReaderArray<int> K = {fReader, "K"};
   TTreeReaderArray<float> x = {fReader, "x"};
   TTreeReaderArray<float> y = {fReader, "y"};
   TTreeReaderArray<float> z = {fReader, "z"};
   TTreeReaderArray<int> thr = {fReader, "thr"};
   TTreeReaderArray<float> time = {fReader, "time"};
   TTreeReaderArray<float> particle1Tag = {fReader, "particle1Tag"};
   TTreeReaderArray<float> particle2Tag = {fReader, "particle2Tag"};
   TTreeReaderArray<float> particle3Tag = {fReader, "particle3Tag"};


   PFA_SDHCAL_TIMING(TTree * /*tree*/ =0) { }
   virtual ~PFA_SDHCAL_TIMING() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(PFA_SDHCAL_TIMING,0);

};

#endif

#ifdef PFA_SDHCAL_TIMING_cxx
void PFA_SDHCAL_TIMING::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t PFA_SDHCAL_TIMING::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef PFA_SDHCAL_TIMING_cxx
