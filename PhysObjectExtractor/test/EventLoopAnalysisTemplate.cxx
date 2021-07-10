//////////////////////////////////////////////////////////////////////
// This template analysis code has been built with fragments from the 
// classes automatically obtained by the TTree MakeClass() method.
//
// The template shows the structre of a potential analysis code
// where more TTree friends can be added with more physics objects.
//
// Done with ROOT version 5.32/00
// from TTree Events/Events
// found on file: myoutput.root
//
//
// Compile me with:
// g++ -g -O3 -Wall -Wextra -o EventLoopAnalysis EventLoopAnaalysisTemplate.cxx $(root-config --cflags --libs)
/////////////////////////////////////////////////////////////////////

//Include ROOT classes
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "TLatex.h"
//Include C++ classes
#include <iostream>
#include <vector>
	
using namespace std;

//book example histograms
TH1F* lumib = new TH1F("lumib","Luminosity block",10000,0,10000);
TH1F* h_nmu = new TH1F("h_nmu","Number of muons",200,0,200);
TH1F* h_mu_e = new TH1F("h_mu_e","Muon Energy",2000,0,2000);
TH1F* h_mu_pt = new TH1F("h_mu_pt","Muon p_{T}",2000,0,2000);


// Fixed size dimensions of array or collections stored in the TTree if any.

class EventLoopAnalysisTemplate {
public :
	
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain	
   TTree	  *tmuons;  
   //Add more trees for friendship

   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of example leaf types
   Int_t           run;
   UInt_t          luminosityBlock;
   ULong64_t	   event;
   vector<float>   *muon_e;
   vector<float>   *muon_pt;

   // List of example branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_muon_e;   //!
   TBranch        *b_muon_pt;   //!

   EventLoopAnalysisTemplate(TTree *tree=0);
   virtual ~EventLoopAnalysisTemplate();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void analysis();
};

EventLoopAnalysisTemplate::EventLoopAnalysisTemplate(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("myoutput.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("myoutput.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("myoutput.root:/myevents");
      dir->GetObject("Events",tree);

      //Get trees for friendship
      tmuons = (TTree*)f->Get("mymuons/Events");
      
      //Make friendship	
      tree->AddFriend(tmuons);
	
   }
   Init(tree);
}

EventLoopAnalysisTemplate::~EventLoopAnalysisTemplate()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventLoopAnalysisTemplate::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}



Long64_t EventLoopAnalysisTemplate::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}


void EventLoopAnalysisTemplate::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   muon_e = 0;
   muon_pt = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("muon_e", &muon_e, &b_muon_e);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   
   Notify();
}


Bool_t EventLoopAnalysisTemplate::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventLoopAnalysisTemplate::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t EventLoopAnalysisTemplate::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void EventLoopAnalysisTemplate::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L events.C
//      Root > events t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
	//Just an informative printout
	if(jentry%1000 == 0) {
	      cout<<"Processed "<<jentry<<" events out of "<<nentries<<endl;
	} 
       //Load the current event	
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       // if (Cut(ientry) < 0) continue;
       //Perform the analysis
       analysis();

    }
}


//-----------------------------------------------------------------
void EventLoopAnalysisTemplate::analysis()
{
//-----------------------------------------------------------------

   //cout<<"Make histogram of lumi blocks"<<endl;
   lumib->Fill(luminosityBlock);	

   //Loop over muons container and fill histo if cut is passed
   Int_t nmuons = muon_pt->size();
   cout<<"Number of muons = "<<nmuons<<endl;
   h_nmu->Fill(nmuons);
   float mu_pt_cut = 20; //in GeV
   if(nmuons>0){
     for (Int_t j=0; j<nmuons;++j){
       cout<<"Muon pT = "<<muon_pt->at(j)<<endl;
       if (muon_pt->at(j)>mu_pt_cut){
	 h_mu_pt->Fill(muon_pt->at(j));
	 h_mu_e->Fill(muon_e->at(j));
       }
     }
   }


}


//-----------------------------------------------------------------
int main()
{
//-----------------------------------------------------------------

    cout<<"Build the analysis object"<<endl;
    EventLoopAnalysisTemplate mytemplate;

    cout<<"Run the event loop"<<endl;
    mytemplate.Loop();

    cout<<"Save the histograms"<<endl;
    TFile* hfile = new TFile("histograms.root","RECREATE");
    lumib->Write();
    h_nmu->Write();
    h_mu_e->Write();
    h_mu_pt->Write();
    hfile->Close();

    return 1;

}


