#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TH1D.h"
#include <iostream>
#include <vector>
//#include "TLorentzVector.h"

//Compile me with:
//g++ -g -O3 -Wall -Wextra -Wpedantic -o PoetAnalysis PoetAnalysisTemplate.cxx $(root-config --cflags --libs)


void EventLoopAnalysis(){

  using namespace std;

  //book a couple of example histograms
  TH1F* h_lb = new TH1F("h_lb","Lumi Block",400,0,400);
  TH1F* h_elec_pt = new TH1F("h_elec_pt","Electron p_{T}",250,0,250);
  TH1F* h_mu_pt = new TH1F("h_mu_pt","Muon p_{T}",250,0,250);
  // ... add more if needed
  
  //open the root file to read the directories and trees of interest
  TFile* infile = new TFile("myoutput.root", "READ");
  
  //read trees of objects of interest
  TTree* tevent = (TTree*)infile->Get("myevents/Events");
  TTree* telectrons = (TTree*)infile->Get("myelectrons/Events");
  TTree* tmuons = (TTree*)infile->Get("mymuons/Events");
  // ...add other trees if needed

  //declare some variables of interest
  //general event
  Int_t run;
  UInt_t lumiblock;
  ULong64_t event;
  //electrons
  Int_t nelec;
  vector<float>*  elec_e = 0;
  vector<float>*  elec_pt = 0;
  //muons
  Int_t nmu;
  vector<float>*  mu_e = 0;
  vector<float>*  mu_pt = 0;
  // ...add variables for other objects if needed
  
  //Set branch addresses of variables of interest
  //for the different objects
  tevent->SetBranchAddress("run",&run);
  tevent->SetBranchAddress("luminosityBlock",&lumiblock);
  tevent->SetBranchAddress("event",&event);
  //electrons
  telectrons->SetBranchAddress("numberelectron",&nelec);
  telectrons->SetBranchAddress("electron_e",&elec_e);
  telectrons->SetBranchAddress("electron_pt",&elec_pt);
  //muons
  tmuons->SetBranchAddress("numbermuon",&nmu);
  tmuons->SetBranchAddress("muon_e",&mu_e);
  tmuons->SetBranchAddress("muon_pt",&mu_pt);
  //... add for other objects if needed
  
  //Add friend trees so we have access to their variables
  //while looping over events
  tevent->AddFriend(telectrons);
  tevent->AddFriend(tmuons);
  //...add other friends if needed
  
  //get number of entries in events Tree (which must be the same
  //as for the other objects trees)
  Long64_t nentries = tevent->GetEntries();
  
  //loop over the events
  for (Long64_t i=0;i<nentries;i++) {
    cout<<"Working on event # "<<i+1<<endl;
    cout<<"Number of muons: "<<nmu<<endl;
    cout<<"Number of electrons: "<<nelec<<endl;

    tevent->GetEntry(i);
    //Fill lumiblock histo
    h_lb->Fill(lumiblock);
    //loop over electrons and fill pT histogram if pT>pt_cut
    float elec_pt_cut = 5;
    for (Int_t j=0; j<nelec;++j){
      if (elec_pt->at(j)>elec_pt_cut) h_elec_pt->Fill(elec_pt->at(j));
    }
    //loop over muons and fill pT histogram if pT>pt_cut
    float mu_pt_cut = 20;
    for (Int_t j=0; j<nmu;++j){
      if (mu_pt->at(j)>mu_pt_cut) h_mu_pt->Fill(mu_pt->at(j));
    }
    
  }//----close loop over vents
  
  //open file to store histograms
  TFile* hfile = new TFile("myhistograms.root","RECREATE");
  h_lb->Write();
  h_elec_pt->Write();
  h_mu_pt->Write();
  hfile->Close();

  delete infile;
  
}//---EventLoopAnalysis()



void RDFrameAnalysis(){
  using namespace ROOT;
  using namespace std;
  //multi-threading does not seem to work with the way RDF are built below
  //EnableImplicitMT();
  
  //open the root file to read the directories and trees of interest
  TFile* infile = new TFile("file:myoutput.root", "READ");
  
  //read trees of objects of interest
  TTree* tevent = (TTree*)infile->Get("myevents/Events");
  TTree* telectrons = (TTree*)infile->Get("myelectrons/Events");
  TTree* tmuons = (TTree*)infile->Get("mymuons/Events");
  // ...add other trees if needed
  
  //Add friend trees so we have access to their variables
  //while building the RDataFrame
  tevent->AddFriend(telectrons);
  tevent->AddFriend(tmuons);
  //...add other friends if needed

  //build the RDFrame
  RDataFrame d(*tevent);


  //auto colNames = d.GetColumnNames();
  // Print columns' names
  //for (auto &&colName : colNames) std::cout << colName << std::endl;

  //Define new filtered objects and fill histograms
  //https://root.cern/doc/master/classROOT_1_1VecOps_1_1RVec.html
  //https://root.cern.ch/doc/v614/vo003__LogicalOperations_8C.html
  auto h_elec_pt = d.Define("felectrons","electron_pt[electron_pt>5]").Histo1D("felectrons");
  auto h_mu_pt = d.Define("fmuons","muon_pt[muon_pt>20]").Histo1D("fmuons");
  auto h_lb = d.Histo1D("luminosityBlock");
  
  //save histograms in a new file
  TFile* hfile = new TFile("myRDFhistograms.root","RECREATE");
  h_lb->Write();
  h_elec_pt->Write();
  h_mu_pt->Write();
  hfile->Close();
  
}//---RDFrameAnalysis()


//The main function like any other C++ code
int main (){
  using namespace std;
  cout<<"Running the EventLoopAnalysis() ...."<<endl;
  //This is the classical analysis using an event loop
  EventLoopAnalysis();
  cout<<"Running the RDFrameAnalysis() ...."<<endl;
  //This is a more sophisticated analyisis using a columnar approach
  RDFrameAnalysis();

  

}
