//////////////////////////////////////////////////////////////////////
// This template analysis code has been built with fragments from the 
// classes automatically obtained by the TTree MakeClass() method.
//
// The template shows the structre of a potential analysis code
// where more TTree friends can be added with more physics objects.
//
// The analysis part is based on the RDataFrame analysis example by Stefan Wunch:
// https://github.com/cms-opendata-analyses/HiggsTauTauNanoAODOutreachAnalysis
//
// Done with ROOT version 5.32/00
// from TTree Events/Events
// found on file: myoutput.root
//
//
// Compile me with:
// g++ -std=c++11 -g -O3 -Wall -Wextra -o EventLoopAnalysis EventLoopAnalysisTemplate.cxx $(root-config --cflags --libs)
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
#include "TStopwatch.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/GenVector_exception.h"
#include "ROOT/RVec.hxx"

//Include C++ classes
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <math.h> 

using namespace std;

/*
 * Base path to local filesystem or to EOS containing the datasets
 */
//const std::string samplesBasePath = "root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/";
const string samplesBasePath = "";

//book example histograms for specific variables
//copy them in the constructor if you add more
const int nhists = 42;

TH1F* SMHiggsToZZTo4L_FourMuons_Higgs_mass = new TH1F("SMHiggsToZZTo4L_FourMuons_Higgs_mass","Higgs Mass",36,70,180);
TH1F* SMHiggsToZZTo4L_FourMuons_Z1_mass = new TH1F("SMHiggsToZZTo4L_FourMuons_Z1_mass","Z1 Mass",36,40,160);
TH1F* SMHiggsToZZTo4L_FourMuons_Z2_mass = new TH1F("SMHiggsToZZTo4L_FourMuons_Z2_mass","Z2 Mass",36,12,160);

TH1F* SMHiggsToZZTo4L_FourElectrons_Higgs_mass  = new TH1F("SMHiggsToZZTo4L_FourElectrons_Higgs_mass ","Higgs Mass",36,70,180);
TH1F* SMHiggsToZZTo4L_FourElectrons_Z1_mass = new TH1F("SMHiggsToZZTo4L_FourElectrons_Z1_mass","Z1 Mass",36,40,160);
TH1F* SMHiggsToZZTo4L_FourElectrons_Z2_mass = new TH1F("SMHiggsToZZTo4L_FourElectrons_Z2_mass","Z2 Mass",36,12,160);

TH1F* SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Higgs_mass  = new TH1F("SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Higgs_mass","Higgs Mass",36,70,180);
TH1F* SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Z1_mass = new TH1F("SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Z1_mass","Z1 Mass",36,40,160);
TH1F* SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Z2_mass = new TH1F("SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Z2_mass","Z2 Mass",36,12,160);

TH1F* ZZTo4mu_FourMuons_Higgs_mass = new TH1F("ZZTo4mu_FourMuons_Higgs_mass","Higgs Mass",36,70,180);
TH1F* ZZTo4mu_FourMuons_Z1_mass = new TH1F("ZZTo4mu_FourMuons_Z1_mass","Z1 Mass",36,40,160);
TH1F* ZZTo4mu_FourMuons_Z2_mass = new TH1F("ZZTo4mu_FourMuons_Z2_mass","Z2 Mass",36,12,160);

TH1F* ZZTo4e_FourElectrons_Higgs_mass  = new TH1F("ZZTo4e_FourElectrons_Higgs_mass ","Higgs Mass",36,70,180);
TH1F* ZZTo4e_FourElectrons_Z1_mass = new TH1F("ZZTo4e_FourElectrons_Z1_mass","Z1 Mass",36,40,160);
TH1F* ZZTo4e_FourElectrons_Z2_mass = new TH1F("ZZTo4e_FourElectrons_Z2_mass","Z2 Mass",36,12,160);

TH1F* ZZTo2e2mu_TwoMuonsTwoElectrons_Higgs_mass  = new TH1F("ZZTo2e2mu_TwoMuonsTwoElectrons_Higgs_mass","Higgs Mass",36,70,180);
TH1F* ZZTo2e2mu_TwoMuonsTwoElectrons_Z1_mass = new TH1F("ZZTo2e2mu_TwoMuonsTwoElectrons_Z1_mass","Z1 Mass",36,40,160);
TH1F* ZZTo2e2mu_TwoMuonsTwoElectrons_Z2_mass = new TH1F("ZZTo2e2mu_TwoMuonsTwoElectrons_Z2_mass","Z2 Mass",36,12,160);

TH1F* dataRunBMu_FourMuons_Higgs_mass = new TH1F("dataRunBMu_FourMuons_Higgs_mass","Higgs Mass",36,70,180);
TH1F* dataRunBMu_FourMuons_Z1_mass = new TH1F("dataRunBMu_FourMuons_Z1_mass","Z1 Mass",36,40,160);
TH1F* dataRunBMu_FourMuons_Z2_mass = new TH1F("dataRunBMu_FourMuons_Z2_mass","Z2 Mass",36,12,160);

TH1F* dataRunBMu_TwoMuonsTwoElectrons_Higgs_mass  = new TH1F("dataRunBMu_TwoMuonsTwoElectrons_Higgs_mass","Higgs Mass",36,70,180);
TH1F* dataRunBMu_TwoMuonsTwoElectrons_Z1_mass = new TH1F("dataRunBMu_TwoMuonsTwoElectrons_Z1_mass","Z1 Mass",36,40,160);
TH1F* dataRunBMu_TwoMuonsTwoElectrons_Z2_mass = new TH1F("dataRunBMu_TwoMuonsTwoElectrons_Z2_mass","Z2 Mass",36,12,160);

TH1F* dataRunCMu_FourMuons_Higgs_mass = new TH1F("dataRunCMu_FourMuons_Higgs_mass","Higgs Mass",36,70,180);
TH1F* dataRunCMu_FourMuons_Z1_mass = new TH1F("dataRunCMu_FourMuons_Z1_mass","Z1 Mass",36,40,160);
TH1F* dataRunCMu_FourMuons_Z2_mass = new TH1F("dataRunCMu_FourMuons_Z2_mass","Z2 Mass",36,12,160);

TH1F* dataRunCMu_TwoMuonsTwoElectrons_Higgs_mass  = new TH1F("dataRunCMu_TwoMuonsTwoElectrons_Higgs_mass","Higgs Mass",36,70,180);
TH1F* dataRunCMu_TwoMuonsTwoElectrons_Z1_mass = new TH1F("dataRunCMu_TwoMuonsTwoElectrons_Z1_mass","Z1 Mass",36,40,160);
TH1F* dataRunCMu_TwoMuonsTwoElectrons_Z2_mass = new TH1F("dataRunCMu_TwoMuonsTwoElectrons_Z2_mass","Z2 Mass",36,12,160);

TH1F* dataRunBElec_FourElectrons_Higgs_mass  = new TH1F("dataRunBElec_FourElectrons_Higgs_mass ","Higgs Mass",36,70,180);
TH1F* dataRunBElec_FourElectrons_Z1_mass = new TH1F("dataRunBElec_FourElectrons_Z1_mass","Z1 Mass",36,40,160);
TH1F* dataRunBElec_FourElectrons_Z2_mass = new TH1F("dataRunBElec_FourElectrons_Z2_mass","Z2 Mass",36,12,160);

TH1F* dataRunBElec_TwoMuonsTwoElectrons_Higgs_mass  = new TH1F("dataRunBElec_TwoMuonsTwoElectrons_Higgs_mass","Higgs Mass",36,70,180);
TH1F* dataRunBElec_TwoMuonsTwoElectrons_Z1_mass = new TH1F("dataRunBElec_TwoMuonsTwoElectrons_Z1_mass","Z1 Mass",36,40,160);
TH1F* dataRunBElec_TwoMuonsTwoElectrons_Z2_mass = new TH1F("dataRunBElec_TwoMuonsTwoElectrons_Z2_mass","Z2 Mass",36,12,160);

TH1F* dataRunCElec_FourElectrons_Higgs_mass  = new TH1F("dataRunCElec_FourElectrons_Higgs_mass ","Higgs Mass",36,70,180);
TH1F* dataRunCElec_FourElectrons_Z1_mass = new TH1F("dataRunCElec_FourElectrons_Z1_mass","Z1 Mass",36,40,160);
TH1F* dataRunCElec_FourElectrons_Z2_mass = new TH1F("dataRunCElec_FourElectrons_Z2_mass","Z2 Mass",36,12,160);

TH1F* dataRunCElec_TwoMuonsTwoElectrons_Higgs_mass  = new TH1F("dataRunCElec_TwoMuonsTwoElectrons_Higgs_mass","Higgs Mass",36,70,180);
TH1F* dataRunCElec_TwoMuonsTwoElectrons_Z1_mass = new TH1F("dataRunCElec_TwoMuonsTwoElectrons_Z1_mass","Z1 Mass",36,40,160);
TH1F* dataRunCElec_TwoMuonsTwoElectrons_Z2_mass = new TH1F("dataRunCElec_TwoMuonsTwoElectrons_Z2_mass","Z2 Mass",36,12,160);

//Requiered trigger
//const string triggerRequest = "HLT_L2DoubleMu23_NoVertex_v*";

// Fixed size dimensions of array or collections stored in the TTree if any.

class EventLoopAnalysisTemplate {
public :
	
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain	
  TTree          *ttrigger;
  TTree          *tevents;
  TTree          *tvertex;
  TTree          *tmuons;
  TTree          *ttaus;
  TTree          *tmets;
  TTree          *telectrons;
  Int_t          fCurrent; //!current Tree number in a TChain

  //for managing files and weights
  TString          labeltag;
  TString          filename;
  Float_t          theweight;

  //array to keep histograms to be written and easily loop over them
  TH1F            *hists[nhists];
  
  // Declaration of example leaf types
  Int_t           run;
  UInt_t          luminosityBlock;
  ULong64_t	   event;
  Int_t           PV_npvs;
  std::map<std::string, int> *triggermap;
  vector<float>   *muon_pt;
  vector<float>   *muon_eta;
  vector<float>   *muon_phi;
  vector<float>   *muon_ch;
  vector<float>   *muon_tightid;
  vector<float>   *muon_pfreliso03all;
  vector<float>   *muon_mass;
  vector<float>   *electron_pt;
  vector<float>   *electron_eta;
  vector<float>   *electron_phi;
  vector<float>   *electron_ch;
  vector<float>   *electron_iso;
  vector<float>   *electron_dxy;
  vector<float>   *electron_dz;
  vector<float>   *electron_dxyError;
  vector<float>   *electron_dzError;
  vector<float>   *electron_mass;
  
  // List of example branches
  TBranch        *b_run;   //!
  TBranch        *b_luminosityBlock;   //!
  TBranch        *b_event;   //!
  TBranch        *b_PV_npvs;   //!
  TBranch        *b_triggermap;   //!
  TBranch        *b_muon_pt;   //!
  TBranch        *b_muon_eta;   //!
  TBranch        *b_muon_phi;   //!
  TBranch        *b_muon_ch;   //!
  TBranch        *b_muon_tightid;   //!
  TBranch        *b_muon_pfreliso03all;   //!
  TBranch        *b_muon_mass;   //!
  TBranch        *b_electron_pt;   //!
  TBranch        *b_electron_eta;   //!
  TBranch        *b_electron_phi;   //!
  TBranch        *b_electron_ch;   //!
  TBranch        *b_electron_iso;   //!
  TBranch        *b_electron_dxy;   //!
  TBranch        *b_electron_dz;   //!
  TBranch        *b_electron_dxyError;   //!
  TBranch        *b_electron_dzError;   //!
  TBranch        *b_electron_mass;   //!
  
  
  EventLoopAnalysisTemplate(TString filename, TString labeltag, Float_t theweight);
  virtual ~EventLoopAnalysisTemplate();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  ROOT::RVec<ROOT::RVec<int>> reconstruct_samekind_electron();
  ROOT::RVec<ROOT::Math::PtEtaPhiMVector> z_fourvectors_samekind_electron(ROOT::RVec<ROOT::RVec<int>> idx);
  ROOT::RVec<ROOT::RVec<int>> reconstruct_samekind_muon();
  ROOT::RVec<ROOT::Math::PtEtaPhiMVector> z_fourvectors_samekind_muon(ROOT::RVec<ROOT::RVec<int>> idx);
  ROOT::RVec<ROOT::Math::PtEtaPhiMVector> z_fourvectors_2el2mu();
  bool filter_deltar_electron(ROOT::RVec<ROOT::RVec<int>> idx);
  bool filter_deltar_muon(ROOT::RVec<ROOT::RVec<int>> idx);
  std::string finalState();
  vector<float> ReconstructHiggs(std::string fs);
  
  void analysis();
  
};

//Constructor of the analysis class
EventLoopAnalysisTemplate::EventLoopAnalysisTemplate(TString thefile, TString thelabel, Float_t sampleweight) : fChain(0)
{
  //Prepare some info for the analysis object:
  filename = thefile;
  labeltag = thelabel;
  theweight = sampleweight;
  

  //Load histograms easy access
  //hast to be in agreement with above definitions.
  hists[0] = SMHiggsToZZTo4L_FourMuons_Higgs_mass;
  hists[1] = SMHiggsToZZTo4L_FourMuons_Z1_mass;
  hists[2] = SMHiggsToZZTo4L_FourMuons_Z2_mass;

  hists[3] = SMHiggsToZZTo4L_FourElectrons_Higgs_mass;
  hists[4] = SMHiggsToZZTo4L_FourElectrons_Z1_mass;
  hists[5] = SMHiggsToZZTo4L_FourElectrons_Z2_mass;

  hists[6] = SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Higgs_mass;
  hists[7] = SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Z1_mass;
  hists[8] = SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Z2_mass;

  hists[9] = ZZTo4mu_FourMuons_Higgs_mass;
  hists[10] = ZZTo4mu_FourMuons_Z1_mass;
  hists[11] = ZZTo4mu_FourMuons_Z2_mass;

  hists[12] = ZZTo4e_FourElectrons_Higgs_mass;
  hists[13] = ZZTo4e_FourElectrons_Z1_mass;
  hists[14] = ZZTo4e_FourElectrons_Z2_mass;

  hists[15] = ZZTo2e2mu_TwoMuonsTwoElectrons_Higgs_mass;
  hists[16] = ZZTo2e2mu_TwoMuonsTwoElectrons_Z1_mass;
  hists[17] = ZZTo2e2mu_TwoMuonsTwoElectrons_Z2_mass;

  hists[18] = dataRunBMu_FourMuons_Higgs_mass;
  hists[19] = dataRunBMu_FourMuons_Z1_mass;
  hists[20] = dataRunBMu_FourMuons_Z2_mass;

  hists[21] = dataRunBMu_TwoMuonsTwoElectrons_Higgs_mass;
  hists[22] = dataRunBMu_TwoMuonsTwoElectrons_Z1_mass;
  hists[23] = dataRunBMu_TwoMuonsTwoElectrons_Z2_mass;

  hists[24] = dataRunCMu_FourMuons_Higgs_mass;
  hists[25] = dataRunCMu_FourMuons_Z1_mass;
  hists[26] = dataRunCMu_FourMuons_Z2_mass;

  hists[27] = dataRunCMu_TwoMuonsTwoElectrons_Higgs_mass;
  hists[28] = dataRunCMu_TwoMuonsTwoElectrons_Z1_mass;
  hists[29] = dataRunCMu_TwoMuonsTwoElectrons_Z2_mass;

  hists[30] = dataRunBElec_FourElectrons_Higgs_mass;
  hists[31] = dataRunBElec_FourElectrons_Z1_mass;
  hists[32] = dataRunBElec_FourElectrons_Z2_mass;

  hists[33] = dataRunBElec_TwoMuonsTwoElectrons_Higgs_mass;
  hists[34] = dataRunBElec_TwoMuonsTwoElectrons_Z1_mass;
  hists[35] = dataRunBElec_TwoMuonsTwoElectrons_Z2_mass;

  hists[36] = dataRunCElec_FourElectrons_Higgs_mass;
  hists[37] = dataRunCElec_FourElectrons_Z1_mass;
  hists[38] = dataRunCElec_FourElectrons_Z2_mass;

  hists[39] = dataRunCElec_TwoMuonsTwoElectrons_Higgs_mass;
  hists[40] = dataRunCElec_TwoMuonsTwoElectrons_Z1_mass;
  hists[41] = dataRunCElec_TwoMuonsTwoElectrons_Z2_mass;


// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TTree* tree = 0;
   TFile *f = TFile::Open(filename);
   //trigger should go first as it is the more complicated one
   tree = (TTree*)f->Get("mytriggers/Events");
   //Get trees for friendship
   tevents = (TTree*)f->Get("myevents/Events");
   tvertex = (TTree*)f->Get("mypvertex/Events");
   tmuons = (TTree*)f->Get("mymuons/Events");
   ttaus = (TTree*)f->Get("mytaus/Events");
   tmets = (TTree*)f->Get("mymets/Events");
   ttaus = (TTree*)f->Get("myelectrons/Events");
   //Make friends so we can have access to friends variables	
   //we may not use all of the available information
   //it is just an example
   tree->AddFriend(tevents);
   tree->AddFriend(tvertex);
   tree->AddFriend(tmuons);
   tree->AddFriend(ttaus);
   tree->AddFriend(telectrons);
   tree->AddFriend(tmets);
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
  //cout<<" Set the environment to read one entry"<<endl;
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
   //In our case, only vectors and maps
   triggermap =0;
   muon_pt = 0;
   muon_eta = 0;
   muon_phi = 0;
   muon_ch = 0;
   muon_tightid = 0;
   muon_pfreliso03all = 0;
   muon_mass = 0;
   electron_pt = 0;
   electron_eta = 0;
   electron_phi = 0;
   electron_ch = 0;
   electron_iso = 0;
   electron_dxy = 0;
   electron_dz = 0;
   electron_dxyError = 0;
   electron_dzError = 0;
   electron_mass = 0;

  
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //Comment out to be able to read map
   //https://root-forum.cern.ch/t/std-map-in-ttree-with-makeclass/14171
   //fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("triggermap",&triggermap,&b_triggermap);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_ch", &muon_ch, &b_muon_ch);
   fChain->SetBranchAddress("muon_tightid", &muon_tightid, &b_muon_tightid);
   fChain->SetBranchAddress("muon_pfreliso03all", &muon_pfreliso03all, &b_muon_pfreliso03all);
   fChain->SetBranchAddress("muon_mass", &muon_mass, &b_muon_mass);
   fChain->SetBranchAddress("electron_pt", &electron_pt, &b_electron_pt);
   fChain->SetBranchAddress("electron_eta", &electron_eta, &b_electron_eta);
   fChain->SetBranchAddress("electron_phi", &electron_phi, &b_electron_phi);
   fChain->SetBranchAddress("electron_ch", &electron_ch, &b_electron_ch);
   fChain->SetBranchAddress("electron_iso", &electron_iso, &b_electron_iso);
   fChain->SetBranchAddress("electron_dxy", &electron_dxy, &b_electron_dxy);
   fChain->SetBranchAddress("electron_dz", &electron_dz, &b_electron_dz);
   fChain->SetBranchAddress("electron_dxyError", &electron_dxyError, &b_electron_dxyError);
   fChain->SetBranchAddress("electron_dzError", &electron_dzError, &b_electron_dzError);
   fChain->SetBranchAddress("electron_mass", &electron_mass, &b_electron_mass);
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

void EventLoopAnalysisTemplate::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //Just an informative printout (change for more or less printing)
    if(jentry%10 == 0) {
      cout<<"Processed "<<jentry<<" events out of "<<nentries<<endl;
    } 
    //cout<<"Load the current event"<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    analysis();
    
  }
}


//-----------------------------------------------------------------
void EventLoopAnalysisTemplate::analysis()
{
//-----------------------------------------------------------------

  //cout<<"analysis() execution"<<endl;

  //minimal selection with trigger requirement
  //if (!MinimalSelection()) return;
  
  const std::string& fs = finalState();
  vector<float> vecm = ReconstructHiggs(fs);
  
  //fill histogram
  
  if(vecm[0] !=0 && vecm[1] !=0 && vecm[2] !=0){
    Int_t histsize = sizeof(hists)/sizeof(hists[0]);
    for(Int_t j=0;j<histsize;++j){
       TString histname = TString(hists[j]->GetName());
       TString thelabel = histname(0,histname.First("_"));
       TString temp = histname(histname.First("_")+1,histname.Sizeof());
       TString tempfs = temp(0,temp.First("_"));
       TString thevar = temp(temp.First("_")+1,temp.Sizeof());
       
       if(fs.compare("FourMuons") == 0){
         if(thelabel == labeltag && tempfs == "FourMuons"){ 
         //Z1
         if(thevar == "Z1_mass"){
           hists[j]->Fill(vecm[0],theweight);
           }
         //Z2
         if(thevar == "Z2_mass"){
           hists[j]->Fill(vecm[1],theweight);
           }
         //Higgs
         if(thevar == "Higgs_mass"){
           hists[j]->Fill(vecm[2],theweight);
           }
         }   
       }
       
       if(fs.compare("FourElectrons") == 0){
         if(thelabel == labeltag && tempfs == "FourElectrons"){ 
         //Z1
         if(thevar == "Z1_mass"){
           hists[j]->Fill(vecm[0],theweight);
           }
         //Z2
         if(thevar == "Z2_mass"){
           hists[j]->Fill(vecm[1],theweight);
           }
         //Higgs
         if(thevar == "Higgs_mass"){
           hists[j]->Fill(vecm[2],theweight);
           }
         }   
       }
       
       if(fs.compare("TwoMuonsTwoElectrons") == 0){
         if(thelabel == labeltag && tempfs == "TwoMuonsTwoElectrons"){ 
         //Z1
         if(thevar == "Z1_mass"){
           hists[j]->Fill(vecm[0],theweight);
           }
         //Z2
         if(thevar == "Z2_mass"){
           hists[j]->Fill(vecm[1],theweight);
           }
         //Higgs
         if(thevar == "Higgs_mass"){
           hists[j]->Fill(vecm[2],theweight);
           }
         }   
       }    
       
    }
  
  }

}//------analysis()

//-----------------------------------------------------------------

ROOT::RVec<ROOT::RVec<int>> EventLoopAnalysisTemplate::reconstruct_samekind_electron(){
 
        ROOT::RVec<ROOT::RVec<int>> idx(2);
        idx[0].reserve(2); idx[1].reserve(2);
        const auto z_mass = 91.2;

        // Find first lepton pair with invariant mass closest to Z mass
        vector< pair<int,int> > idx_cmb ;
        Int_t n1 = electron_pt->size();
        
        for(Int_t midx=0;midx<n1;++midx){
           for(Int_t tidx=0;tidx<2;++tidx){
               idx_cmb.push_back(make_pair(midx,tidx));
              }
        }
        
        auto best_mass = -1;
        size_t best_i1 = 0; size_t best_i2 = 0;
        
        for (size_t i = 0; i < idx_cmb.size(); i++) {
            const auto i1 = idx_cmb.at(i).first;
            const auto i2 = idx_cmb.at(i).second;
            if (electron_ch->at(i1) != electron_ch->at(i2)) {
                ROOT::Math::PtEtaPhiMVector p1(electron_pt->at(i1), electron_eta->at(i1), electron_phi->at(i1), electron_mass->at(i1));
                ROOT::Math::PtEtaPhiMVector p2(electron_pt->at(i2), electron_eta->at(i2), electron_phi->at(i2), electron_mass->at(i2));
                const auto this_mass = (p1 + p2).M();
                if (std::abs(z_mass - this_mass) < std::abs(z_mass - best_mass)) {
                    best_mass = this_mass;
                    best_i1 = i1;
                    best_i2 = i2;
                }
            }
        }
        
        idx[0].emplace_back(best_i1);
        idx[0].emplace_back(best_i2);

        // Reconstruct second Z from remaining lepton pair
        for (size_t i = 0; i < 4; i++) {
            if (i != best_i1 && i != best_i2) {
                idx[1].emplace_back(i);
            }
        }

        // Return indices of the pairs building two Z bosons
        return idx;
}

//-----------------------------------------------------------------

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> EventLoopAnalysisTemplate::z_fourvectors_samekind_electron(ROOT::RVec<ROOT::RVec<int>> idx){

        ROOT::RVec<ROOT::Math::PtEtaPhiMVector> z_fourvecs(2);
        const auto z_mass = 91.2;
        
        for (size_t i = 0; i < 2; i++) {
            const auto i1 = idx[i][0];
            const auto i2 = idx[i][1];
            ROOT::Math::PtEtaPhiMVector p1(electron_pt->at(i1), electron_eta->at(i1), electron_phi->at(i1), electron_mass->at(i1));
            ROOT::Math::PtEtaPhiMVector p2(electron_pt->at(i2), electron_eta->at(i2), electron_phi->at(i2), electron_mass->at(i2));
            z_fourvecs[i] = p1 + p2;
        }
        if (std::abs(z_fourvecs[0].mass() - z_mass) < std::abs(z_fourvecs[1].mass() - z_mass)) {
            return z_fourvecs;
        } else {
            return ROOT::VecOps::Reverse(z_fourvecs);
        }
}

//-----------------------------------------------------------------

ROOT::RVec<ROOT::RVec<int>> EventLoopAnalysisTemplate::reconstruct_samekind_muon(){
 
        ROOT::RVec<ROOT::RVec<int>> idx(2);
        idx[0].reserve(2); idx[1].reserve(2);
        const auto z_mass = 91.2;

        // Find first lepton pair with invariant mass closest to Z mass
        vector< pair<int,int> > idx_cmb ;
        Int_t n1 = muon_pt->size();
        
        for(Int_t midx=0;midx<n1;++midx){
           for(Int_t tidx=0;tidx<2;++tidx){
               idx_cmb.push_back(make_pair(midx,tidx));
              }
        }
        
        auto best_mass = -1;
        size_t best_i1 = 0; size_t best_i2 = 0;
        
        for (size_t i = 0; i < idx_cmb.size(); i++) {
            const auto i1 = idx_cmb.at(i).first;
            const auto i2 = idx_cmb.at(i).second;
            if (muon_ch->at(i1) != muon_ch->at(i2)) {
                ROOT::Math::PtEtaPhiMVector p1(muon_pt->at(i1), muon_eta->at(i1), muon_phi->at(i1), muon_mass->at(i1));
                ROOT::Math::PtEtaPhiMVector p2(muon_pt->at(i2), muon_eta->at(i2), muon_phi->at(i2), muon_mass->at(i2));
                const auto this_mass = (p1 + p2).M();
                if (std::abs(z_mass - this_mass) < std::abs(z_mass - best_mass)) {
                    best_mass = this_mass;
                    best_i1 = i1;
                    best_i2 = i2;
                }
            }
        }
        
        idx[0].emplace_back(best_i1);
        idx[0].emplace_back(best_i2);

        // Reconstruct second Z from remaining lepton pair
        for (size_t i = 0; i < 4; i++) {
            if (i != best_i1 && i != best_i2) {
                idx[1].emplace_back(i);
            }
        }

        // Return indices of the pairs building two Z bosons
        return idx;
}

//-----------------------------------------------------------------

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> EventLoopAnalysisTemplate::z_fourvectors_samekind_muon(ROOT::RVec<ROOT::RVec<int>> idx){

        ROOT::RVec<ROOT::Math::PtEtaPhiMVector> z_fourvecs(2);
        const auto z_mass = 91.2;
        
        for (size_t i = 0; i < 2; i++) {
            const auto i1 = idx[i][0];
            const auto i2 = idx[i][1];
            ROOT::Math::PtEtaPhiMVector p1(muon_pt->at(i1), muon_eta->at(i1), muon_phi->at(i1), muon_mass->at(i1));
            ROOT::Math::PtEtaPhiMVector p2(muon_pt->at(i2), muon_eta->at(i2), muon_phi->at(i2), muon_mass->at(i2));
            z_fourvecs[i] = p1 + p2;
        }
        if (std::abs(z_fourvecs[0].mass() - z_mass) < std::abs(z_fourvecs[1].mass() - z_mass)) {
            return z_fourvecs;
        } else {
            return ROOT::VecOps::Reverse(z_fourvecs);
        }
}

//-----------------------------------------------------------------

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> EventLoopAnalysisTemplate::z_fourvectors_2el2mu(){

        const auto z_mass = 91.2;
                                      
        ROOT::Math::PtEtaPhiMVector p1(muon_pt->at(0), muon_eta->at(0), muon_phi->at(0), muon_mass->at(0));
        ROOT::Math::PtEtaPhiMVector p2(muon_pt->at(1), muon_eta->at(1), muon_phi->at(1), muon_mass->at(1));
        ROOT::Math::PtEtaPhiMVector p3(electron_pt->at(0), electron_eta->at(0), electron_phi->at(0), electron_mass->at(0));
        ROOT::Math::PtEtaPhiMVector p4(electron_pt->at(1), electron_eta->at(1), electron_phi->at(1), electron_mass->at(1));
        ROOT::RVec<ROOT::Math::PtEtaPhiMVector> z_fourvecs = {p1 + p2, p3 + p4};
        if (std::abs(z_fourvecs[0].mass() - z_mass) < std::abs(z_fourvecs[1].mass() - z_mass)) {
            return z_fourvecs;
        } else {
            return ROOT::VecOps::Reverse(z_fourvecs);
        }
}

//-----------------------------------------------------------------

bool EventLoopAnalysisTemplate::filter_deltar_electron(ROOT::RVec<ROOT::RVec<int>>  idx){

        for (size_t i = 0; i < 2; i++) {
            const auto i1 = idx[i][0];
            const auto i2 = idx[i][1];
            const auto dr = sqrt(pow(electron_eta->at(i1) - electron_eta->at(i2), 2) + pow(electron_phi->at(i1) - electron_phi->at(i2), 2));
            if (dr < 0.02) return false;
        }
        
        return true;
}

//-----------------------------------------------------------------

bool EventLoopAnalysisTemplate::filter_deltar_muon(ROOT::RVec<ROOT::RVec<int>>  idx){

        for (size_t i = 0; i < 2; i++) {
            const auto i1 = idx[i][0];
            const auto i2 = idx[i][1];
            const auto dr = sqrt(pow(muon_eta->at(i1) - muon_eta->at(i2), 2) + pow(muon_phi->at(i1) - muon_phi->at(i2), 2));
            if (dr < 0.02) return false;
        }
        
        return true;
}


//-----------------------------------------------------------------

std::string EventLoopAnalysisTemplate::finalState(){
     
     int nmuons = muon_pt->size();
     int nelec = electron_pt->size();
     
     if(nmuons == 4){
       return "FourMuons";
     }
     
     if(nelec == 4){
       return "FourElectrons";
     }
     
     if(nmuons == 2 && nelec == 2){
       return "TwoMuonsTwoElectrons";
     }
     
     return "";

}

//-----------------------------------------------------------------

vector<float> EventLoopAnalysisTemplate::ReconstructHiggs(std::string fs){
    
    bool fil = false;
    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> Z_fourvecs(2);
    vector<float> vecmass(3,0); //Vector of 3 components with 0 value.
    
    // Reconstruct the ZZ system for all final states
    
    if (fs.compare("FourMuons") == 0) {
        auto Z_idx = reconstruct_samekind_muon();
        fil = filter_deltar_muon(Z_idx); //Delta R separation of particles building the Z systems
        if(fil == true){
          Z_fourvecs = z_fourvectors_samekind_muon(Z_idx);
        }
        
    } else if (fs.compare("FourElectrons") == 0) {
        auto Z_idx = reconstruct_samekind_electron();
        fil = filter_deltar_electron(Z_idx); //Delta R separation of particles building the Z systems
        if(fil == true){
          Z_fourvecs = z_fourvectors_samekind_electron(Z_idx);
        }
        
    } else if (fs.compare("TwoMuonsTwoElectrons") == 0) {
        // With exactly two muons and two electrons, the reconstruction is trivial (each Z is built from two of the same kind).
        Z_fourvecs = z_fourvectors_2el2mu();
        
    } else {
        throw std::runtime_error("Unknown final state " + fs);
     
    }

    // Apply cut on the reconstructed Z masses
    if(Z_fourvecs[0].mass() > 40 && Z_fourvecs[0].mass() < 120){ //Mass of first Z candidate in [40, 120]
      if(Z_fourvecs[1].mass() > 12 && Z_fourvecs[1].mass() < 120){ //Mass of second Z candidate in [12, 120]
        // Combine the fourvectors of the two Z bosons to the fourvector of the reconstructed Higgs boson
        auto Higgs_fourvec = Z_fourvecs[0] + Z_fourvecs[1];
        vecmass[0] = Z_fourvecs[0].mass(); //1st Z Boson mass
        vecmass[1] = Z_fourvecs[1].mass(); //2nd Z Boson mass
        vecmass[2] = Higgs_fourvec.mass(); //Higgs mass
        return vecmass;
        
        }
    }
    
    return vecmass; 

}

//-----------------------------------------------------------------

int main()
{
//-----------------------------------------------------------------

  //To be able to read the trigger map
  gROOT->ProcessLine("#include<map>");

  
  //Compute event weights to be used for the respective datasets
  //const float integratedLuminosity = 4.412 * 1000.0; // Run2012B only
  //const float integratedLuminosity = 7.055 * 1000.0; // Run2012C only
  //const float integratedLuminosity = 11.467 * 1000.0; // Run2012B+C
  //const float mc_w = <xsec> / <#events> * integratedLuminosity;
  
  const float integratedLuminosity = 11.58 * 1000.0;
  const float scaleFactorZZTo4l = 1.386; // Correction of the simulation

  const float SMHiggsToZZTo4L_w = 0.0065 / 299973.0 * integratedLuminosity;
  const float ZZTo4mu_w = 0.077 / 1499064.0 * scaleFactorZZTo4l * integratedLuminosity;
  const float ZZTo4e_w = 0.077 / 1499093.0 * scaleFactorZZTo4l * integratedLuminosity;
  const float ZZTo2e2mu_w = 0.18 / 1497445.0 * scaleFactorZZTo4l * integratedLuminosity;
  const float dataRunBMu_w = 1.0;
  const float dataRunCMu_w = 1.0;
  const float dataRunBElec_w = 1.0;
  const float dataRunCElec_w = 1.0;
  
  map<string, pair<string,float> > sampleNames;
  //sampleNames.insert(make_pair("MCSample",make_pair("mc",mc_w)));
  sampleNames.insert(make_pair("SMHiggsToZZTo4L",make_pair("SMHiggsToZZTo4L",SMHiggsToZZTo4L_w)));
  sampleNames.insert(make_pair("ZZTo4mu",make_pair("ZZTo4mu",ZZTo4mu_w)));
  sampleNames.insert(make_pair("ZZTo4e",make_pair("ZZTo4e",ZZTo4e_w)));
  sampleNames.insert(make_pair("ZZTo2e2mu",make_pair("ZZTo2e2mu",ZZTo2e2mu_w)));
  sampleNames.insert(make_pair("Run2012B_DoubleMuParked",make_pair("dataRunBMu",dataRunBMu_w)));
  sampleNames.insert(make_pair("Run2012C_DoubleMuParked",make_pair("dataRunCMu",dataRunCMu_w)));
  sampleNames.insert(make_pair("Run2012B_DoubleElectron",make_pair("dataRunBElec",dataRunBElec_w)));
  sampleNames.insert(make_pair("Run2012C_DoubleElectron",make_pair("dataRunCElec",dataRunCElec_w)));
  
  //loop over sample files with names  defined above
  for(map< string,pair<string,float> >::iterator it=sampleNames.begin();
      it!=sampleNames.end();it++){
    
    TString samplename = it->first;
    TString thelabel = it->second.first;
    Float_t sampleweight = it->second.second;
    
    TStopwatch time;
    time.Start();
 
    cout << ">>> Processing sample " << samplename <<" with label "<<thelabel<<" and weight "<<sampleweight<<":" <<endl;
    
    TString filename = samplesBasePath+samplename+".root";
    
    cout<<"Build the analysis object with file "<<filename<<endl;
    EventLoopAnalysisTemplate mytemplate(filename,thelabel,sampleweight);
    
    cout<<"Run the event loop"<<endl;
    mytemplate.Loop();
    
    time.Stop();
    time.Print();
    
  }

  TFile* hfile = new TFile("histograms.root","RECREATE");

  //Save signal region histos

  SMHiggsToZZTo4L_FourMuons_Higgs_mass->Write();
  SMHiggsToZZTo4L_FourMuons_Z1_mass->Write();
  SMHiggsToZZTo4L_FourMuons_Z2_mass->Write();

  SMHiggsToZZTo4L_FourElectrons_Higgs_mass->Write();
  SMHiggsToZZTo4L_FourElectrons_Z1_mass->Write();
  SMHiggsToZZTo4L_FourElectrons_Z2_mass->Write();

  SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Higgs_mass->Write();
  SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Z1_mass->Write();
  SMHiggsToZZTo4L_TwoMuonsTwoElectrons_Z2_mass->Write();

  ZZTo4mu_FourMuons_Higgs_mass->Write();
  ZZTo4mu_FourMuons_Z1_mass->Write();
  ZZTo4mu_FourMuons_Z2_mass->Write();

  ZZTo4e_FourElectrons_Higgs_mass->Write();
  ZZTo4e_FourElectrons_Z1_mass->Write();
  ZZTo4e_FourElectrons_Z2_mass->Write();

  ZZTo2e2mu_TwoMuonsTwoElectrons_Higgs_mass->Write();
  ZZTo2e2mu_TwoMuonsTwoElectrons_Z1_mass->Write();
  ZZTo2e2mu_TwoMuonsTwoElectrons_Z2_mass->Write();

  dataRunBMu_FourMuons_Higgs_mass->Write();
  dataRunBMu_FourMuons_Z1_mass->Write();
  dataRunBMu_FourMuons_Z2_mass->Write();

  dataRunBMu_TwoMuonsTwoElectrons_Higgs_mass->Write();
  dataRunBMu_TwoMuonsTwoElectrons_Z1_mass->Write();
  dataRunBMu_TwoMuonsTwoElectrons_Z2_mass->Write();

  dataRunCMu_FourMuons_Higgs_mass->Write();
  dataRunCMu_FourMuons_Z1_mass->Write();
  dataRunCMu_FourMuons_Z2_mass->Write();

  dataRunCMu_TwoMuonsTwoElectrons_Higgs_mass->Write();
  dataRunCMu_TwoMuonsTwoElectrons_Z1_mass->Write();
  dataRunCMu_TwoMuonsTwoElectrons_Z2_mass->Write();

  dataRunBElec_FourElectrons_Higgs_mass->Write();
  dataRunBElec_FourElectrons_Z1_mass->Write();
  dataRunBElec_FourElectrons_Z2_mass->Write();

  dataRunBElec_TwoMuonsTwoElectrons_Higgs_mass->Write();
  dataRunBElec_TwoMuonsTwoElectrons_Z1_mass->Write();
  dataRunBElec_TwoMuonsTwoElectrons_Z2_mass->Write();

  dataRunCElec_FourElectrons_Higgs_mass->Write();
  dataRunCElec_FourElectrons_Z1_mass->Write();
  dataRunCElec_FourElectrons_Z2_mass->Write();

  dataRunCElec_TwoMuonsTwoElectrons_Higgs_mass->Write();
  dataRunCElec_TwoMuonsTwoElectrons_Z1_mass->Write();
  dataRunCElec_TwoMuonsTwoElectrons_Z2_mass->Write();
  
  hfile->Close();

  return 0;

}
