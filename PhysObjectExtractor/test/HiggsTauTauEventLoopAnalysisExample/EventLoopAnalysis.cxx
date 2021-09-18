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
// g++ -std=c++11 -g -O3 -Wall -Wextra -o EventLoopAnalysis EventLoopAnalysis.cxx $(root-config --cflags --libs) -lGenVector
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
//Include C++ classes
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>

using namespace std;


/*
 * Base path to local filesystem or to EOS containing the datasets
 */
const std::string samplesBasePath = "root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/";
//const std::string samplesBasePath = "skim5/";


//book example histograms for specific variables
//copy them in the constructor if you add more
const int nhists = 27;

//Histograms for signal region
TH1F* dataRunB_npv = new TH1F("dataRunB_npv","Number of primary vertices",25,5,30);
TH1F* dataRunB_m_vis = new TH1F("dataRunB_m_vis","Visible di-tau mass / GeV",30, 20, 140);
TH1F* dataRunB_eta_2 = new TH1F("dataRunB_eta_2","Tau #eta",30, -2.3, 2.3);


TH1F* dataRunC_npv  = new TH1F("dataRunC_npv","Number of primary vertices",25,5,30);
TH1F* dataRunC_m_vis = new TH1F("dataRunC_m_vis","Visible di-tau mass / GeV",30, 20, 140);
TH1F* dataRunC_eta_2 = new TH1F("dataRunC_eta_2","Tau #eta",30, -2.3, 2.3);


TH1F* ZLL_npv = new TH1F("ZLL_npv","Number of primary vertices",25,5,30);
TH1F* ZLL_m_vis= new TH1F("ZLL_m_vis","Visible di-tau mass / GeV",30, 20, 140);
TH1F* ZLL_eta_2= new TH1F("ZLL_eta_2","Tau #eta",30, -2.3, 2.3);

TH1F* TT_npv = new TH1F("TT_npv","Number of primary vertices",25,5,30);
TH1F* TT_m_vis= new TH1F("TT_m_vis","Visible di-tau mass / GeV",30, 20, 140);
TH1F* TT_eta_2= new TH1F("TT_eta_2","Tau #eta",30, -2.3, 2.3);

TH1F* W3J_npv = new TH1F("W3J_npv","Number of primary vertices",25,5,30);
TH1F* W3J_m_vis= new TH1F("W3J_m_vis","Visible di-tau mass / GeV",30, 20, 140);
TH1F* W3J_eta_2= new TH1F("W3J_eta_2","Tau #eta",30, -2.3, 2.3);

TH1F* W2J_npv = new TH1F("W2J_npv","Number of primary vertices",25,5,30);
TH1F* W2J_m_vis= new TH1F("W2J_m_vis","Visible di-tau mass / GeV",30, 20, 140);
TH1F* W2J_eta_2= new TH1F("W2J_eta_2","Tau #eta",30, -2.3, 2.3);

TH1F* W1J_npv = new TH1F("W1J_npv","Number of primary vertices",25,5,30);
TH1F* W1J_m_vis= new TH1F("W1J_m_vis","Visible di-tau mass / GeV",30, 20, 140);
TH1F* W1J_eta_2= new TH1F("W1J_eta_2","Tau #eta",30, -2.3, 2.3);

TH1F* qqH_npv = new TH1F("qqH_npv","Number of primary vertices",25,5,30);
TH1F* qqH_m_vis= new TH1F("qqH_m_vis","Visible di-tau mass / GeV",30, 20, 140);
TH1F* qqH_eta_2= new TH1F("qqH_eta_2","Tau #eta",30, -2.3, 2.3);

TH1F* ggH_npv = new TH1F("ggH_npv","Number of primary vertices",25,5,30);
TH1F* ggH_m_vis= new TH1F("ggH_m_vis","Visible di-tau mass / GeV",30, 20, 140);
TH1F* ggH_eta_2= new TH1F("ggH_eta_2","Tau #eta",30, -2.3, 2.3);


//Histograms for control region
TH1F* dataRunB_npv_cr = new TH1F("dataRunB_npv_cr","Number of primary vertices",25,5,30);
TH1F* dataRunB_m_vis_cr = new TH1F("dataRunB_m_vis_cr","Visible di-tau mass / GeV",30, 20, 140);
TH1F* dataRunB_eta_2_cr = new TH1F("dataRunB_eta_2_cr","Tau #eta",30, -2.3, 2.3);

TH1F* dataRunC_npv_cr = new TH1F("dataRunC_npv_cr","Number of primary vertices",25,5,30);
TH1F* dataRunC_m_vis_cr = new TH1F("dataRunC_m_vis_cr","Visible di-tau mass / GeV",30, 20, 140);
TH1F* dataRunC_eta_2_cr = new TH1F("dataRunC_eta_2_cr","Tau #eta",30, -2.3, 2.3);

TH1F* ZLL_npv_cr = new TH1F("ZLL_npv_cr","Number of primary vertices",25,5,30);
TH1F* ZLL_m_vis_cr = new TH1F("ZLL_m_vis_cr","Visible di-tau mass / GeV",30, 20, 140);
TH1F* ZLL_eta_2_cr = new TH1F("ZLL_eta_2_cr","Tau #eta",30, -2.3, 2.3);

TH1F* TT_npv_cr = new TH1F("TT_npv_cr","Number of primary vertices",25,5,30);
TH1F* TT_m_vis_cr = new TH1F("TT_m_vis_cr","Visible di-tau mass / GeV",30, 20, 140);
TH1F* TT_eta_2_cr = new TH1F("TT_eta_2_cr","Tau #eta",30, -2.3, 2.3);

TH1F* W3J_npv_cr = new TH1F("W3J_npv_cr","Number of primary vertices",25,5,30);
TH1F* W3J_m_vis_cr = new TH1F("W3J_m_vis_cr","Visible di-tau mass / GeV",30, 20, 140);
TH1F* W3J_eta_2_cr = new TH1F("W3J_eta_2_cr","Tau #eta",30, -2.3, 2.3);

TH1F* W2J_npv_cr = new TH1F("W2J_npv_cr","Number of primary vertices",25,5,30);
TH1F* W2J_m_vis_cr = new TH1F("W2J_m_vis_cr","Visible di-tau mass / GeV",30, 20, 140);
TH1F* W2J_eta_2_cr = new TH1F("W2J_eta_2_cr","Tau #eta",30, -2.3, 2.3);

TH1F* W1J_npv_cr = new TH1F("W1J_npv_cr","Number of primary vertices",25,5,30);
TH1F* W1J_m_vis_cr = new TH1F("W1J_m_vis_cr","Visible di-tau mass / GeV",30, 20, 140);
TH1F* W1J_eta_2_cr = new TH1F("W1J_eta_2_cr","Tau #eta",30, -2.3, 2.3);

TH1F* qqH_npv_cr = new TH1F("qqH_npv_cr","Number of primary vertices",25,5,30);
TH1F* qqH_m_vis_cr = new TH1F("qqH_m_vis_cr","Visible di-tau mass / GeV",30, 20, 140);
TH1F* qqH_eta_2_cr = new TH1F("qqH_eta_2_cr","Tau #eta",30, -2.3, 2.3);

TH1F* ggH_npv_cr = new TH1F("ggH_npv_cr","Number of primary vertices",25,5,30);
TH1F* ggH_m_vis_cr = new TH1F("ggH_m_vis_cr","Visible di-tau mass / GeV",30, 20, 140);
TH1F* ggH_eta_2_cr = new TH1F("ggH_eta_2_cr","Tau #eta",30, -2.3, 2.3);





//Requiered trigger
string triggerRequest = "HLT_IsoMu17_eta2p1_LooseIsoPFTau20";

// Fixed size dimensions of array or collections stored in the TTree if any.


class EventLoopAnalysisTemplate {
public :

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  TTree          *tevents;
  TTree          *tvertex;
  TTree          *ttrigger;
  TTree           *tmuons;
  TTree           *ttaus;
  TTree           *tmets;
  Int_t           fCurrent; //!current Tree number in a TChain

  TString          labeltag;
  TString         filename;
  Float_t          theweight;

  //array to keep histograms to be written and easily loop over them
  TH1F            *hists[nhists];
  TH1F            *hists_cr[nhists];

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
  vector<float>   *tau_pt;
  vector<float>   *tau_eta;
  vector<float>   *tau_phi;
  vector<float>   *tau_ch;
  vector<float>   *tau_iddecaymode;
  vector<float>   *tau_idisotight;
  vector<float>   *tau_idantieletight;
  vector<float>   *tau_idantimutight;
  vector<float>   *tau_reliso_all;
  vector<float>   *tau_mass;
  Float_t         met_pt;
  Float_t         met_phi;

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
  TBranch        *b_tau_pt;   //!
  TBranch        *b_tau_eta;   //!
  TBranch        *b_tau_phi;   //!
  TBranch        *b_tau_ch;   //!
  TBranch        *b_tau_iddecaymode;   //!
  TBranch        *b_tau_idisotight;   //!
  TBranch        *b_tau_idantieletight;   //!
  TBranch        *b_tau_idantimutight;   //!
  TBranch        *b_tau_reliso_all;   //!
  TBranch        *b_tau_mass;   //!
  TBranch        *b_met_pt;   //!
  TBranch        *b_met_phi;   //!

  EventLoopAnalysisTemplate(TString filename, TString labeltag, Float_t theweight);
  virtual ~EventLoopAnalysisTemplate();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  void analysis();
  bool MinimalSelection();
  bool isGoodMuon(Int_t idx);
  bool isGoodTau(Int_t idx);
  std::vector<int> FindMuonTauPair();
  float compute_mt(float pt_1, float phi_1,float pt_met, float phi_met);

};

/*
 * Helper function to compute the difference in the azimuth coordinate taking
 * the boundary conditions at 2 * pi into account.
 */
//-------------------------------------------------------------------------
namespace Helper {
  template <typename T>
  float DeltaPhi(T v1, T v2, const T c = M_PI)
  {
//-------------------------------------------------------------------------
    auto r = std::fmod(v2 - v1, 2.0 * c);
    if (r < -c) {
      r += 2.0 * c;
    }
    else if (r > c) {
      r -= 2.0 * c;
    }
    return r;
  }
}//-----------------Helper



//transverse mass computation
//-----------------------------------------------------------------
float EventLoopAnalysisTemplate::compute_mt(float pt_1, float phi_1,
					    float pt_met, float phi_met)
{
//-----------------------------------------------------------------

  const auto dphi = Helper::DeltaPhi(phi_1, phi_met);
  return sqrt(2.0 * pt_1 * pt_met * (1.0 - cos(dphi)));

}//-----compute_mt



EventLoopAnalysisTemplate::EventLoopAnalysisTemplate(TString thefile, TString thelabel, Float_t sampleweight) : fChain(0)
{
  //Prepare some info for the object:
  filename = thefile;
  labeltag = thelabel;
  theweight = sampleweight;


  //Load histograms for signal region
  hists[0] = dataRunB_npv;
  hists[1] = dataRunB_m_vis;
  hists[2] = dataRunB_eta_2;

  hists[3] = dataRunC_npv;
  hists[4] = dataRunC_m_vis;
  hists[5] = dataRunC_eta_2;

  hists[6] = ZLL_npv;
  hists[7] = ZLL_m_vis;
  hists[8] = ZLL_eta_2;

  hists[9] = TT_npv;
  hists[10] = TT_m_vis;
  hists[11] = TT_eta_2;

  hists[12] = W3J_npv;
  hists[13] = W3J_m_vis;
  hists[14] = W3J_eta_2;

  hists[15] = W2J_npv;
  hists[16] = W2J_m_vis;
  hists[17] = W2J_eta_2;

  hists[18] = W1J_npv;
  hists[19] = W1J_m_vis;
  hists[20] = W1J_eta_2;

  hists[21] = qqH_npv;
  hists[22] = qqH_m_vis;
  hists[23] = qqH_eta_2;

  hists[24] = ggH_npv;
  hists[25] = ggH_m_vis;
  hists[26] = ggH_eta_2;

  //Load histograms for control region
  hists_cr[0] = dataRunB_npv_cr;
  hists_cr[1] = dataRunB_m_vis_cr;
  hists_cr[2] = dataRunB_eta_2_cr;

  hists_cr[3] = dataRunC_npv_cr;
  hists_cr[4] = dataRunC_m_vis_cr;
  hists_cr[5] = dataRunC_eta_2_cr;

  hists_cr[6] = ZLL_npv_cr;
  hists_cr[7] = ZLL_m_vis_cr;
  hists_cr[8] = ZLL_eta_2_cr;

  hists_cr[9] = TT_npv_cr;
  hists_cr[10] = TT_m_vis_cr;
  hists_cr[11] = TT_eta_2_cr;

  hists_cr[12] = W3J_npv_cr;
  hists_cr[13] = W3J_m_vis_cr;
  hists_cr[14] = W3J_eta_2_cr;

  hists_cr[15] = W2J_npv_cr;
  hists_cr[16] = W2J_m_vis_cr;
  hists_cr[17] = W2J_eta_2_cr;

  hists_cr[18] = W1J_npv_cr;
  hists_cr[19] = W1J_m_vis_cr;
  hists_cr[20] = W1J_eta_2_cr;

  hists_cr[21] = qqH_npv_cr;
  hists_cr[22] = qqH_m_vis_cr;
  hists_cr[23] = qqH_eta_2_cr;

  hists_cr[24] = ggH_npv_cr;
  hists_cr[25] = ggH_m_vis_cr;
  hists_cr[26] = ggH_eta_2_cr;


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
  //Make friends so we can have access to friends variables
  //we may not use all of the available information
  //it is just an example
  tree->AddFriend(tevents);
  tree->AddFriend(tvertex);
  tree->AddFriend(tmuons);
  tree->AddFriend(ttaus);
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
   triggermap =0;
   muon_pt = 0;
   muon_eta = 0;
   muon_phi = 0;
   muon_ch = 0;
   muon_tightid = 0;
   muon_pfreliso03all = 0;
   muon_mass = 0;
   tau_pt = 0;
   tau_eta = 0;
   tau_phi = 0;
   tau_ch = 0;
   tau_iddecaymode = 0;
   tau_idisotight = 0;
   tau_idantieletight = 0;
   tau_idantimutight = 0;
   tau_reliso_all = 0;
   tau_mass = 0;

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
   fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_eta", &tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_phi", &tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("tau_ch", &tau_ch, &b_tau_ch);
   fChain->SetBranchAddress("tau_iddecaymode", &tau_iddecaymode, &b_tau_iddecaymode);
   fChain->SetBranchAddress("tau_idisotight", &tau_idisotight, &b_tau_idisotight);
   fChain->SetBranchAddress("tau_idantieletight", &tau_idantieletight, &b_tau_idantieletight);
   fChain->SetBranchAddress("tau_idantimutight", &tau_idantimutight, &b_tau_idantimutight);
   fChain->SetBranchAddress("tau_reliso_all", &tau_reliso_all, &b_tau_reliso_all);
   fChain->SetBranchAddress("tau_mass", &tau_mass, &b_tau_mass);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
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
	//Just an informative printout
	if(jentry%1000 == 0) {
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

  //minimal selection including trigger requirement
  if (!MinimalSelection()) return;

  //Find the best muon-tau pair and get indexes (1 is muon, 2 is tau)
  vector<int> GoodMuonTauPair = FindMuonTauPair();
  int idx_1 = GoodMuonTauPair[0];
  int idx_2 = GoodMuonTauPair[1];
  if (!(idx_1!=-1 && idx_2!=-1)) return;

  //Muon transverse mass cut for W+jets suppression
  if (!(compute_mt(muon_pt->at(idx_1),muon_phi->at(idx_1),met_pt,met_phi)<30)) return;

  //Require isolated muon for signal region
  if (!(muon_pfreliso03all->at(idx_1)<0.1)) return;


  //fill histograms for control region
  if(muon_ch->at(idx_1)*tau_ch->at(idx_2)>0){
      Int_t hists_crsize = sizeof(hists_cr)/sizeof(hists_cr[0]);
      for (Int_t j=0;j<hists_crsize;++j){

	TString histname = TString(hists_cr[j]->GetName());
	TString thelabel = histname(0,histname.First("_"));
	TString thevar = histname(histname.First("_")+1,histname.Sizeof());

	if (thelabel == labeltag){
	  //primary vertices
	  if(thevar == "npv_cr"){hists_cr[j]->Fill(PV_npvs,theweight);}
	  //eta of taus
	  if(thevar == "eta_2_cr"){hists_cr[j]->Fill(tau_eta->at(idx_2),theweight);}
	  //visible mass
	  if(thevar == "m_vis_cr"){
	    ROOT::Math::PtEtaPhiMVector p4_1(muon_pt->at(idx_1),muon_eta->at(idx_1),
					     muon_phi->at(idx_1), muon_mass->at(idx_1));
	    ROOT::Math::PtEtaPhiMVector p4_2(tau_pt->at(idx_2),tau_eta->at(idx_2),
					     tau_phi->at(idx_2), tau_mass->at(idx_2));
	    hists_cr[j]->Fill(float((p4_1+p4_2).M()),theweight);
	  }
	}
      }
  }

  //fill histograms for signal region
  if(muon_ch->at(idx_1)*tau_ch->at(idx_2)<0){
    Int_t histsize = sizeof(hists)/sizeof(hists[0]);
    for (Int_t j=0;j<histsize;++j){

      TString histname = TString(hists[j]->GetName());
      TString thelabel = histname(0,histname.First("_"));
      TString thevar = histname(histname.First("_")+1,histname.Sizeof());

      if (thelabel == labeltag){
	//primary vertices
	if(thevar == "npv"){hists[j]->Fill(PV_npvs,theweight);}
	//eta of taus
	if(thevar == "eta_2"){hists[j]->Fill(tau_eta->at(idx_2),theweight);}
	//visible mass
	if(thevar == "m_vis"){
	  ROOT::Math::PtEtaPhiMVector p4_1(muon_pt->at(idx_1),muon_eta->at(idx_1),
					   muon_phi->at(idx_1), muon_mass->at(idx_1));
	  ROOT::Math::PtEtaPhiMVector p4_2(tau_pt->at(idx_2),tau_eta->at(idx_2),
					   tau_phi->at(idx_2), tau_mass->at(idx_2));
	  hists[j]->Fill(float((p4_1+p4_2).M()),theweight);
	}
      }
    }
  }


}//------analysis()

/*
 * Perform a selection on the minimal requirements of an event
 */
//-----------------------------------------------------------------
bool EventLoopAnalysisTemplate::MinimalSelection()
{
//-----------------------------------------------------------------

 //cout<<"Applying minimal selection"<<endl;
  bool isTrigger = false;

  //Check trigger and acceptance bit
  for (map<string, int>::iterator it=triggermap->begin();it!=triggermap->end();it++){
    if(it->first.find(triggerRequest)!=string::npos &&
       it->second!=0){
	 //cout<<it->first<<"  "<<it->second<<endl;
      isTrigger = true;
    }
  }


  return isTrigger;
}//------MinimalSelection


// Give index of muon and check if passes the good muon check
//-----------------------------------------------------------------
bool EventLoopAnalysisTemplate::isGoodMuon(Int_t idx)
{
//-----------------------------------------------------------------
  bool isGoodMuon = false;

  float mu_eta_cut = 2.1;
  float mu_pt_cut = 17; //in GeV
  if (abs(muon_eta->at(idx))<mu_eta_cut &&
	muon_pt->at(idx)>mu_pt_cut &&
      bool(muon_tightid->at(idx)) ){
    isGoodMuon = true;
  }

  return isGoodMuon;

}//----------isGoodMuon





/*
 * Find the interesting taus in the tau collection
 *
 * The tau candidates in this collection represent hadronic decays of taus, which
 * means that the tau decays to combinations of pions and neutrinos in the final
 * state.
 */
//-----------------------------------------------------------------
bool EventLoopAnalysisTemplate::isGoodTau(Int_t idx)
{
//-----------------------------------------------------------------
 bool isGoodTau = false;

 float tau_eta_cut = 2.3;
 float tau_pt_cut = 20; //in GeV

 if (tau_ch->at(idx)!=0 &&
	abs(tau_eta->at(idx))<tau_eta_cut &&
	tau_pt->at(idx)>tau_pt_cut &&
	bool(tau_iddecaymode->at(idx)) &&
	bool(tau_idisotight->at(idx)) &&
	bool(tau_idantieletight->at(idx)) &&
	bool(tau_idantimutight->at(idx))){
      isGoodTau = true;
    }

 return isGoodTau;

}//-----------isGoodTau





//-----------------------------------------------------------------
std::vector<int> EventLoopAnalysisTemplate::FindMuonTauPair()
{
//-----------------------------------------------------------------

  //Find all possible pairs of muons and taus
  vector< pair<int,int> > comb;
  Int_t nmuons = muon_pt->size();
  Int_t ntaus = tau_pt->size();
  for(Int_t midx=0;midx<nmuons;++midx){
    for(Int_t tidx=0;tidx<ntaus;++tidx){
      comb.push_back(make_pair(midx,tidx));
    }
  }
  const size_t numComb= comb.size();

  //Find valid pairs based on delta r
  vector<int> validPair(numComb, 0);
  for(size_t i = 0; i < numComb; i++) {
    const int i1 = comb.at(i).first;
    const int i2 = comb.at(i).second;
    if(isGoodMuon(i1) && isGoodTau(i2)) {
      const float deltar = sqrt(
				pow(muon_eta->at(i1) - tau_eta->at(i2), 2) +
				pow(Helper::DeltaPhi(muon_phi->at(i1), tau_phi->at(i2)), 2));
      if (deltar > 0.5) {
	validPair[i] = 1;
      }
    }
  }

// Find best muon based on pt
  int idx_1 = -1;
  float maxPt = -1;
  for(size_t i = 0; i < numComb; i++) {
    if(validPair[i] == 0) continue;
    const int tmp = comb.at(i).first;
    if(maxPt < muon_pt->at(tmp)) {
      maxPt = muon_pt->at(tmp);
      idx_1 = tmp;
    }
  }

// Find best tau based on iso
  int idx_2 = -1;
  float minIso = 999;
  for(size_t i = 0; i < numComb; i++) {
    if(validPair[i] == 0) continue;
    if(int(comb.at(i).first) != idx_1) continue;
    const int tmp = comb.at(i).second;
    if(minIso > tau_reliso_all->at(tmp)) {
      minIso = tau_reliso_all->at(tmp);
      idx_2 = tmp;
    }
  }

  vector<int> thegoodidx;
  thegoodidx.push_back(idx_1);
  thegoodidx.push_back(idx_2);

  return thegoodidx;


}//---------FindMuonTauPair










//-----------------------------------------------------------------
int main()
{
//-----------------------------------------------------------------

  gROOT->ProcessLine("#include<map>");

  /*
 * Compute event weights to be used for the respective datasets
 *
 * The event weight reweights the full dataset so that the sum of the weights
 * is equal to the expected number of events in data. The expectation is given by
 * multiplying the integrated luminosity of the data with the cross-section of
 * the process in the datasets divided by the number of simulated events.
 */
  //const float integratedLuminosity = 4.412 * 1000.0; // Run2012B only
  //const float integratedLuminosity = 7.055 * 1000.0; // Run2012C only
  const float integratedLuminosity = 11.467 * 1000.0; // Run2012B+C

  const float ggH_w = 19.6 / 476963.0 * integratedLuminosity;
  const float qqH_w = 1.55 / 491653.0 * integratedLuminosity;
  const float W1J_w =  6381.2 / 29784800.0 * integratedLuminosity;
  const float W2J_w =  2039.8 / 30693853.0 * integratedLuminosity;
  const float W3J_w =  612.5 / 15241144.0 * integratedLuminosity;
  const float TT_w = 225.2 / 6423106.0 * integratedLuminosity;
  const float ZLL_w = 3503.7 / 30458871.0 * integratedLuminosity;
  const float dataRunB_w = 1.0;
  const float dataRunC_w = 1.0;


 map<string, pair<string,float> > sampleNames;
 sampleNames.insert(make_pair("GluGluToHToTauTau",make_pair("ggH",ggH_w)));
 sampleNames.insert(make_pair("VBF_HToTauTau",make_pair("qqH",qqH_w)));
 sampleNames.insert(make_pair("W1JetsToLNu",make_pair("W1J",W1J_w)));
 sampleNames.insert(make_pair("W2JetsToLNu",make_pair("W2J",W2J_w)));
 sampleNames.insert(make_pair("W3JetsToLNu",make_pair("W3J",W3J_w)));
 sampleNames.insert(make_pair("TTbar",make_pair("TT",TT_w)));
 sampleNames.insert(make_pair("DYJetsToLL",make_pair("ZLL",ZLL_w)));
 sampleNames.insert(make_pair("Run2012B_TauPlusX",make_pair("dataRunB",dataRunB_w)));
 sampleNames.insert(make_pair("Run2012C_TauPlusX",make_pair("dataRunC",dataRunC_w)));




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
  dataRunB_npv->Write();
  dataRunB_eta_2->Write();
  dataRunB_m_vis->Write();

  dataRunC_npv->Write();
  dataRunC_eta_2->Write();
  dataRunC_m_vis->Write();

  ZLL_npv->Write();
  ZLL_eta_2->Write();
  ZLL_m_vis->Write();

  TT_npv->Write();
  TT_eta_2->Write();
  TT_m_vis->Write();

  W3J_npv->Write();
  W3J_eta_2->Write();
  W3J_m_vis->Write();

  W2J_npv->Write();
  W2J_eta_2->Write();
  W2J_m_vis->Write();

  W1J_npv->Write();
  W1J_eta_2->Write();
  W1J_m_vis->Write();

  qqH_npv->Write();
  qqH_eta_2->Write();
  qqH_m_vis->Write();

  ggH_npv->Write();
  ggH_eta_2->Write();
  ggH_m_vis->Write();

  //Save control region histos
  dataRunB_npv_cr->Write();
  dataRunB_eta_2_cr->Write();
  dataRunB_m_vis_cr->Write();

  dataRunC_npv_cr->Write();
  dataRunC_eta_2_cr->Write();
  dataRunC_m_vis_cr->Write();

  ZLL_npv_cr->Write();
  ZLL_eta_2_cr->Write();
  ZLL_m_vis_cr->Write();

  TT_npv_cr->Write();
  TT_eta_2_cr->Write();
  TT_m_vis_cr->Write();

  W3J_npv_cr->Write();
  W3J_eta_2_cr->Write();
  W3J_m_vis_cr->Write();

  W2J_npv_cr->Write();
  W2J_eta_2_cr->Write();
  W2J_m_vis_cr->Write();

  W1J_npv_cr->Write();
  W1J_eta_2_cr->Write();
  W1J_m_vis_cr->Write();

  qqH_npv_cr->Write();
  qqH_eta_2_cr->Write();
  qqH_m_vis_cr->Write();

  ggH_npv_cr->Write();
  ggH_eta_2_cr->Write();
  ggH_m_vis_cr->Write();


  hfile->Close();


  return 0;

}
