#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TH1D.h"
#include <iostream>
#include <vector>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <chrono>

//Compile with:
//g++ -o test POET_test.cxx $(root-config --cflags --libs)


void EventLoopAnalysis(){

  using namespace std;

  TH1F* h_GenPart_pt = new TH1F("GenPart_pt","Generator particle transverse momentum",15,0,15);
  h_GenPart_pt->GetXaxis()->SetTitle("p_{T}/GeV");
  TH1F* h_GenPart_num = new TH1F("GenPart_num","Number of generator particles",100,0,1000);
  h_GenPart_num->GetXaxis()->SetTitle("total number");
  TH1F* h_GenPart_phi = new TH1F("GenPart_phi","Generator particle azimuthal angle",30,-3,3);
  h_GenPart_phi->GetXaxis()->SetTitle("phi");
  TH1F* h_GenPart_status = new TH1F("GenPart_status","Generator particle status",5,0,5);
  h_GenPart_status->GetXaxis()->SetTitle("status");
   
  TFile* infile = new TFile("myoutput.root", "READ");
  
  TTree* tevent = (TTree*)infile->Get("myevents/Events");
  TTree* tgenparticles = (TTree*)infile->Get("mygenparticle/Events");

  Int_t run;
  ULong64_t event;
   
  Int_t numGenPart=0;
  vector<float>*  GenPart_pt=0;
  vector<float>*  GenPart_phi=0;
  vector<int>*  GenPart_status=0;
   
  tevent->SetBranchAddress("run",&run);
  tevent->SetBranchAddress("event",&event);
   
  tgenparticles->SetBranchAddress("numGenPart",&numGenPart);
  tgenparticles->SetBranchAddress("GenPart_pt",&GenPart_pt);
  tgenparticles->SetBranchAddress("GenPart_phi",&GenPart_phi);
  tgenparticles->SetBranchAddress("GenPart_status",&GenPart_status);
   
  tevent->AddFriend(tgenparticles);
   
  Long64_t nentries = tevent->GetEntries();

  for (Long64_t i=0;i<nentries;i++) {
    
    tevent->GetEntry(i);
     
    h_GenPart_num->Fill(numGenPart);
    for (Int_t j=0; j<numGenPart;++j)
    {
      h_GenPart_pt->Fill(GenPart_pt->at(j));
      h_GenPart_phi->Fill(GenPart_phi->at(j));
      h_GenPart_status->Fill(GenPart_status->at(j));
    }
      
  }
  
  gStyle->SetOptStat(0);
	TCanvas *GenPart__Canvas=new TCanvas("GenPart","GenPart",1600,800);
  GenPart__Canvas->Divide(2,2);
  GenPart__Canvas->cd(1);
  h_GenPart_pt->Draw();
  GenPart__Canvas->cd(2);
  h_GenPart_num->Draw();
  GenPart__Canvas->cd(3);
  h_GenPart_phi->Draw();
  GenPart__Canvas->cd(4);
  h_GenPart_status->Draw();
  GenPart__Canvas->SaveAs("GenPart_loop.png");

  delete infile;
  
}



void RDFrameAnalysis(){
  using namespace ROOT;
  using namespace std;
  
  TFile* infile = new TFile("file:myoutput.root", "READ");
  
  TTree* tevent = (TTree*)infile->Get("myevents/Events");
  TTree* tgenparticles = (TTree*)infile->Get("mygenparticle/Events");
  
  tevent->AddFriend(tgenparticles);
  
  RDataFrame d(*tevent);

  auto h_GenPart_pt = d.Define("fgenparticles","GenPart_pt").Histo1D({"GenPart_pt", "Generator particle transverse momentum", 15u, 0., 15.},"fgenparticles");
  h_GenPart_pt->GetXaxis()->SetTitle("p_{T}/GeV");
  auto h_GenPart_num = d.Define("fgenparticles","numGenPart").Histo1D({"numGenPart", "Number of generator particles", 100u, 0., 1000.},"fgenparticles");
  h_GenPart_num->GetXaxis()->SetTitle("total number");
  auto h_GenPart_phi = d.Define("fgenparticles","GenPart_phi").Histo1D({"GenPart_phi", "Generator particle azimuthal angle", 30u, -3., 3.},"fgenparticles");
  h_GenPart_phi->GetXaxis()->SetTitle("phi");
  auto h_GenPart_status = d.Define("fgenparticles","GenPart_status").Histo1D({"GenPart_pt", "Generator particle status", 5u, 0., 5.},"fgenparticles");
  h_GenPart_status->GetXaxis()->SetTitle("status");

  gStyle->SetOptStat(0);
	TCanvas *GenPart__Canvas=new TCanvas("GenPart","GenPart",1600,800);
  GenPart__Canvas->Divide(2,2);
  GenPart__Canvas->cd(1);
  h_GenPart_pt->Draw();
  GenPart__Canvas->cd(2);
  h_GenPart_num->Draw();
  GenPart__Canvas->cd(3);
  h_GenPart_phi->Draw();
  GenPart__Canvas->cd(4);
  h_GenPart_status->Draw();
  GenPart__Canvas->SaveAs("GenPart_RDFrame.png");

} 

int main (){
  using namespace std;
  chrono::steady_clock sc;

  cout<<"Running the EventLoopAnalysis() ...."<<endl;
  auto start1 = sc.now();
  EventLoopAnalysis();          
  auto end1 = sc.now();        
  auto time_span1 = static_cast<chrono::duration<double>>(end1 - start1);    
  cout<<"Event loop analysis took: "<<time_span1.count()<<" seconds "<<endl;
  
  cout<<"Running the RDFrameAnalysis() ...."<<endl;
  auto start2 = sc.now();
  RDFrameAnalysis();      
  auto end2 = sc.now();        
  auto time_span2 = static_cast<chrono::duration<double>>(end2 - start2);    
  cout<<"RDFrame analysis took: "<<time_span2.count()<<" seconds "<<endl;
  
}
