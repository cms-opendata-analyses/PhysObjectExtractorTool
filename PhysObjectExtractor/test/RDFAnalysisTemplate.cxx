#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TH1D.h"
#include <iostream>
#include <vector>
//#include "TLorentzVector.h"

//ATENTION!!!: It needs ROOT version > 6.22

//Compile me with:
//g++ -g -O3 -Wall -Wextra -Wpedantic -o RDFAnalysis RDFAnalysisTemplate.cxx $(root-config --cflags --libs)

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
  cout<<"Running the RDFrameAnalysis() ...."<<endl;
  //This is a more sophisticated analyisis using a columnar approach
  RDFrameAnalysis();

  

}
