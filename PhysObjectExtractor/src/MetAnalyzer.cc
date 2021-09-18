// -*- C++ -*-
//
// Package:    MetAnalyzer
// Class:      MetAnalyzer
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//classes to extract PFMET information
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class MetAnalyzer : public edm::EDAnalyzer {
public:
  explicit MetAnalyzer(const edm::ParameterSet&);
  ~MetAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  //declare the input tag for PFMETCollection
  edm::InputTag metInput;
  bool doPat = false;
  edm::InputTag metInputPat;
  
  // ----------member data ---------------------------
  
  TTree *mtree;
  float met_e;
  float met_pt;
  float met_px;
  float met_py;
  float met_phi;
  float met_significance;
  float met_rawe;
  float met_rawpt;
  float met_rawphi;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

MetAnalyzer::MetAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  metInput = iConfig.getParameter<edm::InputTag>("InputCollection");
  doPat = iConfig.getParameter<bool>("doPat");
  if (doPat) metInputPat = iConfig.getParameter<edm::InputTag>("InputCollectionPat");
  edm::Service<TFileService> fs;
  mtree = fs->make<TTree>("Events", "Events");
  
  mtree->Branch("met_e",&met_e,"met_e/F");	
  mtree->GetBranch("met_e")->SetTitle("Sum of transverse energy (corrected if using PAT) [GeV]");
  mtree->Branch("met_pt",&met_pt,"met_pt/F");
  mtree->GetBranch("met_pt")->SetTitle("Missing transverse momentum (corrected if using PAT) [GeV]");
  mtree->Branch("met_px",&met_px,"met_px/F");
  mtree->GetBranch("met_px")->SetTitle("Missing x momentum (corrected if using PAT) [GeV]");
  mtree->Branch("met_py",&met_py,"met_py/F");
  mtree->GetBranch("met_py")->SetTitle("Missing y momentum (corrected if using PAT) [GeV]");
  mtree->Branch("met_phi",&met_phi,"met_phi/F");
  mtree->GetBranch("met_phi")->SetTitle("Missing momentum azimuthal angle (corrected if using PAT)");
  mtree->Branch("met_significance",&met_significance,"met_significance/F");
  mtree->GetBranch("met_significance")->SetTitle("Missing transverse momentum significance");
  if(doPat) { 
    mtree->Branch("met_rawpt",&met_rawpt,"met_rawpt/F");
    mtree->GetBranch("met_rawpt")->SetTitle("Missing transverse momentum (uncorrected) [GeV]");
    mtree->Branch("met_rawphi",&met_rawphi,"met_rawphi/F");
    mtree->GetBranch("met_phi")->SetTitle("Missing momentum azimuthal angle (uncorrected)");
    mtree->Branch("met_rawe",&met_rawe,"met_rawe/F");
    mtree->GetBranch("met_rawe")->SetTitle("Sum of transverse energy (uncorrected) [GeV]");
  }
}

MetAnalyzer::~MetAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
MetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<reco::PFMETCollection> mymets;
   iEvent.getByLabel(metInput, mymets);
   
  if(mymets.isValid()){
    met_e = mymets->begin()->sumEt();
    met_pt = mymets->begin()->pt();
    met_px = mymets->begin()->px();
    met_py = mymets->begin()->py();
    met_phi = mymets->begin()->phi();
    met_significance = mymets->begin()->significance();

    if(doPat){
      met_rawe = met_e;
      met_rawpt = met_pt;
      met_rawphi = met_phi;
    }
  }
  if(doPat){
    Handle<reco::PFMETCollection> patmets;
    iEvent.getByLabel(metInputPat, patmets);

    if(patmets.isValid()){
      met_e = patmets->begin()->sumEt();
      met_pt = patmets->begin()->pt();
      met_px = patmets->begin()->px();
      met_py = patmets->begin()->py();
      met_phi = patmets->begin()->phi();
    }
  }
  
  mtree->Fill();
  return;
   
}

// ------------ method called once each job just before starting event loop  ------------
void
MetAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
MetAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
MetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
MetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
MetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
MetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MetAnalyzer);
