// -*- C++ -*-
//
// Package:    ElectronAnalyzer
// Class:      ElectronAnalyzer
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//classes to extract electron information
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class ElectronAnalyzer : public edm::EDAnalyzer {
public:
  explicit ElectronAnalyzer(const edm::ParameterSet&);
  ~ElectronAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);

  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //declare a function to do the electron analysis
  void analyzeElectrons(const edm::Event& iEvent, const edm::Handle<reco::GsfElectronCollection> &electrons);

  //declare the input tag for GsfElectronCollection
  edm::InputTag electronInput;

  // ----------member data ---------------------------

  int numelectron; //number of electrons in the event

  TFile *mfile;
  TTree *mtree;
  
  std::vector<float> electron_e;
  std::vector<float> electron_pt;
  std::vector<float> electron_px;
  std::vector<float> electron_py;
  std::vector<float> electron_pz;
  std::vector<float> electron_eta;
  std::vector<float> electron_phi;
  std::vector<float> electron_ch;
  std::vector<bool> electron_isLoose;
  std::vector<bool> electron_isMedium;
  std::vector<bool> electron_isTight;
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

ElectronAnalyzer::ElectronAnalyzer(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
	electronInput = iConfig.getParameter<edm::InputTag>("InputCollection");
}

ElectronAnalyzer::~ElectronAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<reco::GsfElectronCollection> myelectrons;
   iEvent.getByLabel(electronInput, myelectrons);

   analyzeElectrons(iEvent,myelectrons);

   mtree->Fill();
   return;
}

void
ElectronAnalyzer::analyzeElectrons(const edm::Event& iEvent, const edm::Handle<reco::GsfElectronCollection> &electrons)
{
  using namespace edm;
  using namespace std;
  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);
  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel("offlineBeamSpot", bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();
  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
  math::XYZPoint pv(vertices->begin()->position());

  numelectron = 0;
  electron_e.clear();
  electron_pt.clear();
  electron_px.clear();
  electron_py.clear();
  electron_pz.clear();
  electron_eta.clear();
  electron_phi.clear();
  electron_ch.clear();
  electron_isLoose.clear();
  electron_isMedium.clear();
  electron_isTight.clear();
  
  if(electrons.isValid()){
     // get the number of electrons in the event
    numelectron=electrons->size();
    float pfIso = -999;
    for (reco::GsfElectronCollection::const_iterator itElec=electrons->begin(); itElec!=electrons->end(); ++itElec){
      int missing_hits = itElec->gsfTrack()->trackerExpectedHitsInner().numberOfHits()-itElec->gsfTrack()->hitPattern().numberOfHits();
      bool passelectronveto = !ConversionTools::hasMatchedConversion(*itElec, hConversions, beamspot.position());
      if (itElec->passingPflowPreselection()) {
        auto iso03 = itElec->pfIsolationVariables();
        float pfIso = (iso03.chargedHadronIso + iso03.neutralHadronIso + iso03.photonIso)/itElec->pt();
      } 
      auto trk = itElec->gsfTrack();
      bool el_isLoose = false;
      bool el_isMedium = false;
      bool el_isTight = false;
      if ( abs(itElec->eta()) <= 1.479 ) {   
	if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.007 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.15 && 
	     itElec->sigmaIetaIeta()<.01 && itElec->hadronicOverEm()<.12 && 
	     abs(trk->dxy(pv))<.02 && abs(trk->dz(pv))<.2 && 
	          missing_hits<=1 && pfIso<.15 && passelectronveto==true &&
	     abs(1/itElec->ecalEnergy()-1/(itElec->ecalEnergy()/itElec->eSuperClusterOverP()))<.05 ){
	    
          el_isLoose = true;
	    
	  if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.004 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.06 && abs(trk->dz(pv))<.1 ){
	    el_isMedium = true;
	        
	    if (abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.03 && missing_hits<=0 && pfIso<.10 ){
	      el_isTight = true;
	    }
	  }
	}
      }
      else if ( abs(itElec->eta()) > 1.479 && abs(itElec->eta()) < 2.5 ) {
        if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.009 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.1 && 
	     itElec->sigmaIetaIeta()<.03 && itElec->hadronicOverEm()<.1 && 
	     abs(trk->dxy(pv))<.02 && abs(trk->dz(pv))<.2 && 
	          missing_hits<=1 && pfIso<.15 && passelectronveto==true &&
             abs(1/itElec->ecalEnergy()-1/(itElec->ecalEnergy()/itElec->eSuperClusterOverP()))<.05) {
	    
          el_isLoose = true;
	    
	  if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.007 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.03 && abs(trk->dz(pv))<.1 ){
	    el_isMedium = true;
	        
	    if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.005 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.02 && missing_hits<=0 && pfIso<.10 ){
	      el_isTight = true;
	    }
	  }
        }
      }

      electron_e.push_back(itElec->energy());
      electron_pt.push_back(itElec->pt());
      electron_px.push_back(itElec->px());
      electron_py.push_back(itElec->py());
      electron_pz.push_back(itElec->pz());
      electron_eta.push_back(itElec->eta());
      electron_phi.push_back(itElec->phi());
      electron_ch.push_back(itElec->charge());
      electron_isLoose.push_back(el_isLoose);
      electron_isMedium.push_back(el_isMedium);
      electron_isTight.push_back(el_isTight);
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void
ElectronAnalyzer::beginJob()
{

mfile = new TFile("ElectronInfo.root","RECREATE");
mtree = new TTree("mtree","Electron information");

  mtree->Branch("numberelectron",&numelectron);
  mtree->Branch("electron_e",&electron_e);
  mtree->Branch("electron_pt",&electron_pt);
  mtree->Branch("electron_px",&electron_px);
  mtree->Branch("electron_py",&electron_py);
  mtree->Branch("electron_pz",&electron_pz);
  mtree->Branch("electron_eta",&electron_eta);
  mtree->Branch("electron_phi",&electron_phi);
  mtree->Branch("electron_ch",&electron_ch);
}

// ------------ method called once each job just after ending the event loop  ------------
void
ElectronAnalyzer::endJob()
{
  mfile->Write();
}

// ------------ method called when starting to processes a run  ------------
void
ElectronAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
ElectronAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void
ElectronAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
ElectronAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronAnalyzer);
