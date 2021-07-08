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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

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

  //declare the input tag for GsfElectronCollection
  edm::InputTag electronInput;
  
  // ----------member data ---------------------------
  
  
  TTree *mtree;
  int numelectron; //number of electrons in the event
  std::vector<float> electron_e;
  std::vector<float> electron_pt;
  std::vector<float> electron_px;
  std::vector<float> electron_py;
  std::vector<float> electron_pz;
  std::vector<float> electron_eta;
  std::vector<float> electron_phi;
  std::vector<float> electron_ch;
  std::vector<float> electron_iso;
  std::vector<bool> electron_isLoose;
  std::vector<bool> electron_isMedium;
  std::vector<bool> electron_isTight;
  std::vector<float> electron_dxy;
  std::vector<float> electron_dz;
  std::vector<float> electron_dxyError;
  std::vector<float> electron_dzError;
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
  edm::Service<TFileService> fs;
  mtree = fs->make<TTree>("Events", "Events");
  
  
  mtree->Branch("numberelectron",&numelectron);
  mtree->GetBranch("numberelectron")->SetTitle("number of electrons");
  mtree->Branch("electron_e",&electron_e);
  mtree->GetBranch("electron_e")->SetTitle("electron energy");
  mtree->Branch("electron_pt",&electron_pt);
  mtree->GetBranch("electron_pt")->SetTitle("electron transverse momentum");
  mtree->Branch("electron_px",&electron_px);
  mtree->GetBranch("electron_px")->SetTitle("electron momentum x-component");
  mtree->Branch("electron_py",&electron_py);
  mtree->GetBranch("electron_py")->SetTitle("electron momentum y-component");
  mtree->Branch("electron_pz",&electron_pz);
  mtree->GetBranch("electron_pz")->SetTitle("electron momentum z-component");
  mtree->Branch("electron_eta",&electron_eta);
  mtree->GetBranch("electron_eta")->SetTitle("electron pseudorapidity");
  mtree->Branch("electron_phi",&electron_phi);
  mtree->GetBranch("electron_phi")->SetTitle("electron polar angle");
  mtree->Branch("electron_ch",&electron_ch);
  mtree->GetBranch("electron_ch")->SetTitle("electron charge");
  mtree->Branch("electron_iso",&electron_iso);
  mtree->GetBranch("electron_iso")->SetTitle("electron isolation");
  mtree->Branch("electron_isLoose",&electron_isLoose);
  mtree->GetBranch("electron_isLoose")->SetTitle("electron tagged loose");
  mtree->Branch("electron_isMedium",&electron_isMedium);
  mtree->GetBranch("electron_isMedium")->SetTitle("electron tagged medium");
  mtree->Branch("electron_isTight",&electron_isTight);
  mtree->GetBranch("electron_isTight")->SetTitle("electron tagged tight");
  mtree->Branch("electron_dxy",&electron_dxy);
  mtree->GetBranch("electron_dxy")->SetTitle("electron transverse plane impact parameter (mm)");
  mtree->Branch("electron_dz",&electron_dz);
  mtree->GetBranch("electron_dz")->SetTitle("electron longitudinal impact parameter (mm)");
  mtree->Branch("electron_dxyError",&electron_dxyError);
  mtree->GetBranch("electron_dxyError")->SetTitle("electron transverse impact parameter uncertainty (mm)");
  mtree->Branch("electron_dzError",&electron_dzError);
  mtree->GetBranch("electron_dzError")->SetTitle("electron longitudinal impact parameter uncertainty (mm)");
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
   Handle<reco::ConversionCollection> hConversions;
   iEvent.getByLabel("allConversions", hConversions);
   Handle<reco::BeamSpot> bsHandle;
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
   electron_iso.clear();
   electron_isLoose.clear();
   electron_isMedium.clear();
   electron_isTight.clear();
   electron_dxy.clear();
   electron_dz.clear();
   electron_dxyError.clear();
   electron_dzError.clear();

   if(myelectrons.isValid()){
     // get the number of electrons in the event
     numelectron=myelectrons->size();
     for (reco::GsfElectronCollection::const_iterator itElec=myelectrons->begin(); itElec!=myelectrons->end(); ++itElec){
	     
       int missing_hits = itElec->gsfTrack()->trackerExpectedHitsInner().numberOfHits()-itElec->gsfTrack()->hitPattern().numberOfHits();
       bool passelectronveto = !ConversionTools::hasMatchedConversion(*itElec, hConversions, beamspot.position());
	     
       float el_pfIso = 999;
       if (itElec->passingPflowPreselection()) {
	 auto iso03 = itElec->pfIsolationVariables();
	 el_pfIso = (iso03.chargedHadronIso + iso03.neutralHadronIso + iso03.photonIso)/itElec->pt();
       } 
       auto trk = itElec->gsfTrack();
       bool isLoose = false, isMedium = false, isTight = false;
       if ( abs(itElec->eta()) <= 1.479 ) {   
	 if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.007 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.15 && 
	      itElec->sigmaIetaIeta()<.01 && itElec->hadronicOverEm()<.12 && 
	      abs(trk->dxy(pv))<.02 && abs(trk->dz(pv))<.2 && 
	      missing_hits<=1 && passelectronveto==true &&
	      abs(1/itElec->ecalEnergy()-1/(itElec->ecalEnergy()/itElec->eSuperClusterOverP()))<.05 && 
	      el_pfIso<.15){
	       
	   isLoose = true;
	       
	   if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.004 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.06 && abs(trk->dz(pv))<.1 ){
	     isMedium = true;
	             
	     if (abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.03 && missing_hits<=0 && el_pfIso<.10 ){
	       isTight = true;
	     }
	   }
	 }
       }
       else if ( abs(itElec->eta()) > 1.479 && abs(itElec->eta()) < 2.5 ) {
	 if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.009 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.1 && 
	      itElec->sigmaIetaIeta()<.03 && itElec->hadronicOverEm()<.1 && 
	      abs(trk->dxy(pv))<.02 && abs(trk->dz(pv))<.2 && 
	      missing_hits<=1 && el_pfIso<.15 && passelectronveto==true &&
	      abs(1/itElec->ecalEnergy()-1/(itElec->ecalEnergy()/itElec->eSuperClusterOverP()))<.05) {
	       
	   isLoose = true;
	       
	   if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.007 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.03 && abs(trk->dz(pv))<.1 ){
	     isMedium = true;
	             
	     if ( abs(itElec->deltaEtaSuperClusterTrackAtVtx())<.005 && abs(itElec->deltaPhiSuperClusterTrackAtVtx())<.02 && missing_hits<=0 && el_pfIso<.10 ){
	       isTight = true;
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
       electron_iso.push_back(el_pfIso);
       electron_isLoose.push_back(isLoose);
       electron_isMedium.push_back(isMedium);
       electron_isTight.push_back(isTight);
       electron_dxy.push_back(trk->dxy(pv));
       electron_dz.push_back(trk->dz(pv));
       electron_dxyError.push_back(trk->d0Error());
       electron_dzError.push_back(trk->dzError());
     }
  }

  mtree->Fill();
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void
ElectronAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
ElectronAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
ElectronAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
ElectronAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
ElectronAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

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
