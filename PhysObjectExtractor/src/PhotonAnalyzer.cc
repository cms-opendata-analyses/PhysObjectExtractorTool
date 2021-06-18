// -*- C++ -*-
//
// Package:    PhotonAnalyzer
// Class:      PhotonAnalyzer
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

//classes to extract Photon information
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class PhotonAnalyzer : public edm::EDAnalyzer {
public:
  explicit PhotonAnalyzer(const edm::ParameterSet&);
  ~PhotonAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  //declare the input tag for PhotonCollection
  edm::InputTag photonInput;

  // ----------member data ---------------------------
  
  TTree *mtree;
  int numphoton; //number of photons in the event
  std::vector<float> photon_e;
  std::vector<float> photon_pt;
  std::vector<float> photon_px;
  std::vector<float> photon_py;
  std::vector<float> photon_pz;
  std::vector<float> photon_eta;
  std::vector<float> photon_phi;
  std::vector<float> photon_ch;
  std::vector<float> photon_iso;
  std::vector<bool>  photon_isLoose;
  std::vector<bool>  photon_isMedium;
  std::vector<bool>  photon_isTight;
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

PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& iConfig)

{
//now do what ever initialization is needed
  photonInput = iConfig.getParameter<edm::InputTag>("InputCollection");
  edm::Service<TFileService> fs;
  mtree = fs->make<TTree>("Events", "Events");
  
  
  mtree->Branch("numberphoton",&numphoton);
  mtree->Branch("photon_e",&photon_e);
  mtree->Branch("photon_pt",&photon_pt);
  mtree->Branch("photon_px",&photon_px);
  mtree->Branch("photon_py",&photon_py);
  mtree->Branch("photon_pz",&photon_pz);
  mtree->Branch("photon_eta",&photon_eta);
  mtree->Branch("photon_phi",&photon_phi);
  mtree->Branch("photon_ch",&photon_ch);
  mtree->Branch("photon_iso",&photon_iso);
  mtree->Branch("photon_isLoose",&photon_isLoose);
  mtree->Branch("photon_isMedium",&photon_isMedium);
  mtree->Branch("photon_isTight",&photon_isTight);
}


PhotonAnalyzer::~PhotonAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<reco::PhotonCollection> myphotons;
   iEvent.getByLabel(photonInput, myphotons);
   Handle<reco::GsfElectronCollection> electrons;
   iEvent.getByLabel(InputTag("gsfElectrons"), electrons);
   Handle<reco::ConversionCollection> hConversions;
   iEvent.getByLabel("allConversions", hConversions);
   Handle<reco::BeamSpot> bsHandle;
   iEvent.getByLabel("offlineBeamSpot", bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();
   Handle<double> rhoHandle;
   iEvent.getByLabel(InputTag("fixedGridRhoAll"), rhoHandle);
   double rhoIso = std::max(*(rhoHandle.product()), 0.0);   

   numphoton = 0;
   photon_e.clear();
   photon_pt.clear();
   photon_px.clear();
   photon_py.clear();
   photon_pz.clear();
   photon_eta.clear();
   photon_phi.clear();
   photon_ch.clear();
   photon_iso.clear();
   photon_isLoose.clear();
   photon_isMedium.clear();
   photon_isTight.clear();

   if(myphotons.isValid()){
     // get the number of photons in the event
     numphoton=myphotons->size();
     for (reco::PhotonCollection::const_iterator itphoton=myphotons->begin(); itphoton!=myphotons->end(); ++itphoton){
       bool passelectronveto = !ConversionTools::hasMatchedPromptElectron(itphoton->superCluster(), electrons, hConversions, beamspot.position());
       double scEta = (itphoton)->superCluster()->eta();
       double CH_AEff, NH_AEff, Ph_AEff;
       if(fabs(scEta) >2.4) {
	 CH_AEff = 0.012;
	 NH_AEff = 0.072;
	 Ph_AEff = 0.266;
       }
       else if(fabs(scEta) >2.3) {
	 CH_AEff = 0.020;
	 NH_AEff = 0.039;
	 Ph_AEff = 0.260;
       } 
       else if(fabs(scEta) >2.2) {
	 CH_AEff = 0.016;
	 NH_AEff = 0.024;
	 Ph_AEff = 0.262;
       } 
       else if(fabs(scEta) >2.0) {
	 CH_AEff = 0.012;
	 NH_AEff = 0.015;
	 Ph_AEff = 0.216;
       } 
       else if(fabs(scEta) >1.479) {
	 CH_AEff = 0.014;
	 NH_AEff = 0.039;
	 Ph_AEff = 0.112;
       } 
       else if(fabs(scEta) >0.1) {
	 CH_AEff = 0.010;
	 NH_AEff = 0.057;
	 Ph_AEff = 0.130;
       } 
       else {
	 CH_AEff = 0.012;
	 NH_AEff = 0.030;
	 Ph_AEff = 0.148;
       } 
       double ph_hOverEm = itphoton->hadTowOverEm();
       double ph_sigIetaIeta = itphoton->sigmaIetaIeta();
       double ph_photonIso = itphoton->photonIso();
       double corrPFCHIso = std::max(itphoton->chargedHadronIso() - rhoIso * CH_AEff, 0.);
       double corrPFNHIso = std::max(itphoton->neutralHadronIso() - rhoIso * NH_AEff, 0.);
       double corrPFPhIso = std::max(itphoton->photonIso() - rhoIso * Ph_AEff, 0.);
       bool isLoose = false, isMedium = false, isTight = false;
       if ( itphoton->eta() <= 1.479 ){
	 if ( ph_hOverEm<.05 && ph_sigIetaIeta<.012 && 
	      corrPFCHIso<2.6 && corrPFNHIso<(3.5+.04*itphoton->pt()) && 
	      corrPFPhIso<(1.3+.005*itphoton->pt()) && passelectronveto==true) {
	   isLoose = true;

	   if ( ph_sigIetaIeta<.011 && corrPFCHIso<1.5 && corrPFNHIso<(1.0+.04*itphoton->pt()) && corrPFPhIso<(.7+.005*itphoton->pt())){
	     isMedium = true;

	     if ( corrPFCHIso<.7 && corrPFNHIso<(.4+.04*itphoton->pt()) && corrPFPhIso<(.5+0.005*itphoton->pt()) ){
	       isTight = true;
	     }
	   }
	 }
       }
       else if ( itphoton->eta() > 1.479 && itphoton->eta() < 2.5 ) {
	 if ( ph_hOverEm<.05 && ph_sigIetaIeta<.034 && corrPFCHIso<2.3 && corrPFNHIso<(2.9+.04*itphoton->pt()) && passelectronveto==true ){
	   isLoose = true;
	           
	   if ( ph_sigIetaIeta<.033 && corrPFCHIso<1.2 && corrPFNHIso<(1.5+.04*itphoton->pt()) && corrPFPhIso<(1.0+.005*itphoton->pt())) {
	     isMedium = true;

	     if ( ph_sigIetaIeta<0.031 && corrPFCHIso<0.5){
	       isTight = true;
	     }
	   }
	 }
       }
       
       photon_e.push_back(itphoton->energy());
       photon_pt.push_back(itphoton->pt());
       photon_px.push_back(itphoton->px());
       photon_py.push_back(itphoton->py());
       photon_pz.push_back(itphoton->pz());
       photon_eta.push_back(itphoton->eta());
       photon_phi.push_back(itphoton->phi());
       photon_ch.push_back(itphoton->charge());
       photon_iso.push_back(ph_photonIso);
       photon_isLoose.push_back(isLoose);
       photon_isMedium.push_back(isMedium);
       photon_isTight.push_back(isTight);
     }
  }
  
  mtree->Fill();
  return;
  
}

// ------------ method called once each job just before starting event loop  ------------
void
PhotonAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
PhotonAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
PhotonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
PhotonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
PhotonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
PhotonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonAnalyzer);
