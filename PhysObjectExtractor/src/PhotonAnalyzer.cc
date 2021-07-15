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
#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>
#include<tuple>

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

  struct AEff {
    float CH_AEff;
    float NH_AEff;
    float Ph_AEff;
  };
  virtual struct AEff effectiveArea0p3cone(float eta);  
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
  std::vector<float> photon_chIso;
  std::vector<float> photon_nhIso;
  std::vector<float> photon_phIso;
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
  mtree->GetBranch("numberphoton")->SetTitle("number of photons");
  mtree->Branch("photon_e",&photon_e);
  mtree->GetBranch("photon_e")->SetTitle("photon energy");
  mtree->Branch("photon_pt",&photon_pt);
  mtree->GetBranch("photon_pt")->SetTitle("photon transverse momentum");
  mtree->Branch("photon_px",&photon_px);
  mtree->GetBranch("photon_px")->SetTitle("photon momentum x-component");
  mtree->Branch("photon_py",&photon_py);
  mtree->GetBranch("photon_py")->SetTitle("photon momentum y-component");
  mtree->Branch("photon_pz",&photon_pz);
  mtree->GetBranch("photon_pz")->SetTitle("photon momentum z-component");
  mtree->Branch("photon_eta",&photon_eta);
  mtree->GetBranch("photon_eta")->SetTitle("photon pseudorapidity");
  mtree->Branch("photon_phi",&photon_phi);
  mtree->GetBranch("photon_phi")->SetTitle("photon polar angle");
  mtree->Branch("photon_ch",&photon_ch);
  mtree->GetBranch("photon_ch")->SetTitle("photon charge");
  mtree->Branch("photon_chIso",&photon_chIso);
  mtree->GetBranch("photon_chIso")->SetTitle("corrected photon charged hadron isolation");
  mtree->Branch("photon_nhIso",&photon_nhIso);
  mtree->GetBranch("photon_nhIso")->SetTitle("corrected photon neutral hadron isolation");  
  mtree->Branch("photon_phIso",&photon_phIso);
  mtree->GetBranch("photon_phIso")->SetTitle("corrected photon isolation from other photons");
  mtree->Branch("photon_isLoose",&photon_isLoose);
  mtree->GetBranch("photon_isLoose")->SetTitle("photon tagged loose");
  mtree->Branch("photon_isMedium",&photon_isMedium);
  mtree->GetBranch("photon_isMedium")->SetTitle("photon tagged medium");
  mtree->Branch("photon_isTight",&photon_isTight);
  mtree->GetBranch("photon_isTight")->SetTitle("photon tagged tight");
}


PhotonAnalyzer::~PhotonAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

struct PhotonAnalyzer::AEff
PhotonAnalyzer::effectiveArea0p3cone(float eta)
{
  struct AEff A;
  if(fabs(eta) >2.4) {
    A.CH_AEff = 0.012;
    A.NH_AEff = 0.072;
    A.Ph_AEff = 0.266;
  }
  else if(fabs(eta) >2.3) {
    A.CH_AEff = 0.020;
    A.NH_AEff = 0.039;
    A.Ph_AEff = 0.260;
  }
  else if(fabs(eta) >2.2) {
    A.CH_AEff = 0.016;
    A.NH_AEff = 0.024;
    A.Ph_AEff = 0.262;
  }
  else if(fabs(eta) >2.0) {
    A.CH_AEff = 0.012;
    A.NH_AEff = 0.015;
    A.Ph_AEff = 0.216;
  }
  else if(fabs(eta) >1.479) {
    A.CH_AEff = 0.014;
    A.NH_AEff = 0.039;
    A.Ph_AEff = 0.112;
  }
  else if(fabs(eta) >0.1) {
    A.CH_AEff = 0.010;
    A.NH_AEff = 0.057;
    A.Ph_AEff = 0.130;
  }
  else {
    A.CH_AEff = 0.012;
    A.NH_AEff = 0.030;
    A.Ph_AEff = 0.148;
  }
  return A;
}


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
   Handle<reco::PFCandidateCollection> pfCands;
   iEvent.getByLabel("particleFlow", pfCands);
   Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel("offlinePrimaryVertices", vertices);

   numphoton = 0;
   photon_e.clear();
   photon_pt.clear();
   photon_px.clear();
   photon_py.clear();
   photon_pz.clear();
   photon_eta.clear();
   photon_phi.clear();
   photon_ch.clear();
   photon_chIso.clear();
   photon_nhIso.clear();
   photon_phIso.clear();
   photon_isLoose.clear();
   photon_isMedium.clear();
   photon_isTight.clear();

   if(myphotons.isValid()){
     // get the number of photons in the event
     numphoton=myphotons->size();
     for (reco::PhotonCollection::const_iterator itphoton=myphotons->begin(); itphoton!=myphotons->end(); ++itphoton){
       bool passelectronveto = !ConversionTools::hasMatchedPromptElectron(itphoton->superCluster(), electrons, hConversions, beamspot.position());
       double scEta = (itphoton)->superCluster()->eta();
       struct PhotonAnalyzer::AEff aEff = effectiveArea0p3cone(scEta);
       double ph_hOverEm = itphoton->hadTowOverEm();
       double ph_sigIetaIeta = itphoton->sigmaIetaIeta();
       PFIsolationEstimator isolator;
       isolator.initializePhotonIsolation(kTRUE);
       isolator. setConeSize(0.3);
       const reco::VertexRef vertex(vertices, 0);
       //const reco::PFCandidate pfph(itphoton->charge(),itphoton->p4(),gamma);
       const reco::Photon &thephoton = *itphoton;
       isolator.fGetIsolation(&thephoton, pfCands.product(), vertex, vertices);
       double corrPFCHIso = std::max(isolator.getIsolationCharged() - rhoIso * aEff.CH_AEff, 0.)/itphoton->pt();
       double corrPFNHIso = std::max(isolator.getIsolationNeutral() - rhoIso * aEff.NH_AEff, 0.)/itphoton->pt();
       double corrPFPhIso = std::max(isolator.getIsolationPhoton() - rhoIso * aEff.Ph_AEff, 0.)/itphoton->pt();
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
       photon_chIso.push_back(corrPFCHIso);
       photon_nhIso.push_back(corrPFNHIso);
       photon_phIso.push_back(corrPFPhIso);
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
