// -*- C++ -*-
//
// Package:    PhysicsObjectsInfo
// Class:      PhysicsObjectsInfo
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

//classes to extract jet information
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//classes to extract met information
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

//classes to extract muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

//classes to extract photon information
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

//classes to extract tau information
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

//classes to extract trigger information
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class PhysicsObjectsInfo : public edm::EDAnalyzer {
   public:
      explicit PhysicsObjectsInfo(const edm::ParameterSet&);
      ~PhysicsObjectsInfo();

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

//declare a function to do the jet analysis
      void analyzeJets(const edm::Event& iEvent, const edm::Handle<reco::PFJetCollection> &jets);

//declare a function to do the met analysis
      void analyzeMets(const edm::Event& iEvent, const edm::Handle<reco::PFMETCollection> &mets);

//declare a function to do the muon analysis
      void analyzeMuons(const edm::Event& iEvent, const edm::Handle<reco::MuonCollection> &muons);

//declare a function to do the photon analysis
      void analyzePhotons(const edm::Event& iEvent, const edm::Handle<reco::PhotonCollection> &photons);

//declare a function to do the tau analysis
      void analyzeTaus(const edm::Event& iEvent, const edm::Handle<reco::PFTauCollection> &taus);

//declare a function to do the trigger analysis
      void analyzeTriggObject(const edm::Event& iEvent, const edm::Handle<trigger::TriggerEvent> &trigEvent, const edm::InputTag &trigEventTag_);

//declare the input tag for the electron collection
      edm::InputTag electronInput;
//declare the input tag for the jet collection
      edm::InputTag jetInput;
//declare the input tag for the met collection
      edm::InputTag metInput;
//declare the input tag for the muon collection
      edm::InputTag muonInput;
//declare the input tag for the photon collection
      edm::InputTag photonInput;
//declare the input tag for the tau collection
      edm::InputTag tauInput;
//declare de filter (module) of the trigger
      std::string   filterName_;

//declare the variables for save data
	TFile *mfile;
	TTree *mtree;

//declare variables that will go into the root files
// ----------electron member data ---------------------------
	int numelectron; //number of electrons in the event

	  std::vector<float> electron_e;
  	std::vector<float> electron_pt;
  	std::vector<float> electron_px;
  	std::vector<float> electron_py;
  	std::vector<float> electron_pz;
  	std::vector<float> electron_eta;
  	std::vector<float> electron_phi;
  	std::vector<float> electron_ch;

// ----------jet member data ---------------------------
	int numjet; //number of jets in the event

	  std::vector<float> jet_e;
  	std::vector<float> jet_pt;
  	std::vector<float> jet_px;
  	std::vector<float> jet_py;
  	std::vector<float> jet_pz;
  	std::vector<float> jet_eta;
  	std::vector<float> jet_phi;
  	std::vector<float> jet_ch;

// ----------met member data ---------------------------
	int nummet; //number of mets in the event

	  std::vector<float> met_e;
  	std::vector<float> met_pt;
  	std::vector<float> met_px;
  	std::vector<float> met_py;
  	std::vector<float> met_phi;
  	std::vector<float> met_ch;

// ----------muon member data ---------------------------
	int nummuon; //number of muons in the event

	  std::vector<float> muon_e;
  	std::vector<float> muon_pt;
  	std::vector<float> muon_px;
  	std::vector<float> muon_py;
  	std::vector<float> muon_pz;
  	std::vector<float> muon_eta;
  	std::vector<float> muon_phi;
  	std::vector<float> muon_ch;

// ----------photon member data ---------------------------
	int numphoton; //number of photons in the event

	  std::vector<float> photon_e;
  	std::vector<float> photon_pt;
  	std::vector<float> photon_px;
  	std::vector<float> photon_py;
  	std::vector<float> photon_pz;
  	std::vector<float> photon_eta;
  	std::vector<float> photon_phi;
  	std::vector<float> photon_ch;

// ----------tau member data -------------------------------
  int numtau; //number of taus in the event

    std::vector<float> tau_e;
    std::vector<float> tau_pt;
    std::vector<float> tau_px;
    std::vector<float> tau_py;
    std::vector<float> tau_pz;
    std::vector<float> tau_eta;
    std::vector<float> tau_phi;
    std::vector<float> tau_ch;

// ----------trigger member data ---------------------------
	int numtrigobj; //number of trigger objects in the event

	  std::vector<float> trigobj_e;
  	std::vector<float> trigobj_pt;
  	std::vector<float> trigobj_px;
  	std::vector<float> trigobj_py;
  	std::vector<float> trigobj_pz;
  	std::vector<float> trigobj_eta;
  	std::vector<float> trigobj_phi;
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

PhysicsObjectsInfo::PhysicsObjectsInfo(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
//this takes the type of object (input tag) specified in the configuration python file, corresponding to each container.
//https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable
	electronInput = iConfig.getParameter<edm::InputTag>("ElectronInputCollection");
	jetInput = iConfig.getParameter<edm::InputTag>("JetInputCollection");
	metInput = iConfig.getParameter<edm::InputTag>("MetInputCollection");
	muonInput = iConfig.getParameter<edm::InputTag>("MuonInputCollection");
	photonInput = iConfig.getParameter<edm::InputTag>("PhotonInputCollection");
  tauInput = iConfig.getParameter<edm::InputTag>("TauInputCollection");
//this take the filter (module) of the trigger
	filterName_ = iConfig.getParameter<std::string>("filterName");
}

PhysicsObjectsInfo::~PhysicsObjectsInfo()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhysicsObjectsInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

//Declare the handle (container) to store electrons.
   Handle<reco::GsfElectronCollection> myelectrons;
//This is where the desired object is extracted from the EDM file,
//since it is composed of several branches with all the resulting objects of each event.
//Then, objects of the input tag type placed in the configuration file, in this case "gsfElectrons",
//are stored in the variable "myelectrons". It must be ensured that the input tag goes with its
//corresponding container (handle). In addition, a container can accept different types of input tag,
//for this reason, this tool is structured in this way, since it allows changing the input tag
//from the python file (configuration file) without the need to recompile the code.
   iEvent.getByLabel(electronInput, myelectrons);
//Here, pass the object collection and their event to new function, in order to get some information about them.
   analyzeElectrons(iEvent,myelectrons);

//The following functions follow the same process:
//Declare the handle (container) to store jets.
   Handle<reco::PFJetCollection> myjets;
   iEvent.getByLabel(jetInput, myjets);
   analyzeJets(iEvent,myjets);

//Declare the handle (container) to store mets.
   Handle<reco::PFMETCollection> mymets;
   iEvent.getByLabel(metInput, mymets);
   analyzeMets(iEvent,mymets);

//Declare the handle (container) to store muons.
   Handle<reco::MuonCollection> mymuons;
   iEvent.getByLabel(muonInput, mymuons);
   analyzeMuons(iEvent,mymuons);

//Declare the handle (container) to store photons.
   Handle<reco::PhotonCollection> myphotons;
   iEvent.getByLabel(photonInput, myphotons);
   analyzePhotons(iEvent,myphotons);

//Declare the handle (container) to store taus.
   Handle<reco::PFTauCollection> mytaus;
   iEvent.getByLabel(tauInput, mytaus);
   analyzeTaus(iEvent,mytaus);

//Declare the handle (container) to store trigger objects.
  InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); //make sure have correct process on MC
//data process=HLT, MC depends, Spring11 is REDIGI311X
	Handle<trigger::TriggerEvent> mytrigEvent;
	iEvent.getByLabel(trigEventTag,mytrigEvent);
	analyzeTriggObject(iEvent,mytrigEvent,trigEventTag);

//Now, the information is stored.
   mtree->Fill();
   return;
}

void
PhysicsObjectsInfo::analyzeElectrons(const edm::Event& iEvent, const edm::Handle<reco::GsfElectronCollection> &electrons)
{
	//Since there are several events, each variable is cleaned at the beginning of the analysis.
	  numelectron = 0;
	  electron_e.clear();
	  electron_pt.clear();
	  electron_px.clear();
	  electron_py.clear();
	  electron_pz.clear();
	  electron_eta.clear();
	  electron_phi.clear();
	  electron_ch.clear();
// Check if the object collection is valid.
  if(electrons.isValid()){
     // get the number of electrons in the event
     numelectron=(*electrons).size();
	  //A cycle is created that sweep over all the electrons of the event.
        for (reco::GsfElectronCollection::const_iterator itElec=electrons->begin(); itElec!=electrons->end(); ++itElec){
            //get general information
      	    electron_e.push_back(itElec->energy());
      	    electron_pt.push_back(itElec->pt());
      	    electron_px.push_back(itElec->px());
      	    electron_py.push_back(itElec->py());
      	    electron_pz.push_back(itElec->pz());
      	    electron_eta.push_back(itElec->eta());
      	    electron_phi.push_back(itElec->phi());
      	    electron_ch.push_back(itElec->charge());
            //Here we can choose electron with specific cuts and get more specific information.
        }
  }
}

void
PhysicsObjectsInfo::analyzeJets(const edm::Event& iEvent, const edm::Handle<reco::PFJetCollection> &jets)
{
	  numjet = 0;
	  jet_e.clear();
	  jet_pt.clear();
	  jet_px.clear();
	  jet_py.clear();
	  jet_pz.clear();
	  jet_eta.clear();
	  jet_phi.clear();
	  jet_ch.clear();

  if(jets.isValid()){
     // get the number of jets in the event
     numjet=(*jets).size();
        for (reco::PFJetCollection::const_iterator itjet=jets->begin(); itjet!=jets->end(); ++itjet){

      	    jet_e.push_back(itjet->energy());
      	    jet_pt.push_back(itjet->pt());
      	    jet_px.push_back(itjet->px());
      	    jet_py.push_back(itjet->py());
      	    jet_pz.push_back(itjet->pz());
      	    jet_eta.push_back(itjet->eta());
      	    jet_phi.push_back(itjet->phi());
      	    jet_ch.push_back(itjet->charge());
        }
  }
}

void
PhysicsObjectsInfo::analyzeMets(const edm::Event& iEvent, const edm::Handle<reco::PFMETCollection> &mets)
{
	  nummet = 0;
	  met_e.clear();
	  met_pt.clear();
	  met_px.clear();
	  met_py.clear();
	  met_phi.clear();
	  met_ch.clear();

  if(mets.isValid()){
     // get the number of mets in the event
     nummet=(*mets).size();
        for (reco::PFMETCollection::const_iterator itmet=mets->begin(); itmet!=mets->end(); ++itmet){

      	    met_e.push_back(itmet->energy());
      	    met_pt.push_back(itmet->pt());
      	    met_px.push_back(itmet->px());
      	    met_py.push_back(itmet->py());
      	    met_phi.push_back(itmet->phi());
      	    met_ch.push_back(itmet->charge());
        }
  }
}

void
PhysicsObjectsInfo::analyzeMuons(const edm::Event& iEvent, const edm::Handle<reco::MuonCollection> &muons)
{
	  nummuon = 0;
	  muon_e.clear();
	  muon_pt.clear();
	  muon_px.clear();
	  muon_py.clear();
	  muon_pz.clear();
	  muon_eta.clear();
	  muon_phi.clear();
	  muon_ch.clear();

  if(muons.isValid()){
     // get the number of muons in the event
     nummuon=(*muons).size();
        for (reco::MuonCollection::const_iterator itmuon=muons->begin(); itmuon!=muons->end(); ++itmuon){

      	    muon_e.push_back(itmuon->energy());
      	    muon_pt.push_back(itmuon->pt());
      	    muon_px.push_back(itmuon->px());
      	    muon_py.push_back(itmuon->py());
      	    muon_pz.push_back(itmuon->pz());
      	    muon_eta.push_back(itmuon->eta());
      	    muon_phi.push_back(itmuon->phi());
      	    muon_ch.push_back(itmuon->charge());
        }
  }
}

void
PhysicsObjectsInfo::analyzePhotons(const edm::Event& iEvent, const edm::Handle<reco::PhotonCollection> &photons)
{
	  numphoton = 0;
	  photon_e.clear();
	  photon_pt.clear();
	  photon_px.clear();
	  photon_py.clear();
	  photon_pz.clear();
	  photon_eta.clear();
	  photon_phi.clear();
	  photon_ch.clear();

  if(photons.isValid()){
     // get the number of photons in the event
     numphoton=(*photons).size();
        for (reco::PhotonCollection::const_iterator itphoton=photons->begin(); itphoton!=photons->end(); ++itphoton){

      	    photon_e.push_back(itphoton->energy());
      	    photon_pt.push_back(itphoton->pt());
      	    photon_px.push_back(itphoton->px());
      	    photon_py.push_back(itphoton->py());
      	    photon_pz.push_back(itphoton->pz());
      	    photon_eta.push_back(itphoton->eta());
      	    photon_phi.push_back(itphoton->phi());
      	    photon_ch.push_back(itphoton->charge());
        }
  }
}

void
PhysicsObjectsInfo::analyzeTaus(const edm::Event& iEvent, const edm::Handle<reco::PFTauCollection> &taus)
{
	  numtau = 0;
	  tau_e.clear();
	  tau_pt.clear();
	  tau_px.clear();
	  tau_py.clear();
	  tau_pz.clear();
	  tau_eta.clear();
	  tau_phi.clear();
	  tau_ch.clear();

  if(taus.isValid()){
     // get the number of taus in the event
     numtau=taus->size();
        for (reco::PFTauCollection::const_iterator itTau=taus->begin(); itTau!=taus->end(); ++itTau){

    	        tau_e.push_back(itTau->energy());
    	        tau_pt.push_back(itTau->pt());
    	        tau_px.push_back(itTau->px());
    	        tau_py.push_back(itTau->py());
    	        tau_pz.push_back(itTau->pz());
    	        tau_eta.push_back(itTau->eta());
    	        tau_phi.push_back(itTau->phi());
    	        tau_ch.push_back(itTau->charge());
        }
  }
}

void
PhysicsObjectsInfo::analyzeTriggObject(const edm::Event& iEvent, const edm::Handle<trigger::TriggerEvent> &trigEvent, const edm::InputTag &trigEventTag_)
{
	  numtrigobj = 0;
	  trigobj_e.clear();
	  trigobj_pt.clear();
	  trigobj_px.clear();
	  trigobj_py.clear();
	  trigobj_pz.clear();
	  trigobj_eta.clear();
	  trigobj_phi.clear();

    trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterName_,"",trigEventTag_.process()));
    if(filterIndex<trigEvent->sizeFilters()){
        const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex);
        const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());

        //now loop of the trigger objects passing filter
        for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){
              const trigger::TriggerObject trigobj = trigObjColl[*keyIt];

              //do what you want with the trigger objects, you have
              //eta,phi,pt,mass,p,px,py,pz,et,energy accessors
        	    trigobj_e.push_back(trigobj.energy());
        	    trigobj_pt.push_back(trigobj.pt());
        	    trigobj_px.push_back(trigobj.px());
        	    trigobj_py.push_back(trigobj.py());
        	    trigobj_pz.push_back(trigobj.pz());
        	    trigobj_eta.push_back(trigobj.eta());
        	    trigobj_phi.push_back(trigobj.phi());

        	    numtrigobj=numtrigobj+1;
        }
  }//end filter size check
}

// ------------ method called once each job just before starting event loop  ------------
void
PhysicsObjectsInfo::beginJob()
{

mfile = new TFile("ObjectInfoNtuple.root","RECREATE");
mtree = new TTree("mtree","Objects information");

  mtree->Branch("numberelectron",&numelectron);
  mtree->Branch("electron_e",&electron_e);
  mtree->Branch("electron_pt",&electron_pt);
  mtree->Branch("electron_px",&electron_px);
  mtree->Branch("electron_py",&electron_py);
  mtree->Branch("electron_pz",&electron_pz);
  mtree->Branch("electron_eta",&electron_eta);
  mtree->Branch("electron_phi",&electron_phi);
  mtree->Branch("electron_ch",&electron_ch);

  mtree->Branch("numberjet",&numjet);
  mtree->Branch("jet_e",&jet_e);
  mtree->Branch("jet_pt",&jet_pt);
  mtree->Branch("jet_px",&jet_px);
  mtree->Branch("jet_py",&jet_py);
  mtree->Branch("jet_pz",&jet_pz);
  mtree->Branch("jet_eta",&jet_eta);
  mtree->Branch("jet_phi",&jet_phi);
  mtree->Branch("jet_ch",&jet_ch);

  mtree->Branch("numbermet",&nummet);
  mtree->Branch("met_e",&met_e);
  mtree->Branch("met_pt",&met_pt);
  mtree->Branch("met_px",&met_px);
  mtree->Branch("met_py",&met_py);
  mtree->Branch("met_phi",&met_phi);
  mtree->Branch("met_ch",&met_ch);

  mtree->Branch("numbermuon",&nummuon);
  mtree->Branch("muon_e",&muon_e);
  mtree->Branch("muon_pt",&muon_pt);
  mtree->Branch("muon_px",&muon_px);
  mtree->Branch("muon_py",&muon_py);
  mtree->Branch("muon_pz",&muon_pz);
  mtree->Branch("muon_eta",&muon_eta);
  mtree->Branch("muon_phi",&muon_phi);
  mtree->Branch("muon_ch",&muon_ch);

  mtree->Branch("numberphoton",&numphoton);
  mtree->Branch("photon_e",&photon_e);
  mtree->Branch("photon_pt",&photon_pt);
  mtree->Branch("photon_px",&photon_px);
  mtree->Branch("photon_py",&photon_py);
  mtree->Branch("photon_pz",&photon_pz);
  mtree->Branch("photon_eta",&photon_eta);
  mtree->Branch("photon_phi",&photon_phi);
  mtree->Branch("photon_ch",&photon_ch);

  mtree->Branch("numbertau",&numtau);
  mtree->Branch("tau_e",&tau_e);
  mtree->Branch("tau_pt",&tau_pt);
  mtree->Branch("tau_px",&tau_px);
  mtree->Branch("tau_py",&tau_py);
  mtree->Branch("tau_pz",&tau_pz);
  mtree->Branch("tau_eta",&tau_eta);
  mtree->Branch("tau_phi",&tau_phi);
  mtree->Branch("tau_ch",&tau_ch);

  mtree->Branch("numbertrigobj",&numtrigobj);
  mtree->Branch("trigobj_e",&trigobj_e);
  mtree->Branch("trigobj_pt",&trigobj_pt);
  mtree->Branch("trigobj_px",&trigobj_px);
  mtree->Branch("trigobj_py",&trigobj_py);
  mtree->Branch("trigobj_pz",&trigobj_pz);
  mtree->Branch("trigobj_eta",&trigobj_eta);
  mtree->Branch("trigobj_phi",&trigobj_phi);
}

// ------------ method called once each job just after ending the event loop  ------------
void
PhysicsObjectsInfo::endJob()
{
  //In order to save the file:
  mfile->Write();
}

// ------------ method called when starting to processes a run  ------------
void
PhysicsObjectsInfo::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
PhysicsObjectsInfo::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void
PhysicsObjectsInfo::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
PhysicsObjectsInfo::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhysicsObjectsInfo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhysicsObjectsInfo);
