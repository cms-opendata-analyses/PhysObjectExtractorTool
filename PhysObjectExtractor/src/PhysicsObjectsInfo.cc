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

//classes to extract Muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h" 

//classes to extract Photon information
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"

//class to save the histograms
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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
	
	//declare the variables for save data
	TFile *mfile;
	TTree *mtree;
	
//declare histograms and variables that will go into the root files
// ----------electron member data ---------------------------
	int numelectron; //number of electrons in the event
	TH1D *electronhisto;
	TH1D *electronhist_e;
	TH1D *electronhist_pt;
	TH1D *electronhist_px;
	TH1D *electronhist_py;
	TH1D *electronhist_pz;
	TH1D *electronhist_eta;
	TH1D *electronhist_phi;
	TH1D *electronhist_ch;
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
	TH1D *jethisto;
	TH1D *jethist_e;
	TH1D *jethist_pt;
	TH1D *jethist_px;
	TH1D *jethist_py;
	TH1D *jethist_pz;
	TH1D *jethist_eta;
	TH1D *jethist_phi;
	TH1D *jethist_ch;
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
	TH1D *methisto;
	TH1D *methist_e;
	TH1D *methist_pt;
	TH1D *methist_px;
	TH1D *methist_py;
	TH1D *methist_phi;
	TH1D *methist_ch;
	std::vector<float> met_e;
  	std::vector<float> met_pt;
  	std::vector<float> met_px;
  	std::vector<float> met_py;
  	std::vector<float> met_phi;
  	std::vector<float> met_ch;	

// ----------muon member data ---------------------------
	int nummuon; //number of muons in the event
	TH1D *muonhisto;
	TH1D *muonhist_e;
	TH1D *muonhist_pt;
	TH1D *muonhist_px;
	TH1D *muonhist_py;
	TH1D *muonhist_pz;
	TH1D *muonhist_eta;
	TH1D *muonhist_phi;
	TH1D *muonhist_ch;
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
	TH1D *photonhisto;
	TH1D *photonhist_e;
	TH1D *photonhist_pt;
	TH1D *photonhist_px;
	TH1D *photonhist_py;
	TH1D *photonhist_pz;
	TH1D *photonhist_eta;
	TH1D *photonhist_phi;
	TH1D *photonhist_ch;
	std::vector<float> photon_e;
  	std::vector<float> photon_pt;
  	std::vector<float> photon_px;
  	std::vector<float> photon_py;
  	std::vector<float> photon_pz;
  	std::vector<float> photon_eta;
  	std::vector<float> photon_phi;
  	std::vector<float> photon_ch;
	
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
	//to access the TFileService object in a framework module
	edm::Service<TFileService> fs;

//create electron histograms using the function template "make"	
	electronhist_e = fs->make <TH1D>("Electron energy", "Electron energy", 100, 0, 5000);
	electronhist_pt = fs->make <TH1D>("Electron pt", "Electron pt ", 100,0,5000 );
	electronhist_px = fs->make <TH1D>("Electron px", "Electron px ", 100, 0, 5000 );
	electronhist_py = fs->make <TH1D>("Electron py", "Electron py ", 100, 0, 5000 );
	electronhist_pz = fs->make <TH1D>("Electron pz", "Electron pz ", 100, 0, 5000 );
	electronhist_eta = fs->make <TH1D>("Electron eta", "Electron eta ", 100, 0, 5000 );
	electronhist_phi = fs->make <TH1D>("Electron phi", "Electron phi ", 100, 0, 5000 );
	electronhist_ch =  fs->make <TH1D>("Electron ch", "Electron ch ", 100,0,5000 );
	electronhisto = fs->make <TH1D>("Electron histogram", "Electron histo", 100, 0, 5000);
//create jet histograms
	jethist_e = fs->make <TH1D>("Jet energy", "Jet energy", 100, 0, 5000);
	jethist_pt = fs->make <TH1D>("Jet pt", "Jet pt ", 100,0,5000 );
	jethist_px = fs->make <TH1D>("Jet px", "Jet px ", 100, 0, 5000 );
	jethist_py = fs->make <TH1D>("Jet py", "Jet py ", 100, 0, 5000 );
	jethist_pz = fs->make <TH1D>("Jet pz", "Jet pz ", 100, 0, 5000 );
	jethist_eta = fs->make <TH1D>("Jet eta", "Jet eta ", 100, 0, 5000 );
	jethist_phi = fs->make <TH1D>("Jet phi", "Jet phi ", 100, 0, 5000 );
	jethist_ch =  fs->make <TH1D>("Jet ch", "Jet ch ", 100,0,5000 );
	jethisto = fs->make <TH1D>("Jet histogram", "Jet histo", 100, 0, 5000);
//create met histograms
	methist_e = fs->make <TH1D>("Met energy", "Met energy", 100, 0, 5000);
	methist_pt = fs->make <TH1D>("Met pt", "Met pt ", 100,0,5000 );
	methist_px = fs->make <TH1D>("Met px", "Met px ", 100, 0, 5000 );
	methist_py = fs->make <TH1D>("Met py", "Met py ", 100, 0, 5000 );
	methist_phi = fs->make <TH1D>("Met phi", "Met phi ", 100, 0, 5000 );
	methist_ch =  fs->make <TH1D>("Met ch", "Met ch ", 100,0,5000 );
	methisto = fs->make <TH1D>("Met histogram", "Met histo", 100, 0, 5000);
//create muon histograms
	muonhist_e = fs->make <TH1D>("Muon energy", "Muon energy", 100, 0, 5000);
	muonhist_pt = fs->make <TH1D>("Muon pt", "Muon pt ", 100,0,5000 );
	muonhist_px = fs->make <TH1D>("Muon px", "Muon px ", 100, 0, 5000 );
	muonhist_py = fs->make <TH1D>("Muon py", "Muon py ", 100, 0, 5000 );
	muonhist_pz = fs->make <TH1D>("Muon pz", "Muon pz ", 100, 0, 5000 );
	muonhist_eta = fs->make <TH1D>("Muon eta", "Muon eta ", 100, 0, 5000 );
	muonhist_phi = fs->make <TH1D>("Muon phi", "Muon phi ", 100, 0, 5000 );
	muonhist_ch =  fs->make <TH1D>("Muon ch", "Muon ch ", 100,0,5000 );
	muonhisto = fs->make <TH1D>("Muon histogram", "Muon histo", 100, 0, 5000);
//create photon histograms
	photonhist_e = fs->make <TH1D>("Photon energy", "Photon energy", 100, 0, 5000);
	photonhist_pt = fs->make <TH1D>("Photon pt", "Photon pt ", 100,0,5000 );
	photonhist_px = fs->make <TH1D>("Photon px", "Photon px ", 100, 0, 5000 );
	photonhist_py = fs->make <TH1D>("Photon py", "Photon py ", 100, 0, 5000 );
	photonhist_pz = fs->make <TH1D>("Photon pz", "Photon pz ", 100, 0, 5000 );
	photonhist_eta = fs->make <TH1D>("Photon eta", "Photon eta ", 100, 0, 5000 );
	photonhist_phi = fs->make <TH1D>("Photon phi", "Photon phi ", 100, 0, 5000 );
	photonhist_ch =  fs->make <TH1D>("Photon ch", "Photon ch ", 100,0,5000 );
	photonhisto = fs->make <TH1D>("Photon histogram", "Photon histo", 100, 0, 5000);	

//this takes the type of object (input tag) specified in the configuration python file, corresponding to each container.
//https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable	
	electronInput = iConfig.getParameter<edm::InputTag>("ElectronInputCollection");
	jetInput = iConfig.getParameter<edm::InputTag>("JetInputCollection");
	metInput = iConfig.getParameter<edm::InputTag>("MetInputCollection");
	muonInput = iConfig.getParameter<edm::InputTag>("MuonInputCollection");
	photonInput = iConfig.getParameter<edm::InputTag>("PhotonInputCollection");
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
//Then, objects of the input tag type placed in the configuration file, in this case "electrons", 
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

//Now, the information is stored.
   mtree->Fill();
   return;

}

//************************************************************************
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
     elechisto->Fill(electrons->size());
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
            //fill each histogram
	    electronhist_e->Fill(itElec->energy());
	    electronhist_pt->Fill(itElec->pt());
	    electronhist_px->Fill(itElec->px());
	    electronhist_py->Fill(itElec->py());
	    electronhist_pz->Fill(itElec->pz());
	    electronhist_eta->Fill(itElec->eta());
	    electronhist_phi->Fill(itElec->phi());
	    electronhist_ch->Fill(itElec->charge());
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
     jethisto->Fill(jets->size());
        for (reco::PFJetCollection::const_iterator itjet=jets->begin(); itjet!=jets->end(); ++itjet){

	    jet_e.push_back(itjet->energy());
	    jet_pt.push_back(itjet->pt());
	    jet_px.push_back(itjet->px());
	    jet_py.push_back(itjet->py());
	    jet_pz.push_back(itjet->pz());
	    jet_eta.push_back(itjet->eta());
	    jet_phi.push_back(itjet->phi());
	    jet_ch.push_back(itjet->charge());

	    jethist_e->Fill(itjet->energy());
	    jethist_pt->Fill(itjet->pt());
	    jethist_px->Fill(itjet->px());
	    jethist_py->Fill(itjet->py());
	    jethist_pz->Fill(itjet->pz());
	    jethist_eta->Fill(itjet->eta());
	    jethist_phi->Fill(itjet->phi());
	    jethist_ch->Fill(itjet->charge());

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
     methisto->Fill(mets->size());
        for (reco::PFMETCollection::const_iterator itmet=mets->begin(); itmet!=mets->end(); ++itmet){

	    met_e.push_back(itmet->energy());
	    met_pt.push_back(itmet->pt());
	    met_px.push_back(itmet->px());
	    met_py.push_back(itmet->py());
	    met_phi.push_back(itmet->phi());
	    met_ch.push_back(itmet->charge());

	    methist_e->Fill(itmet->energy());
	    methist_pt->Fill(itmet->pt());
	    methist_px->Fill(itmet->px());
	    methist_py->Fill(itmet->py());
	    methist_phi->Fill(itmet->phi());
	    methist_ch->Fill(itmet->charge());

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
     muonhisto->Fill(muons->size());
        for (reco::MuonCollection::const_iterator itmuon=muons->begin(); itmuon!=muons->end(); ++itmuon){

	    muon_e.push_back(itmuon->energy());
	    muon_pt.push_back(itmuon->pt());
	    muon_px.push_back(itmuon->px());
	    muon_py.push_back(itmuon->py());
	    muon_pz.push_back(itmuon->pz());
	    muon_eta.push_back(itmuon->eta());
	    muon_phi.push_back(itmuon->phi());
	    muon_ch.push_back(itmuon->charge());

	    muonhist_e->Fill(itmuon->energy());
	    muonhist_pt->Fill(itmuon->pt());
	    muonhist_px->Fill(itmuon->px());
	    muonhist_py->Fill(itmuon->py());
	    muonhist_pz->Fill(itmuon->pz());
	    muonhist_eta->Fill(itmuon->eta());
	    muonhist_phi->Fill(itmuon->phi());
	    muonhist_ch->Fill(itmuon->charge());

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
     photonhisto->Fill(photons->size());
        for (reco::PhotonCollection::const_iterator itphoton=photons->begin(); itphoton!=photons->end(); ++itphoton){

	    photon_e.push_back(itphoton->energy());
	    photon_pt.push_back(itphoton->pt());
	    photon_px.push_back(itphoton->px());
	    photon_py.push_back(itphoton->py());
	    photon_pz.push_back(itphoton->pz());
	    photon_eta.push_back(itphoton->eta());
	    photon_phi.push_back(itphoton->phi());
	    photon_ch.push_back(itphoton->charge());

	    photonhist_e->Fill(itphoton->energy());
	    photonhist_pt->Fill(itphoton->pt());
	    photonhist_px->Fill(itphoton->px());
	    photonhist_py->Fill(itphoton->py());
	    photonhist_pz->Fill(itphoton->pz());
	    photonhist_eta->Fill(itphoton->eta());
	    photonhist_phi->Fill(itphoton->phi());
	    photonhist_ch->Fill(itphoton->charge());

        }
  }
}
//*************************************************************************


// ------------ method called once each job just before starting event loop  ------------
void
PhysicsObjectsInfo::beginJob()
{

mfile = new TFile("ObjectInfoNtuple.root","RECREATE");
mtree = new TTree("mtree","Objects information");
	
  mtree->Branch("electron_e",&electron_e);
  mtree->Branch("electron_pt",&electron_pt);
  mtree->Branch("electron_px",&electron_px);
  mtree->Branch("electron_py",&electron_py);
  mtree->Branch("electron_pz",&electron_pz);
  mtree->Branch("electron_eta",&electron_eta);
  mtree->Branch("electron_phi",&electron_phi);
  mtree->Branch("electron_ch",&electron_ch);
	
  mtree->Branch("jet_e",&jet_e);
  mtree->Branch("jet_pt",&jet_pt);
  mtree->Branch("jet_px",&jet_px);
  mtree->Branch("jet_py",&jet_py);
  mtree->Branch("jet_pz",&jet_pz);
  mtree->Branch("jet_eta",&jet_eta);
  mtree->Branch("jet_phi",&jet_phi);
  mtree->Branch("jet_ch",&jet_ch);

  mtree->Branch("met_e",&met_e);
  mtree->Branch("met_pt",&met_pt);
  mtree->Branch("met_px",&met_px);
  mtree->Branch("met_py",&met_py);
  mtree->Branch("met_phi",&met_phi);
  mtree->Branch("met_ch",&met_ch);	
	
  mtree->Branch("muon_e",&muon_e);
  mtree->Branch("muon_pt",&muon_pt);
  mtree->Branch("muon_px",&muon_px);
  mtree->Branch("muon_py",&muon_py);
  mtree->Branch("muon_pz",&muon_pz);
  mtree->Branch("muon_eta",&muon_eta);
  mtree->Branch("muon_phi",&muon_phi);
  mtree->Branch("muon_ch",&muon_ch);	
	
  mtree->Branch("photon_e",&photon_e);
  mtree->Branch("photon_pt",&photon_pt);
  mtree->Branch("photon_px",&photon_px);
  mtree->Branch("photon_py",&photon_py);
  mtree->Branch("photon_pz",&photon_pz);
  mtree->Branch("photon_eta",&photon_eta);
  mtree->Branch("photon_phi",&photon_phi);
  mtree->Branch("photon_ch",&photon_ch);	
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
