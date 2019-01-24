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
#include "DataFormats/JetReco/interface/PFJetFwd.h"

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
      void analyzeMuons(const edm::Event& iEvent, const edm::Handle<reco::GsfElectronCollection> &electrons);	

//declare a function to do the jet analysis
      void analyzeJets(const edm::Event& iEvent, const edm::Handle<reco::PFJetCollection> &jets);
	
//declare a function to do the met analysis
      void analyzeMets(const edm::Event& iEvent, const edm::Handle<reco::PFMETCollection> &mets);	

//declare a function to do the muon analysis
      void analyzeMuons(const edm::Event& iEvent, const edm::Handle<reco::MuonCollection> &muons);

//declare a function to do the photon analysis
      void analyzePhotons(const edm::Event& iEvent, const edm::Handle<reco::PhotonCollection> &photons);
	
//se declara el input tag de tipo GsfElectronCollection         
      edm::InputTag electronInput;
//se declara el input tag de tipo PFJetCollection         
      edm::InputTag jetInput;
//se declara el input tag de tipo PFMETCollection         
      edm::InputTag metInput;
//se declara el input tag de tipo MuonCollection         
      edm::InputTag muonInput;
//se declara el input tag de tipo PhotonCollection         
      edm::InputTag photonInput;	

// ----------electron member data ---------------------------
	int numelectron; //number of electrons in the event
	TH1D *elechisto;
	TH1D *electronhist_e;
	TH1D *electronhist_pt;
	TH1D *electronhist_px;
	TH1D *electronhist_py;
	TH1D *electronhist_pz;
	TH1D *electronhist_eta;
	TH1D *electronhist_phi;
	TH1D *electronhist_ch;
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
	//TFile *mfile;
	//TTree *mtree;
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
	//TH1D *methist_pz;
	//TH1D *methist_eta;
	TH1D *methist_phi;
	TH1D *methist_ch;
	//TFile *mfile;
	//TTree *mtree;
	std::vector<float> met_e;
  	std::vector<float> met_pt;
  	std::vector<float> met_px;
  	std::vector<float> met_py;
  	//std::vector<float> met_pz;
  	//std::vector<float> met_eta;
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
	TH1D *hist_ch;
	//TFile *mfile;
	//TTree *mtree;
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
	//TFile *mfile;
	//TTree *mtree;
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
	edm::Service<TFileService> fs;

// se crean los electron histogramas
	hist_e = fs->make <TH1D>("hist_energy", "Electron energy", 100, 0, 5000);
	hist_pt = fs->make <TH1D>("hist_pt", "Electron pt ", 100,0,5000 );
	hist_px = fs->make <TH1D>("hist_px", "Electron px ", 100, 0, 5000 );
	hist_py = fs->make <TH1D>("hist_py", "Electron py ", 100, 0, 5000 );
	hist_pz = fs->make <TH1D>("hist_pz", "Electron pz ", 100, 0, 5000 );
	hist_eta = fs->make <TH1D>("hist_eta", "Electron eta ", 100, 0, 5000 );
	hist_phi = fs->make <TH1D>("hist_phi", "Electron phi ", 100, 0, 5000 );
	hist_ch =  fs->make <TH1D>("hist_ch", "Electron ch ", 100,0,5000 );
	elechisto = fs->make <TH1D>("elechisto", "Electron histo", 100, 0, 5000);
// se crean los jet histogramas
	hist_e = fs->make <TH1D>("hist_energy", "Jet energy", 100, 0, 5000);
	hist_pt = fs->make <TH1D>("hist_pt", "Jet pt ", 100,0,5000 );
	hist_px = fs->make <TH1D>("hist_px", "Jet px ", 100, 0, 5000 );
	hist_py = fs->make <TH1D>("hist_py", "Jet py ", 100, 0, 5000 );
	hist_pz = fs->make <TH1D>("hist_pz", "Jet pz ", 100, 0, 5000 );
	hist_eta = fs->make <TH1D>("hist_eta", "Jet eta ", 100, 0, 5000 );
	hist_phi = fs->make <TH1D>("hist_phi", "Jet phi ", 100, 0, 5000 );
	hist_ch =  fs->make <TH1D>("hist_ch", "Jet ch ", 100,0,5000 );
	jethisto = fs->make <TH1D>("jethisto", "Jet histo", 100, 0, 5000);
// se crean los met histogramas
	hist_e = fs->make <TH1D>("hist_energy", "Met energy", 100, 0, 5000);
	hist_pt = fs->make <TH1D>("hist_pt", "Met pt ", 100,0,5000 );
	hist_px = fs->make <TH1D>("hist_px", "Met px ", 100, 0, 5000 );
	hist_py = fs->make <TH1D>("hist_py", "Met py ", 100, 0, 5000 );
	//hist_pz = fs->make <TH1D>("hist_pz", "Met pz ", 100, 0, 5000 );
	//hist_eta = fs->make <TH1D>("hist_eta", "Met eta ", 100, 0, 5000 );
	hist_phi = fs->make <TH1D>("hist_phi", "Met phi ", 100, 0, 5000 );
	hist_ch =  fs->make <TH1D>("hist_ch", "Met ch ", 100,0,5000 );
	methisto = fs->make <TH1D>("methisto", "Met histo", 100, 0, 5000);
// se crean los muon histogramas
	hist_e = fs->make <TH1D>("hist_energy", "Muon energy", 100, 0, 5000);
	hist_pt = fs->make <TH1D>("hist_pt", "Muon pt ", 100,0,5000 );
	hist_px = fs->make <TH1D>("hist_px", "Muon px ", 100, 0, 5000 );
	hist_py = fs->make <TH1D>("hist_py", "Muon py ", 100, 0, 5000 );
	hist_pz = fs->make <TH1D>("hist_pz", "Muon pz ", 100, 0, 5000 );
	hist_eta = fs->make <TH1D>("hist_eta", "Muon eta ", 100, 0, 5000 );
	hist_phi = fs->make <TH1D>("hist_phi", "Muon phi ", 100, 0, 5000 );
	hist_ch =  fs->make <TH1D>("hist_ch", "Muon ch ", 100,0,5000 );
	muonhisto = fs->make <TH1D>("muonhisto", "Muon histo", 100, 0, 5000);
// se crean los photon histogramas
	hist_e = fs->make <TH1D>("hist_energy", "Photon energy", 100, 0, 5000);
	hist_pt = fs->make <TH1D>("hist_pt", "Photon pt ", 100,0,5000 );
	hist_px = fs->make <TH1D>("hist_px", "Photon px ", 100, 0, 5000 );
	hist_py = fs->make <TH1D>("hist_py", "Photon py ", 100, 0, 5000 );
	hist_pz = fs->make <TH1D>("hist_pz", "Photon pz ", 100, 0, 5000 );
	hist_eta = fs->make <TH1D>("hist_eta", "Photon eta ", 100, 0, 5000 );
	hist_phi = fs->make <TH1D>("hist_phi", "Photon phi ", 100, 0, 5000 );
	hist_ch =  fs->make <TH1D>("hist_ch", "Photon ch ", 100,0,5000 );
	photonhisto = fs->make <TH1D>("photonhisto", "Photon histo", 100, 0, 5000);	

	electronInput = iConfig.getParameter<edm::InputTag>("InputCollection");
	jetInput = iConfig.getParameter<edm::InputTag>("InputCollection");
	metInput = iConfig.getParameter<edm::InputTag>("InputCollection");
	muonInput = iConfig.getParameter<edm::InputTag>("InputCollection");
	photonInput = iConfig.getParameter<edm::InputTag>("InputCollection");
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

   Handle<reco::GsfElectronCollection> myelectrons;
   iEvent.getByLabel(electronInput, myelectrons);
   analyzeElectrons(iEvent,myelectrons);
	
   Handle<reco::PFJetCollection> myjets;
   iEvent.getByLabel(jetInput, myjets);
   analyzeJets(iEvent,myjets);

   Handle<reco::PFMETCollection> mymets;
   iEvent.getByLabel(metInput, mymets);
   analyzeMets(iEvent,mymets);	

   Handle<reco::MuonCollection> mymuons;
   iEvent.getByLabel(muonInput, mymuons);
   analyzeMuons(iEvent,mymuons);
	        
	Handle<reco::PhotonCollection> myphotons;
   iEvent.getByLabel(photonInput, myphotons);
   analyzePhotons(iEvent,myphotons);

   mtree->Fill();
   return;

}

//************************************************************************

void 
PhysicsObjectsInfo::analyzeJets(const edm::Event& iEvent, const edm::Handle<reco::PFJetCollection> &jets)
{
	  numjet = 0;
	  _e.clear();
	  _pt.clear();
	  _px.clear();
	  _py.clear();
	  _pz.clear();
	  _eta.clear();
	  _phi.clear();
	  _ch.clear();

  if(jets.isValid()){
     // get the number of jets in the event
     numjet=(*jets).size();
     jethisto->Fill(jets->size());
        for (reco::PFJetCollection::const_iterator itjet=jets->begin(); itjet!=jets->end(); ++itjet){

	    _e.push_back(itjet->energy());
	    _pt.push_back(itjet->pt());
	    _px.push_back(itjet->px());
	    _py.push_back(itjet->py());
	    _pz.push_back(itjet->pz());
	    _eta.push_back(itjet->eta());
	    _phi.push_back(itjet->phi());
	    _ch.push_back(itjet->charge());

	    hist_e->Fill(itjet->energy());
	    hist_pt->Fill(itjet->pt());
	    hist_px->Fill(itjet->px());
	    hist_py->Fill(itjet->py());
	    hist_pz->Fill(itjet->pz());
	    hist_eta->Fill(itjet->eta());
	    hist_phi->Fill(itjet->phi());
	    hist_ch->Fill(itjet->charge());

        }
  }
}

//*************************************************************************


// ------------ method called once each job just before starting event loop  ------------
void
PhysicsObjectsInfo::beginJob()
{

mfile = new TFile("JetInfo.root","RECREATE");
mtree = new TTree("mtree","Jet information");

  mytree->Branch("jet_e",&_e);
  mytree->Branch("jet_pt",&_pt);
  mytree->Branch("jet_px",&_px);
  mytree->Branch("jet_py",&_py);
  mytree->Branch("jet_pz",&_pz);
  mytree->Branch("jet_eta",&_eta);
  mytree->Branch("jet_phi",&_phi);
  mytree->Branch("jet_ch",&_ch);

}

// ------------ method called once each job just after ending the event loop  ------------
void
PhysicsObjectsInfo::endJob()
{
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
