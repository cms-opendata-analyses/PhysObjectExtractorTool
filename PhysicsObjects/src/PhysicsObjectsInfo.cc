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
//declare a function to do the jet analysis
      void analyzeJets(const edm::Event& iEvent, const edm::Handle<reco::PFJetCollection> &jets);


//se declara el input tag de tipo PFJetCollection         
      edm::InputTag jetInput;

	  // ----------member data ---------------------------

	int numjet; //number of jets in the event
	TH1D *jethisto;
	TH1D *hist_e;
	TH1D *hist_pt;
	TH1D *hist_px;
	TH1D *hist_py;
	TH1D *hist_pz;
	TH1D *hist_eta;
	TH1D *hist_phi;
	TH1D *hist_ch;
	TFile *mfile;
	TTree *mtree;

	std::vector<float> _e;
  	std::vector<float> _pt;
  	std::vector<float> _px;
  	std::vector<float> _py;
  	std::vector<float> _pz;
  	std::vector<float> _eta;
  	std::vector<float> _phi;
  	std::vector<float> _ch;
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

// se crean los histogramas
	hist_e = fs->make <TH1D>("hist_energy", "Electron energy", 100, 0, 5000);
	hist_pt = fs->make <TH1D>("hist_pt", "Electron pt ", 100,0,5000 );
	hist_px = fs->make <TH1D>("hist_px", "Electron px ", 100, 0, 5000 );
	hist_py = fs->make <TH1D>("hist_py", "Electron py ", 100, 0, 5000 );
	hist_pz = fs->make <TH1D>("hist_pz", "Electron pz ", 100, 0, 5000 );
	hist_eta = fs->make <TH1D>("hist_eta", "Electron eta ", 100, 0, 5000 );
	hist_phi = fs->make <TH1D>("hist_phi", "Electron phi ", 100, 0, 5000 );
	hist_ch =  fs->make <TH1D>("hist_ch", "Electron ch ", 100,0,5000 );
	jethisto = fs->make <TH1D>("jethisto", "Jet histo", 100, 0, 5000);

	jetInput = iConfig.getParameter<edm::InputTag>("InputCollection");

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


   Handle<reco::PFJetCollection> myjets;
   iEvent.getByLabel(jetInput, myjets);

   analyzeJets(iEvent,myjets);

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
