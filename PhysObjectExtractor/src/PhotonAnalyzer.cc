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

//classes to extract Photon information
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

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
//declare a function to do the photon analysis
      void analyzePhotons(const edm::Event& iEvent, const edm::Handle<reco::PhotonCollection> &photons);


//declare the input tag for PhotonCollection
      edm::InputTag photonInput;

	  // ----------member data ---------------------------

	int numphoton; //number of photons in the event

	TFile *mfile;
	TTree *mtree;

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

PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& iConfig)

{
//now do what ever initialization is needed
	photonInput = iConfig.getParameter<edm::InputTag>("InputCollection");
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

   analyzePhotons(iEvent,myphotons);

   mtree->Fill();
   return;
}

void
PhotonAnalyzer::analyzePhotons(const edm::Event& iEvent, const edm::Handle<reco::PhotonCollection> &photons)
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

// ------------ method called once each job just before starting event loop  ------------
void
PhotonAnalyzer::beginJob()
{

mfile = new TFile("PhotonInfo.root","RECREATE");
mtree = new TTree("mtree","Photon information");

  mtree->Branch("numberphoton",&numphoton);
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
PhotonAnalyzer::endJob()
{
  mfile->Write();
}

// ------------ method called when starting to processes a run  ------------
void
PhotonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
PhotonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void
PhotonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

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
