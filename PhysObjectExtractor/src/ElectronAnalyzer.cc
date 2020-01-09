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

//se declara el input tag de tipo GsfElectronCollection
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
	  numelectron = 0;
	  electron_e.clear();
	  electron_pt.clear();
	  electron_px.clear();
	  electron_py.clear();
	  electron_pz.clear();
	  electron_eta.clear();
	  electron_phi.clear();
	  electron_ch.clear();

  if(electrons.isValid()){
     // get the number of electrons in the event
     numelectron=electrons->size();
        for (reco::GsfElectronCollection::const_iterator itElec=electrons->begin(); itElec!=electrons->end(); ++itElec){

    	        electron_e.push_back(itElec->energy());
    	        electron_pt.push_back(itElec->pt());
    	        electron_px.push_back(itElec->px());
    	        electron_py.push_back(itElec->py());
    	        electron_pz.push_back(itElec->pz());
    	        electron_eta.push_back(itElec->eta());
    	        electron_phi.push_back(itElec->phi());
    	        electron_ch.push_back(itElec->charge());
        }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void
ElectronAnalyzer::beginJob()
{

mfile = new TFile("ElectronInfo.root","RECREATE");
mtree = new TTree("mtree","Electron information");

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
