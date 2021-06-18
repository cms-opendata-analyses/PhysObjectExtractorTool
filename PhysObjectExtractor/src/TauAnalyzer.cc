// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      TauAnalyzer
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

//classes to extract tau information
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class TauAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TauAnalyzer(const edm::ParameterSet&);
      ~TauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      
      //declare the input tag for PFTauCollection
      edm::InputTag tauInput;
      
      // ----------member data ---------------------------
      
      
      TTree *mtree;
      int numtau; //number of taus in the event
      std::vector<float> tau_e;
      std::vector<float> tau_pt;
      std::vector<float> tau_px;
      std::vector<float> tau_py;
      std::vector<float> tau_pz;
      std::vector<float> tau_eta;
      std::vector<float> tau_phi;
      std::vector<float> tau_ch;
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

TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
	tauInput = iConfig.getParameter<edm::InputTag>("InputCollection");
	edm::Service<TFileService> fs;
	mtree = fs->make<TTree>("Events", "Events");
	
	
	mtree->Branch("numbertau",&numtau);
	mtree->Branch("tau_e",&tau_e);
	mtree->Branch("tau_pt",&tau_pt);
	mtree->Branch("tau_px",&tau_px);
	mtree->Branch("tau_py",&tau_py);
	mtree->Branch("tau_pz",&tau_pz);
	mtree->Branch("tau_eta",&tau_eta);
	mtree->Branch("tau_phi",&tau_phi);
	mtree->Branch("tau_ch",&tau_ch);
	
}

TauAnalyzer::~TauAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<reco::PFTauCollection> mytaus;
   iEvent.getByLabel(tauInput, mytaus);
   
   numtau = 0;
   tau_e.clear();
   tau_pt.clear();
   tau_px.clear();
   tau_py.clear();
   tau_pz.clear();
   tau_eta.clear();
   tau_phi.clear();
   tau_ch.clear();

  if(mytaus.isValid()){
     // get the number of taus in the event
     numtau=mytaus->size();
        for (reco::PFTauCollection::const_iterator itTau=mytaus->begin(); itTau!=mytaus->end(); ++itTau){

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
  
  mtree->Fill();
  return;
   
}

// ------------ method called once each job just before starting event loop  ------------
void
TauAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
TauAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
TauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
TauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
TauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
