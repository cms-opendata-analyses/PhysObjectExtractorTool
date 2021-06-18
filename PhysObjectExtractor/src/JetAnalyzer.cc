// -*- C++ -*-
//
// Package:    JetAnalyzer
// Class:      JetAnalyzer
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

//classes to extract PFJet information
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class JetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit JetAnalyzer(const edm::ParameterSet&);
      ~JetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      
      //declare the input tag for PFJetCollection
      edm::InputTag jetInput;

      // ----------member data ---------------------------
      
      
      int numjet; //number of jets in the event
      TTree *mtree;
      std::vector<float> jet_e;
      std::vector<float> jet_pt;
      std::vector<float> jet_px;
      std::vector<float> jet_py;
      std::vector<float> jet_pz;
      std::vector<float> jet_eta;
      std::vector<float> jet_phi;
      std::vector<float> jet_ch;
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

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
	jetInput = iConfig.getParameter<edm::InputTag>("InputCollection");
	edm::Service<TFileService> fs;
	mtree = fs->make<TTree>("Events", "Events");
	
	
	mtree->Branch("numberjet",&numjet);
	mtree->Branch("jet_e",&jet_e);
	mtree->Branch("jet_pt",&jet_pt);
	mtree->Branch("jet_px",&jet_px);
	mtree->Branch("jet_py",&jet_py); 
	mtree->Branch("jet_pz",&jet_pz);
	mtree->Branch("jet_eta",&jet_eta);
	mtree->Branch("jet_phi",&jet_phi);
	mtree->Branch("jet_ch",&jet_ch);
	
}

JetAnalyzer::~JetAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
JetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<reco::PFJetCollection> myjets;
   iEvent.getByLabel(jetInput, myjets);
   
   numjet = 0;
   jet_e.clear();
   jet_pt.clear();
   jet_px.clear();
   jet_py.clear();
   jet_pz.clear();
   jet_eta.clear();
   jet_phi.clear();
   jet_ch.clear();

  if(myjets.isValid()){
     // get the number of jets in the event
     numjet=myjets->size();
        for (reco::PFJetCollection::const_iterator itjet=myjets->begin(); itjet!=myjets->end(); ++itjet){

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
  
  mtree->Fill();
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void
JetAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
JetAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
JetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
JetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
JetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
JetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetAnalyzer);
