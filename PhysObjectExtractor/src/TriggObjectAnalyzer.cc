// -*- C++ -*-
//
// Package:    TriggObject/TriggObjectAnalyzer
// Class:      TriggObjectAnalyzer
//
 
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//class to extract Trigger Object information
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class TriggObjectAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TriggObjectAnalyzer(const edm::ParameterSet&);
      ~TriggObjectAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

      // ----------member data ---------------------------
      
     TTree *mtree;
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
TriggObjectAnalyzer::TriggObjectAnalyzer(const edm::ParameterSet& iConfig): 
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")))
   
{
   //now do what ever initialization is needed
   
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
  
   mtree->Branch("numbertrigobj",&numtrigobj);
   mtree->Branch("trigobj_e",&trigobj_e);
   mtree->Branch("trigobj_pt",&trigobj_pt);
   mtree->Branch("trigobj_px",&trigobj_px);
   mtree->Branch("trigobj_py",&trigobj_py);
   mtree->Branch("trigobj_pz",&trigobj_pz);
   mtree->Branch("trigobj_eta",&trigobj_eta);
   mtree->Branch("trigobj_phi",&trigobj_phi);
}


TriggObjectAnalyzer::~TriggObjectAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggObjectAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerObjects_, triggerObjects);

   numtrigobj = 0;
   trigobj_e.clear();
   trigobj_pt.clear();
   trigobj_px.clear();
   trigobj_py.clear();
   trigobj_pz.clear();
   trigobj_eta.clear();
   trigobj_phi.clear();

    for (pat::TriggerObjectStandAlone obj : *triggerObjects)
    {
        trigobj_e.push_back(obj.energy());
	    trigobj_pt.push_back(obj.pt());
	    trigobj_px.push_back(obj.px());
	    trigobj_py.push_back(obj.py());
	    trigobj_pz.push_back(obj.pz());
	    trigobj_eta.push_back(obj.eta());
	    trigobj_phi.push_back(obj.phi());

	    numtrigobj++;
    } 

  mtree->Fill();
  return;      
          
 
}


// ------------ method called once each job just before starting event loop  ------------
void
TriggObjectAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TriggObjectAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggObjectAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggObjectAnalyzer);
