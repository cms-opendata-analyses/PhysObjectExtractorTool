// -*- C++ -*-
//
// Package:    MetAnalyzer
// Class:      MetAnalyzer
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

//classes to extract PFMET information
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class MetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MetAnalyzer(const edm::ParameterSet&);
      ~MetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      
      //declare the input tag for PFMETCollection
      edm::InputTag metInput;

      // ----------member data ---------------------------
      
      TTree *mtree;
      int nummet; //number of mets in the event
      std::vector<float> met_e;
      std::vector<float> met_pt;
      std::vector<float> met_px;
      std::vector<float> met_py;
    //std::vector<float> met_pz;
    //std::vector<float> met_eta;
      std::vector<float> met_phi;
      std::vector<float> met_ch;
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

MetAnalyzer::MetAnalyzer(const edm::ParameterSet& iConfig)

{
//now do what ever initialization is needed
	metInput = iConfig.getParameter<edm::InputTag>("InputCollection");
	edm::Service<TFileService> fs;
	mtree = fs->make<TTree>("Events", "Events");
	
	
	mtree->Branch("numbermet",&nummet);
	mtree->Branch("met_e",&met_e);
	mtree->Branch("met_pt",&met_pt);
	mtree->Branch("met_px",&met_px);
	mtree->Branch("met_py",&met_py);
        //mtree->Branch("met_pz",&met_pz);
        //mtree->Branch("met_eta",&met_eta);
        mtree->Branch("met_phi",&met_phi);
        mtree->Branch("met_ch",&met_ch);
        
}

MetAnalyzer::~MetAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
MetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<reco::PFMETCollection> mymets;
   iEvent.getByLabel(metInput, mymets);
   
   nummet = 0;
   met_e.clear();
   met_pt.clear();
   met_px.clear();
   met_py.clear();
   //met_pz.clear();
   //met_eta.clear();
   met_phi.clear();
   met_ch.clear();

  if(mymets.isValid()){
     // get the number of mets in the event
     nummet=mymets->size();
        for (reco::PFMETCollection::const_iterator itmet=mymets->begin(); itmet!=mymets->end(); ++itmet){

        	    met_e.push_back(itmet->energy());
        	    met_pt.push_back(itmet->pt());
        	    met_px.push_back(itmet->px());
        	    met_py.push_back(itmet->py());
        	    //met_pz.push_back(itmet->pz());
        	    //met_eta.push_back(itmet->eta());
        	    met_phi.push_back(itmet->phi());
        	    met_ch.push_back(itmet->charge());
        }
  }
  
  mtree->Fill();
  return;
   
}

// ------------ method called once each job just before starting event loop  ------------
void
MetAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
MetAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
MetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
MetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
MetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
MetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MetAnalyzer);
