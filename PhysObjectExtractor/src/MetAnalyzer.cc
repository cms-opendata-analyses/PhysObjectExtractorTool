// -*- C++ -*-
//
// Package:    Met/MetAnalyzer
// Class:      MetAnalyzer
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

//class to extract MET information
#include "DataFormats/PatCandidates/interface/MET.h"


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



class MetAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MetAnalyzer(const edm::ParameterSet&);
      ~MetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::METCollection> rawToken_;

      // ----------member data ---------------------------
      
     TTree *mtree;
     float met_e;
     float met_pt;
     float met_px;
     float met_py;
     float met_phi;
     float met_significance;
     float met_rawe;
     float met_rawpt;
     float met_rawphi;
      
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
MetAnalyzer::MetAnalyzer(const edm::ParameterSet& iConfig): 
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  rawToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("rawmets")))
{
   //now do what ever initialization is needed
   
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
  
  mtree->Branch("met_e",&met_e,"met_e/F");	
  mtree->GetBranch("met_e")->SetTitle("Sum of transverse energy (corrected) [GeV]");
  mtree->Branch("met_pt",&met_pt,"met_pt/F");
  mtree->GetBranch("met_pt")->SetTitle("Missing transverse momentum (corrected) [GeV]");
  mtree->Branch("met_px",&met_px,"met_px/F");
  mtree->GetBranch("met_px")->SetTitle("Missing x momentum (corrected) [GeV]");
  mtree->Branch("met_py",&met_py,"met_py/F");
  mtree->GetBranch("met_py")->SetTitle("Missing y momentum (corrected) [GeV]");
  mtree->Branch("met_phi",&met_phi,"met_phi/F");
  mtree->GetBranch("met_phi")->SetTitle("Missing momentum azimuthal angle (corrected)");
  mtree->Branch("met_significance",&met_significance,"met_significance/F");
  mtree->GetBranch("met_significance")->SetTitle("Missing transverse momentum significance");

  mtree->Branch("met_rawpt",&met_rawpt,"met_rawpt/F");
  mtree->GetBranch("met_rawpt")->SetTitle("Missing transverse momentum (uncorrected) [GeV]");
  mtree->Branch("met_rawphi",&met_rawphi,"met_rawphi/F");
  mtree->GetBranch("met_phi")->SetTitle("Missing momentum azimuthal angle (uncorrected)");
  mtree->Branch("met_rawe",&met_rawe,"met_rawe/F");
  mtree->GetBranch("met_rawe")->SetTitle("Sum of transverse energy (uncorrected) [GeV]");
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

   Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);

   const pat::MET &met = mets->front();

   met_e = met.sumEt();
   met_pt = met.pt();
   met_px = met.px();
   met_py = met.py();
   met_phi = met.phi();
   met_significance = met.significance();

   Handle<pat::METCollection> rawmets;
   iEvent.getByToken(rawToken_, rawmets);

   const pat::MET &rawmet = rawmets->front();
   
   met_rawe = rawmet.sumEt();
   met_rawpt = rawmet.pt();
   met_rawphi = rawmet.phi();

   mtree->Fill();
   return;      
          
 
}


// ------------ method called once each job just before starting event loop  ------------
void
MetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MetAnalyzer::endJob()
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

