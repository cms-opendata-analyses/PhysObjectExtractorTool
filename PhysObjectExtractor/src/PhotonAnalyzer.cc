// -*- C++ -*-
//
// Package:    Photon/PhotonAnalyzer
// Class:      PhotonAnalyzer
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

//class to extract photon information
#include "DataFormats/PatCandidates/interface/Photon.h"
 
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



class PhotonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PhotonAnalyzer(const edm::ParameterSet&);
      ~PhotonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<pat::PhotonCollection> photonToken_;

      // ----------member data ---------------------------
      
     
   TTree *mtree;
   int numphoton; //number of photons in the event
   std::vector<float> photon_e;
   std::vector<float> photon_pt;
   std::vector<float> photon_px;
   std::vector<float> photon_py;
   std::vector<float> photon_pz;
   std::vector<float> photon_eta;
   std::vector<float> photon_phi;
   std::vector<float> photon_ch;
   std::vector<float> photon_chIso;
   std::vector<float> photon_nhIso;
   std::vector<float> photon_phIso;
   std::vector<bool>  photon_isLoose;
   std::vector<bool>  photon_isMedium;
   std::vector<bool>  photon_isTight;
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
PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& iConfig): 
 photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons")))
   
{
   //now do what ever initialization is needed
   
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
  
   mtree->Branch("numberphoton",&numphoton);
   mtree->GetBranch("numberphoton")->SetTitle("number of photons");
   mtree->Branch("photon_e",&photon_e);
   mtree->GetBranch("photon_e")->SetTitle("photon energy");
   mtree->Branch("photon_pt",&photon_pt);
   mtree->GetBranch("photon_pt")->SetTitle("photon transverse momentum");
   mtree->Branch("photon_px",&photon_px);
   mtree->GetBranch("photon_px")->SetTitle("photon momentum x-component");
   mtree->Branch("photon_py",&photon_py);
   mtree->GetBranch("photon_py")->SetTitle("photon momentum y-component");
   mtree->Branch("photon_pz",&photon_pz);
   mtree->GetBranch("photon_pz")->SetTitle("photon momentum z-component");
   mtree->Branch("photon_eta",&photon_eta);
   mtree->GetBranch("photon_eta")->SetTitle("photon pseudorapidity");
   mtree->Branch("photon_phi",&photon_phi);
   mtree->GetBranch("photon_phi")->SetTitle("photon polar angle");
   mtree->Branch("photon_ch",&photon_ch);
   mtree->GetBranch("photon_ch")->SetTitle("photon charge");
   mtree->Branch("photon_chIso",&photon_chIso);
   mtree->GetBranch("photon_chIso")->SetTitle("corrected photon charged hadron isolation");
   mtree->Branch("photon_nhIso",&photon_nhIso);
   mtree->GetBranch("photon_nhIso")->SetTitle("corrected photon neutral hadron isolation");  
   mtree->Branch("photon_phIso",&photon_phIso);
   mtree->GetBranch("photon_phIso")->SetTitle("corrected photon isolation from other photons");
   mtree->Branch("photon_isLoose",&photon_isLoose);
   mtree->GetBranch("photon_isLoose")->SetTitle("photon tagged loose");
   mtree->Branch("photon_isMedium",&photon_isMedium);
   mtree->GetBranch("photon_isMedium")->SetTitle("photon tagged medium");
   mtree->Branch("photon_isTight",&photon_isTight);
   mtree->GetBranch("photon_isTight")->SetTitle("photon tagged tight");
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

   Handle<pat::PhotonCollection> photons;
   iEvent.getByToken(photonToken_, photons);

   numphoton = 0;
   photon_e.clear();
   photon_pt.clear();
   photon_px.clear();
   photon_py.clear();
   photon_pz.clear();
   photon_eta.clear();
   photon_phi.clear();
   photon_ch.clear();
   photon_chIso.clear();
   photon_nhIso.clear();
   photon_phIso.clear();
   photon_isLoose.clear();
   photon_isMedium.clear();
   photon_isTight.clear();

     for (const pat::Photon &pho : *photons)
    {
      photon_e.push_back(pho.energy());
	   photon_pt.push_back(pho.pt());
	   photon_px.push_back(pho.px());
	   photon_py.push_back(pho.py());
	   photon_pz.push_back(pho.pz());
	   photon_eta.push_back(pho.eta());
	   photon_phi.push_back(pho.phi());
	   photon_ch.push_back(pho.charge());
      photon_chIso.push_back(pho.chargedHadronIso());
      photon_nhIso.push_back(pho.neutralHadronIso());
      photon_phIso.push_back(pho.photonIso());
	   photon_isLoose.push_back(pho.photonID("cutBasedPhotonID-PHYS14-PU20bx25-V2p1-standalone-loose"));
      photon_isMedium.push_back(pho.photonID("cutBasedPhotonID-PHYS14-PU20bx25-V2p1-standalone-medium"));
      photon_isTight.push_back(pho.photonID("cutBasedPhotonID-PHYS14-PU20bx25-V2p1-standalone-tight"));

	   numphoton++;
    } 

  mtree->Fill();
  return;      
          
 
}


// ------------ method called once each job just before starting event loop  ------------
void
PhotonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
PhotonAnalyzer::endJob()
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
