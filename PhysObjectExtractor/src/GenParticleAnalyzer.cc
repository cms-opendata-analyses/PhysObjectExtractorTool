// -*- C++ -*-
//
// Package:    GenParticle/GenParticleAnalyzer
// Class:      GenParticleAnalyzer
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

//class to extract MC information
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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



class GenParticleAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GenParticleAnalyzer(const edm::ParameterSet&);
      ~GenParticleAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      std::vector<std::string>  particle;

      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;

      // ----------member data ---------------------------
      
     TTree *mtree;
      
    int numGenPart;
    std::vector<int> GenPart_status;
    std::vector<float> GenPart_pt;
    std::vector<float> GenPart_eta;
    std::vector<float> GenPart_mass;
    std::vector<int> GenPart_pdgId;
    std::vector<float> GenPart_phi;
    std::vector<float> GenPart_px;
    std::vector<float> GenPart_py;
    std::vector<float> GenPart_pz; 
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
GenParticleAnalyzer::GenParticleAnalyzer(const edm::ParameterSet& iConfig):
 particle(iConfig.getParameter<std::vector<std::string> >("input_particle")), 
 prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned")))
   
{
   //now do what ever initialization is needed
   
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
  
   mtree->Branch("numGenPart",&numGenPart);
   mtree->GetBranch("numGenPart")->SetTitle("number of generator particles");
   mtree->Branch("GenPart_pt",&GenPart_pt);
   mtree->GetBranch("GenPart_pt")->SetTitle("generator particle transverse momentum");
   mtree->Branch("GenPart_eta",&GenPart_eta);
   mtree->GetBranch("GenPart_eta")->SetTitle("generator particle pseudorapidity");
   mtree->Branch("GenPart_mass",&GenPart_mass);
   mtree->GetBranch("GenPart_mass")->SetTitle("generator particle mass");
   mtree->Branch("GenPart_pdgId",&GenPart_pdgId);
   mtree->GetBranch("GenPart_pdgId")->SetTitle("generator particle PDG id");
   mtree->Branch("GenPart_phi",&GenPart_phi);
   mtree->GetBranch("GenPart_phi")->SetTitle("generator particle azimuthal angle of momentum vector");
   mtree->Branch("GenPart_px",&GenPart_px);
   mtree->GetBranch("GenPart_px")->SetTitle("generator particle x coordinate of momentum vector");
   mtree->Branch("GenPart_py",&GenPart_py);
   mtree->GetBranch("GenPart_py")->SetTitle("generator particle y coordinate of momentum vector");
   mtree->Branch("GenPart_pz",&GenPart_pz);
   mtree->GetBranch("GenPart_pz")->SetTitle("generator particle z coordinate of momentum vector");
   mtree->Branch("GenPart_status",&GenPart_status);
   mtree->GetBranch("GenPart_status")->SetTitle("Particle status. 1=stable");

}


GenParticleAnalyzer::~GenParticleAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<edm::View<reco::GenParticle> > pruned;
   iEvent.getByToken(prunedGenToken_,pruned);

   numGenPart=0;
   GenPart_pt.clear();
   GenPart_eta.clear();
   GenPart_mass.clear();
   GenPart_pdgId.clear();
   GenPart_phi.clear();
   GenPart_px.clear();
   GenPart_py.clear();
   GenPart_pz.clear();
   GenPart_status.clear();

   unsigned int i;
   std::string s1,s2;
   std::vector<int> status_parsed;
   std::vector<int> pdgId_parsed;
   std::string delimiter = ":";

   for(i=0;i<particle.size();i++)
   {
       //get status and pgdId from configuration
       s1=particle[i].substr(0,particle[i].find(delimiter));
       s2=particle[i].substr(particle[i].find(delimiter)+1,particle[i].size());
       //parse string to int
       status_parsed.push_back(stoi(s1));
       pdgId_parsed.push_back(stoi(s2));
   }

     for (const reco::GenParticle &gpt : *pruned)
    {
        for(i=0;i<particle.size();i++)
        {
         if((status_parsed[i]==gpt.status() && pdgId_parsed[i]==gpt.pdgId())||(status_parsed[i]==0 && pdgId_parsed[i]==0))
         {
          GenPart_pt.push_back(gpt.pt());
          GenPart_eta.push_back(gpt.eta());
          GenPart_mass.push_back(gpt.mass());
          GenPart_pdgId.push_back(gpt.pdgId());
          GenPart_phi.push_back(gpt.phi());
          GenPart_status.push_back(gpt.status());
          GenPart_px.push_back(gpt.px());
          GenPart_py.push_back(gpt.py());
          GenPart_pz.push_back(gpt.pz());
	      numGenPart++;
         }
        }
    } 

  mtree->Fill();
  return;      
          
 
}


// ------------ method called once each job just before starting event loop  ------------
void
GenParticleAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
GenParticleAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenParticleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(GenParticleAnalyzer);
