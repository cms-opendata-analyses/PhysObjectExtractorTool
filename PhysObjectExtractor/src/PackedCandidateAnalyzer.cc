// -*- C++ -*-
//
// Package:    PackedCandidate/PackedCandidateAnalyzer
// Class:      PackedCandidateAnalyzer
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

//class to extract Packed Candidate information
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/iterator.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include <vector>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class PackedCandidateAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PackedCandidateAnalyzer(const edm::ParameterSet&);
      ~PackedCandidateAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<pat::PackedCandidateCollection> packedToken_;

      // ----------member data ---------------------------
      
     TTree *mtree;
      
    int numCandidates; //number of particles
    std::vector<float> packed_pt; // transverse momentum
    std::vector<float> packed_eta; 
    std::vector<float> packed_mass; // mass of the particle
    std::vector<float> packed_energy;
    std::vector<float> packed_phi; // momentum azimuthal angle
    std::vector<float> packed_ch; // charge of particle
    std::vector<float> packed_px; // x component of momentum
    std::vector<float> packed_py; // y component of momentum
    std::vector<float> packed_pz; // z component of momentum
    std::vector<float> packed_theta; // polar angle momentum
    std::vector<float> packed_vx; // x coordinate of vertex position
    std::vector<float> packed_vy; // y coordinate of vertex position
    std::vector<float> packed_vz; // z coordinate of vertex position
    std::vector<float> packed_lostInnerHits;
    std::vector<float> packed_PuppiWeight;
    std::vector<float> packed_PuppiWeightNoLep;
    std::vector<float> packed_hcalFraction;
    std::vector<int> packed_pdgId;
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
PackedCandidateAnalyzer::PackedCandidateAnalyzer(const edm::ParameterSet& iConfig):
   packedToken_(consumes <pat::PackedCandidateCollection> (iConfig.getParameter<edm::InputTag>("packed")))
{
   //now do what ever initialization is needed
   
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
  
   mtree->Branch("numCandidates",&numCandidates);
   mtree->GetBranch("numCandidates")->SetTitle("number of packed particles");
   mtree->Branch("packed_pt",&packed_pt);
   mtree->GetBranch("packed_pt")->SetTitle("the flowing particle's transverse momentum");
   mtree->Branch("packed_eta",&packed_eta);
   mtree->GetBranch("packed_eta")->SetTitle("generator particle pseudorapidity");
   mtree->Branch("packed_mass",&packed_mass);
   mtree->GetBranch("packed_mass")->SetTitle("packed candidate mass");
   mtree->Branch("packed_energy",&packed_energy);
   mtree->GetBranch("packed_energy")->SetTitle("packed candidate energy");
   mtree->Branch("packed_ch",&packed_ch);
   mtree->GetBranch("packed_ch")->SetTitle("charge of the particle");
   mtree->Branch("packed_phi",&packed_phi);
   mtree->GetBranch("packed_phi")->SetTitle("packed candidate azimuthal angle of momentum vector");
   mtree->Branch("packed_px",&packed_px);
   mtree->GetBranch("packed_px")->SetTitle("packed candidate x coordinate of momentum vector");
   mtree->Branch("packed_py",&packed_py);
   mtree->GetBranch("packed_py")->SetTitle("packed candidate y coordinate of momentum vector");
   mtree->Branch("packed_pz",&packed_pz);
   mtree->GetBranch("packed_pz")->SetTitle("packed candidate z coordinate of momentum vector");
   mtree->Branch("packed_theta", &packed_theta);
   mtree->GetBranch("packed_theta")->SetTitle("polar angle of momentum");
   mtree->Branch("packed_vx",&packed_vx);
   mtree->GetBranch("packed_vx")->SetTitle("x of the vertex");
   mtree->Branch("packed_vy",&packed_vy);
   mtree->GetBranch("packed_vy")->SetTitle("y of the vertex");
   mtree->Branch("packed_vz",&packed_vz);
   mtree->GetBranch("packed_vz")->SetTitle("z of the vertex");
   mtree->Branch("packed_lostInnerHits",&packed_lostInnerHits);
   mtree->GetBranch("packed_lostInnerHits")->SetTitle("Inner hit information");
   mtree->Branch("packed_PuppiWeight",&packed_PuppiWeight);
   mtree->GetBranch("packed_PuppiWeight")->SetTitle("Puppi Weight");
   mtree->Branch("packed_PuppiWeightNoLep",&packed_PuppiWeightNoLep);
   mtree->GetBranch("packed_PuppiWeightNoLep")->SetTitle("Puppi MET");
   mtree->Branch("packed_hcalFraction",&packed_hcalFraction);
   mtree->GetBranch("packed_hcalFraction")->SetTitle("Energy Measured from Hadronic Calorimeter");
   mtree->Branch("packed_pdgId",&packed_pdgId);
   mtree->GetBranch("packed_pdgId")->SetTitle("PDG ID of particle-- 11,13,22 - ele/mu/gamma, 211-charged hadrons, 130-neutral hadrons, 1&2-hadronic and em particles in HF");
}


PackedCandidateAnalyzer::~PackedCandidateAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PackedCandidateAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<pat::PackedCandidateCollection> packed;
   iEvent.getByToken(packedToken_,packed);

   numCandidates=0;
   packed_pt.clear();
   packed_eta.clear();
   packed_mass.clear();
   packed_energy.clear();
   packed_phi.clear();
   packed_ch.clear();
   packed_px.clear();
   packed_py.clear();
   packed_pz.clear();
   packed_theta.clear();
   packed_vx.clear();
   packed_vy.clear();
   packed_vz.clear();
   packed_lostInnerHits.clear();
   packed_PuppiWeight.clear();
   packed_PuppiWeightNoLep.clear();
   packed_hcalFraction.clear();
   packed_pdgId.clear();

   
   for (const pat::PackedCandidate &pack : *packed)
   {
         packed_pt.push_back(pack.pt());
         packed_eta.push_back(pack.eta());
         packed_mass.push_back(pack.mass());
         packed_energy.push_back(pack.energy());
         packed_phi.push_back(pack.phi());
         packed_ch.push_back(pack.charge());
         packed_px.push_back(pack.px());
         packed_py.push_back(pack.py());
         packed_pz.push_back(pack.pz());
         packed_theta.push_back(pack.theta());
         packed_vx.push_back(pack.vx());
         packed_vy.push_back(pack.vy());
         packed_vz.push_back(pack.vz());
         packed_lostInnerHits.push_back(pack.lostInnerHits());
         packed_PuppiWeight.push_back(pack.puppiWeight());
         packed_PuppiWeightNoLep.push_back(pack.puppiWeightNoLep());
         packed_hcalFraction.push_back(pack.hcalFraction());
         packed_pdgId.push_back(pack.pdgId());
	      numCandidates++;
    } 

  mtree->Fill();
  return;      
          
 
}


// ------------ method called once each job just before starting event loop  ------------
void
PackedCandidateAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
PackedCandidateAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PackedCandidateAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(PackedCandidateAnalyzer);