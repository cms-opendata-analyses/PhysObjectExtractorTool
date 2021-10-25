// -*- C++ -*-
//
// Package:    Vertex/VertexAnalyzer
// Class:      VertexAnalyzer
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

//class to extract Vertex information
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"


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



class VertexAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit VertexAnalyzer(const edm::ParameterSet&);
      ~VertexAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamToken_;

      // ----------member data ---------------------------
      
    TTree *mtree; 
    std::vector<float> PV_chi2;
    std::vector<float> PV_ndof;
    int PV_npvs;
    int PV_npvsGood;
    std::vector<float> PV_x;
    std::vector<float> PV_y; 
    std::vector<float> PV_z;
    Float_t Bsp_z;
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
VertexAnalyzer::VertexAnalyzer(const edm::ParameterSet& iConfig): 
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   beamToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beams")))
   
   
{
   //now do what ever initialization is needed
   
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
  
   mtree->Branch("PV_chi2",&PV_chi2);
   mtree->GetBranch("PV_chi2")->SetTitle("main primary vertex chi2");
   mtree->Branch("PV_ndof",&PV_ndof);
   mtree->GetBranch("PV_ndof")->SetTitle("main primary vertex number of degree of freedom");
   mtree->Branch("PV_npvs",&PV_npvs);
   mtree->GetBranch("PV_npvs")->SetTitle("total number of reconstructed primary vertices");
   mtree->Branch("PV_npvsGood",&PV_npvsGood);
   mtree->GetBranch("PV_npvsGood")->SetTitle("number of good reconstructed primary vertices. selection:!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2");
   mtree->Branch("PV_x",&PV_x);
   mtree->GetBranch("PV_x")->SetTitle("main primary vertex x coordinate");
   mtree->Branch("PV_y",&PV_y);
   mtree->GetBranch("PV_y")->SetTitle("main primary vertex y coordinate");
   mtree->Branch("PV_z",&PV_z);
   mtree->GetBranch("PV_z")->SetTitle("main primary vertex z coordinate");
}


VertexAnalyzer::~VertexAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
VertexAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return; // skip the event if no PV found
    

   Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByToken(beamToken_, beamSpotHandle);
   reco::BeamSpot vertexBeamSpot= *beamSpotHandle;

   PV_chi2.clear();
   PV_ndof.clear();
   PV_x.clear();
   PV_y.clear();
   PV_z.clear();
   PV_npvs=0;
   PV_npvsGood=0;
   
   Bsp_z = vertexBeamSpot.z0();
     for (const reco::Vertex &vtx : *vertices)
    {
       PV_chi2.push_back(vtx.chi2());
       PV_ndof.push_back(vtx.ndof());
       PV_x.push_back(vtx.x());
       PV_y.push_back(vtx.y());
       PV_z.push_back(vtx.z());
       ++PV_npvs;
       if (!vtx.isFake() && vtx.isValid() && vtx.ndof()>4 && fabs(vtx.z()-Bsp_z)<24. && vtx.position().Rho() < 2.)
       ++PV_npvsGood;
    } 

  mtree->Fill();
  return;      
          
 
}


// ------------ method called once each job just before starting event loop  ------------
void
VertexAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
VertexAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VertexAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexAnalyzer);
