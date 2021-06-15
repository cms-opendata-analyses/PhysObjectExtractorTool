// -*- C++ -*-
//
// Package:    VertexAnalyzer
// Class:      VertexAnalyzer
// 
/**\class VertexAnalyzer VertexAnalyzer.cc Vertex/VertexAnalyzer/src/VertexAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Sat Jun 12 11:03:58 CEST 2021
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <vector>
//
// class declaration
//

class VertexAnalyzer : public edm::EDAnalyzer {
   public:
      explicit VertexAnalyzer(const edm::ParameterSet&);
      ~VertexAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
     
    std::vector<float> PV_chi2;
    std::vector<float> PV_ndof;
    int PV_npvs;
    int PV_npvsGood;
    std::vector<float> PV_score;
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
VertexAnalyzer::VertexAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

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
   
   PV_chi2.clear();
   PV_ndof.clear();
   PV_score.clear();
   PV_x.clear();
   PV_y.clear();
   PV_z.clear();

   edm::Handle<reco::VertexCollection> Primvertex;
   iEvent.getByLabel("offlinePrimaryVertices",Primvertex);
   
   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel("generalTracks", tracks);
   
   PV_npvsGood=0;
  
   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", beamSpotHandle); 
   reco::BeamSpot vertexBeamSpot= *beamSpotHandle;
   Bsp_z = vertexBeamSpot.z0();

   for (reco::VertexCollection::const_iterator vite = Primvertex->begin(); vite != Primvertex->end(); ++vite)
   { 
     float score=0;
     PV_chi2.push_back(vite->chi2());
     PV_ndof.push_back(vite->ndof());
     PV_x.push_back(vite->x());
     PV_y.push_back(vite->y());
     PV_z.push_back(vite->z());
     for (reco::Vertex::trackRef_iterator iTrack = vite->tracks_begin(); iTrack != vite->tracks_end(); ++iTrack)                                                                                                        {                                                                                                                                                                                                                      const reco::TrackRef trackRef = iTrack->castTo<reco::TrackRef>();  float trackpt = trackRef->pt();                                                                                                                 score += trackpt*trackpt; }                          
  PV_score.push_back(score);
 
 if (!vite->isFake() && vite->isValid() && vite->ndof()>4 && fabs(vite->z()-Bsp_z)<24. && vite->position().Rho() < 2.) 
        ++PV_npvsGood;

  }
  for(unsigned int i=0; i < PV_score.size(); i++){
    std::cout <<"PV_score # "<<i<<"="<<PV_score.at(i)<<std::endl; }
  PV_npvs=Primvertex->size();
  std::cout<<PV_npvs<<std::endl;
  std::cout<<PV_npvsGood<<std::endl;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
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

// ------------ method called when starting to processes a run  ------------
void 
VertexAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
VertexAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
VertexAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
VertexAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
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
