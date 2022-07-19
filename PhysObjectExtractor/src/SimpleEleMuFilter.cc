// -*- C++ -*-
//
// Package:    PhysObjectExtractor/SimpleEleMuFilter
// Class:      SimpleEleMuFilter
// 
/**\class SimpleEleMuFilter SimpleEleMuFilter.cc plugins/SimpleEleMuFilter.cc

 Description: [This is a simple filter example to filter in one energetic electron or muon]]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Mon, 18 Jul 2022 09:02:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


//
// class declaration
//

class SimpleEleMuFilter : public edm::stream::EDFilter<> {
   public:
      explicit SimpleEleMuFilter(const edm::ParameterSet&);
      ~SimpleEleMuFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      
      

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      double mu_minpt_;
      double mu_etacut_; 
      double ele_minpt_;
      double ele_etacut_;

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
SimpleEleMuFilter::SimpleEleMuFilter(const edm::ParameterSet& iConfig):
electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
{
   //now do what ever initialization is needed
   mu_minpt_ = iConfig.getParameter<double>("mu_minpt");
  mu_etacut_ = iConfig.getParameter<double>("mu_etacut");
  ele_minpt_ = iConfig.getParameter<double>("ele_minpt");
  ele_etacut_ = iConfig.getParameter<double>("ele_etacut");

}


SimpleEleMuFilter::~SimpleEleMuFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SimpleEleMuFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   //at least one vertex
   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return false; // skip the event if no PV is found
   const reco::Vertex &PV = vertices->front();

   //check if there is at least a muon passing the conditions
   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   bool aGoodMuon = false;
   for (const pat::Muon &mu : *muons){
     if ((abs(mu.eta()))<mu_etacut_ && 
	 mu.pt()>mu_minpt_ &&
	 mu.isTightMuon(PV)){
       aGoodMuon = true;
       break;   
     }
   }

   //check if here is at least an electron passing the conditions
   //unless we already found a good muon
   bool aGoodElectron = false;
   if (!aGoodMuon){
     Handle<pat::ElectronCollection> electrons;
     iEvent.getByToken(electronToken_, electrons);
     for (const pat::Electron &el : *electrons){
       if ((abs(el.eta()))<ele_etacut_ &&
	   el.pt()>ele_minpt_ &&
	   el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight")){
	 aGoodElectron = true;
	 break;   
       }
     }
   }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   return (aGoodMuon||aGoodElectron);
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
SimpleEleMuFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
SimpleEleMuFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
SimpleEleMuFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
SimpleEleMuFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SimpleEleMuFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SimpleEleMuFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleEleMuFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SimpleEleMuFilter);
