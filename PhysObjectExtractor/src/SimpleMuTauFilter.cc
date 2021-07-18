// -*- C++ -*-
//
// Package:    SimpleMuTauFilter
// Class:      SimpleMuTauFilter
// 
/**\class SimpleMuTauFilter SimpleMuTauFilter.cc PhysObjectExtractorTool/SimpleMuTauFilter/src/SimpleMuTauFilter.cc

 Description: [one line class summary]

This is a simple filter example to filter on at least a muon and at least
a tau of certain characteristics, which are mostly configurable.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Sat Jul 17 22:23:23 CEST 2021
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//classes to extract Muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//classes to extract tau information
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

//
// class declaration
//

class SimpleMuTauFilter : public edm::EDFilter {
   public:
      explicit SimpleMuTauFilter(const edm::ParameterSet&);
      ~SimpleMuTauFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  edm::InputTag muonInput;
  edm::InputTag tauInput;
  double mu_minpt_;
  double mu_etacut_;
  double tau_minpt_;
  double tau_etacut_;
  

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
SimpleMuTauFilter::SimpleMuTauFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  muonInput = iConfig.getParameter<edm::InputTag>("InputCollectionMuons");
  tauInput = iConfig.getParameter<edm::InputTag>("InputCollectionTaus");
  mu_minpt_ = iConfig.getParameter<double>("mu_minpt");
  mu_etacut_ = iConfig.getParameter<double>("mu_etacut");
  tau_minpt_ = iConfig.getParameter<double>("tau_minpt");
  tau_etacut_ = iConfig.getParameter<double>("tau_etacut");

}


SimpleMuTauFilter::~SimpleMuTauFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SimpleMuTauFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;


 //Filter on at least one good muon
  Handle<reco::MuonCollection> mymuons;
  iEvent.getByLabel(muonInput, mymuons);

  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);

  bool isGoodMuon =false;

  if(mymuons.isValid()){
    math::XYZPoint pv(vertices->begin()->position());
    for (reco::MuonCollection::const_iterator itmuon=mymuons->begin(); itmuon!=mymuons->end(); ++itmuon){
      if(abs(itmuon->eta())<mu_etacut_ &&
	 itmuon->pt()>mu_minpt_ &&
	 bool(muon::isTightMuon(*itmuon, *vertices->begin()))){
	isGoodMuon = true;
      }
    }
  }

   //Filter on at least one good tau
  Handle<reco::PFTauCollection> mytaus;
  iEvent.getByLabel(tauInput, mytaus);

  Handle<reco::PFTauDiscriminator> tausLooseIso, tausVLooseIso, tausMediumIso, tausTightIso, tausTightEleRej, tausTightMuonRej, tausDecayMode, tausRawIso;
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByDecayModeFinding"),tausDecayMode);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr"), tausTightIso);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByTightElectronRejection"), tausTightEleRej);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByTightMuonRejection"), tausTightMuonRej); 

  bool isGoodTau = false;
  
  if(mytaus.isValid()){
    for (reco::PFTauCollection::const_iterator itTau=mytaus->begin(); itTau!=mytaus->end(); ++itTau){
      const auto idx = itTau - mytaus->begin();
      if(itTau->charge()!=0 &&
	 abs(itTau->eta())<tau_etacut_ &&
	 itTau->pt()>tau_minpt_ &&
	 bool(tausDecayMode->operator[](idx).second) &&
	 bool(tausTightIso->operator[](idx).second) &&
	 bool(tausTightEleRej->operator[](idx).second) &&
	 bool(tausTightMuonRej->operator[](idx).second) ){
	isGoodTau = true;
      }
    }
  }

 
  


  using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   return isGoodMuon*isGoodTau;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SimpleMuTauFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimpleMuTauFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
SimpleMuTauFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
SimpleMuTauFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
SimpleMuTauFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
SimpleMuTauFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleMuTauFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SimpleMuTauFilter);
