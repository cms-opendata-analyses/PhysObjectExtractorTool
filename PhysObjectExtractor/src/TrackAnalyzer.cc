// -*- C++ -*-
//
// Package:    TrackAnalyzer
// Class:      TrackAnalyzer
// 
/**\class TrackAnalyzer TrackAnalyzer.cc Track/TrackAnalyzer/src/TrackAnalyzer.cc

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

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"


#include "TTree.h"
#include "TFile.h"
#include <vector>
//
// class declaration
//

class TrackAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TrackAnalyzer(const edm::ParameterSet&);
      ~TrackAnalyzer();

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
    
    TTree *mtree;
    int numtracks; 
    std::vector<float> track_pt;
    std::vector<float> track_ptError;
    std::vector<float> track_charge;
    std::vector<float> track_chi2;
    std::vector<float> track_eta;
    std::vector<float> track_etaError;
    std::vector<float> track_lambda;
    std::vector<float> track_lambdaError;
    std::vector<float> track_ndof;
    std::vector<float> track_phi;
    std::vector<float> track_phiError;
    std::vector<float> track_px;
    std::vector<float> track_py;
    std::vector<float> track_pz;
    std::vector<float> track_theta;
    std::vector<float> track_thetaError;
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
TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
   mtree->Branch("numtracks",&numtracks);
   mtree->GetBranch("numtracks")->SetTitle("number of tracks");
   mtree->Branch("track_pt",&track_pt);
   mtree->GetBranch("track_pt")->SetTitle("track transverse momentum");
   mtree->Branch("track_ptError",&track_ptError);
   mtree->GetBranch("track_ptError")->SetTitle("error on track transverse momentum");
   mtree->Branch("track_charge",&track_charge);
   mtree->GetBranch("track_charge")->SetTitle("track charge"); 
   mtree->Branch("track_chi2",&track_chi2);
   mtree->GetBranch("track_chi2")->SetTitle("track chi2");
   mtree->Branch("track_eta",&track_eta);
   mtree->GetBranch("track_eta")->SetTitle("track pseudorapidity of momentum vector");  
   mtree->Branch("track_etaError",&track_etaError);
   mtree->GetBranch("track_etaError")->SetTitle("track error on pseudorapidity of momentum vector"); 
   mtree->Branch("track_lambda",&track_lambda);
   mtree->GetBranch("track_lambda")->SetTitle("track lambda");
   mtree->Branch("track_lambdaError",&track_lambdaError);
   mtree->GetBranch("track_lambdaError")->SetTitle("error on track lambda");
   mtree->Branch("track_ndof",&track_ndof);
   mtree->GetBranch("track_ndof")->SetTitle("track number of degrees of freedom");
   mtree->Branch("track_phi",&track_phi);
   mtree->GetBranch("track_phi")->SetTitle("track azimuthal angle of momentum vector");
   mtree->Branch("track_phiError",&track_phiError);
   mtree->GetBranch("track_phiError")->SetTitle("error on track azimuthal angle of momentum vector"); 
   mtree->Branch("track_px",&track_px);
   mtree->GetBranch("track_px")->SetTitle("track x coordinate of momentum vector");
   mtree->Branch("track_py",&track_py);
   mtree->GetBranch("track_py")->SetTitle("track y coordinate of momentum vector"); 
   mtree->Branch("track_pz",&track_pz);
   mtree->GetBranch("track_pz")->SetTitle("track z coordinate of momentum vector");
   mtree->Branch("track_theta",&track_theta);
   mtree->GetBranch("track_theta")->SetTitle("track polar angle");
   mtree->Branch("track_thetaError",&track_thetaError);
   mtree->GetBranch("track_thetaError")->SetTitle("error on track polar angle");
}


TrackAnalyzer::~TrackAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   numtracks=0;
   track_pt.clear();
   track_charge.clear();
   track_ptError.clear();
   track_chi2.clear();
   track_eta.clear();
   track_etaError.clear();
   track_lambda.clear();
   track_lambdaError.clear();
   track_ndof.clear();
   track_phi.clear();
   track_phiError.clear();
   track_px.clear();
   track_py.clear();
   track_pz.clear();
   track_theta.clear();
   track_thetaError.clear();

   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel("generalTracks", tracks);

   if(tracks.isValid())
   {
      numtracks=tracks->size();
   for (reco::TrackCollection::const_iterator iTrack = tracks->begin(); iTrack != tracks->end(); ++iTrack)
      {
        track_pt.push_back(iTrack->pt());
        track_ptError.push_back(iTrack->ptError());
        track_charge.push_back(iTrack->charge());
        track_chi2.push_back(iTrack->chi2());
        track_eta.push_back(iTrack->eta());
        track_etaError.push_back(iTrack->etaError());
        track_lambda.push_back(iTrack->lambda());
        track_lambdaError.push_back(iTrack->lambdaError());
        track_ndof.push_back(iTrack->ndof());
        track_phi.push_back(iTrack->phi());
        track_phiError.push_back(iTrack->phiError());
        track_px.push_back(iTrack->px());
        track_py.push_back(iTrack->py());
        track_pz.push_back(iTrack->pz());
        track_theta.push_back(iTrack->theta());
        track_thetaError.push_back(iTrack->thetaError());
      }
   }


 mtree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
TrackAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TrackAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void
TrackAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
TrackAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
TrackAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TrackAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
