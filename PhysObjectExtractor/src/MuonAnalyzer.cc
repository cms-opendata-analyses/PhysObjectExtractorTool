// -*- C++ -*-
//
// Package:    MuonAnalyzer
// Class:      MuonAnalyzer
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

//classes to extract Muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class MuonAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MuonAnalyzer(const edm::ParameterSet&);
      ~MuonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


      //declare the input tag for MuonCollection
      edm::InputTag muonInput;

      // ----------member data ---------------------------

      TTree *mtree;
      int nummuon; //number of muons in the event
      std::vector<float> muon_e;
      std::vector<float> muon_pt;
      std::vector<float> muon_px;
      std::vector<float> muon_py;
      std::vector<float> muon_pz;
      std::vector<float> muon_eta;
      std::vector<float> muon_phi;
      std::vector<float> muon_ch;
      std::vector<float> muon_mass;
      std::vector<float> muon_pfreliso03all;
      std::vector<float> muon_pfreliso04all;
      std::vector<float> muon_tightid;
      std::vector<float> muon_softid;
      std::vector<float> muon_dxy;
      std::vector<float> muon_dxyErr;
      std::vector<float> muon_dz;
      std::vector<float> muon_dzErr;
      std::vector<float> muon_genpartidx;
      std::vector<float> muon_jetidx;
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

MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& iConfig)

{
//now do what ever initialization is needed
	muonInput = iConfig.getParameter<edm::InputTag>("InputCollection");
	edm::Service<TFileService> fs;
	mtree = fs->make<TTree>("Events", "Events");

	mtree->Branch("numbermuon",&nummuon);
	mtree->GetBranch("numbermuon")->SetTitle("number of muons");
	mtree->Branch("muon_e",&muon_e);
	mtree->GetBranch("muon_e")->SetTitle("muon energy");
	mtree->Branch("muon_pt",&muon_pt);
	mtree->GetBranch("muon_pt")->SetTitle("muon transverse momentum");
 	mtree->Branch("muon_px",&muon_px);
 	mtree->GetBranch("muon_px")->SetTitle("muon momentum x-component");
 	mtree->Branch("muon_py",&muon_py);
 	mtree->GetBranch("muon_py")->SetTitle("muon momentum y-component");
	mtree->Branch("muon_pz",&muon_pz);
	mtree->GetBranch("muon_pz")->SetTitle("muon momentum z-component");
 	mtree->Branch("muon_eta",&muon_eta);
 	mtree->GetBranch("muon_eta")->SetTitle("muon pseudorapidity");
	mtree->Branch("muon_phi",&muon_phi);
	mtree->GetBranch("muon_phi")->SetTitle("muon polar angle");
	mtree->Branch("muon_ch",&muon_ch);
	mtree->GetBranch("muon_ch")->SetTitle("muon charge");
	mtree->Branch("muon_mass",&muon_mass);
	mtree->GetBranch("muon_mass")->SetTitle("muon mass");
	mtree->Branch("muon_pfreliso03all",&muon_pfreliso03all);
	mtree->GetBranch("muon_pfreliso03all")->SetTitle("muon pfreliso03all");
	mtree->Branch("muon_pfreliso04all",&muon_pfreliso04all);
	mtree->GetBranch("muon_pfreliso04all")->SetTitle("muon pfreliso04all");
	mtree->Branch("muon_tightid",&muon_tightid);
	mtree->GetBranch("muon_tightid")->SetTitle("muon tight id");
	mtree->Branch("muon_softid",&muon_softid);
	mtree->GetBranch("muon_softid")->SetTitle("muon soft id");
	mtree->Branch("muon_dxy",&muon_dxy);
	mtree->GetBranch("muon_dxy")->SetTitle("muon dxy");
	mtree->Branch("muon_dxyErr",&muon_dxyErr);
	mtree->GetBranch("muon_dxyErr")->SetTitle("muon dxyErr");
	mtree->Branch("muon_dz",&muon_dz);
	mtree->GetBranch("muon_dz")->SetTitle("muon dz");
	mtree->Branch("muon_dzErr",&muon_dzErr);
	mtree->GetBranch("muon_dzErr")->SetTitle("muon dzErr");
	mtree->Branch("muon_jetidx",&muon_jetidx);
	mtree->GetBranch("muon_jetidx")->SetTitle("muon jetidx");
	mtree->Branch("muon_genpartidx",&muon_genpartidx);
	mtree->GetBranch("muon_genpartidx")->SetTitle("muon genpartidx");
}

MuonAnalyzer::~MuonAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<reco::MuonCollection> mymuons;
   iEvent.getByLabel(muonInput, mymuons);

   nummuon = 0;
   muon_e.clear();
   muon_pt.clear();
   muon_px.clear();
   muon_py.clear();
   muon_pz.clear();
   muon_eta.clear();
   muon_phi.clear();
   muon_ch.clear();
   muon_mass.clear();
   muon_pfreliso03all.clear();
   muon_pfreliso04all.clear();
   muon_tightid.clear();
   muon_softid.clear();
   muon_dxy.clear();
   muon_dxyErr.clear();
   muon_dz.clear();
   muon_dzErr.clear();
   muon_jetidx.clear();
   muon_genpartidx.clear();
   
   Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);


  if(mymuons.isValid()){
  
     // get the number of muons in the event
     nummuon=(*mymuons).size();
     const float mu_min_pt = 3;
     math::XYZPoint pv(vertices->begin()->position());
     
        for (reco::MuonCollection::const_iterator itmuon=mymuons->begin(); itmuon!=mymuons->end(); ++itmuon){
           if (itmuon->pt() > mu_min_pt) {

        	    muon_e.push_back(itmuon->energy());
        	    muon_pt.push_back(itmuon->pt());
        	    muon_px.push_back(itmuon->px());
        	    muon_py.push_back(itmuon->py());
        	    muon_pz.push_back(itmuon->pz());
        	    muon_eta.push_back(itmuon->eta());
        	    muon_phi.push_back(itmuon->phi());
        	    muon_ch.push_back(itmuon->charge());
        	    muon_mass.push_back(itmuon->mass());
        	    
        	    if (itmuon->isPFMuon() && itmuon->isPFIsolationValid()) {
        	      auto iso03 = itmuon->pfIsolationR03();
        	      muon_pfreliso03all.push_back((iso03.sumChargedHadronPt + iso03.sumNeutralHadronEt + iso03.sumPhotonEt)/itmuon->pt());
        	      auto iso04 = itmuon->pfIsolationR04();
        	      muon_pfreliso04all.push_back((iso04.sumChargedHadronPt + iso04.sumNeutralHadronEt + iso04.sumPhotonEt)/itmuon->pt());
        	    } else {
        	      muon_pfreliso03all.push_back(-999);
        	      muon_pfreliso04all.push_back(-999);
        	    }

                    muon_tightid.push_back(muon::isTightMuon(*itmuon, *vertices->begin()));
        	    muon_softid.push_back(muon::isSoftMuon(*itmuon, *vertices->begin()));
        	    auto trk = itmuon->globalTrack();
        	    if (trk.isNonnull()) {
        	      muon_dxy.push_back(trk->dxy(pv));
        	      muon_dz.push_back(trk->dz(pv));
        	      muon_dxyErr.push_back(trk->d0Error());
        	      muon_dzErr.push_back(trk->dzError());
        	    } else {
        	      muon_dxy.push_back(-999);
        	      muon_dxyErr.push_back(-999);
        	      muon_dz.push_back(-999);
        	      muon_dzErr.push_back(-999);
        	    }
        	    muon_genpartidx.push_back(-1);
        	    muon_jetidx.push_back(-1);
          } 
        }
  }
	
  mtree->Fill();
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void
MuonAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
MuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
MuonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
MuonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
MuonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonAnalyzer);
