// -*- C++ -*-
//
// Package:    Electron/ElectronAnalyzer
// Class:      ElectronAnalyzer
//
 
// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//class to extract electron information
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include "Math/VectorUtil.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


//using reco::TrackCollection;

class ElectronNtuplerVIDwithMVADemo : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ElectronNtuplerVIDwithMVADemo(const edm::ParameterSet&);
      ~ElectronNtuplerVIDwithMVADemo();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   enum ElectronMatchType {UNMATCHED =0,//
			  TRUE_PROMPT_ELECTRON, //
			  TRUE_ELECTRON_FROM_TAU,//
			  TRUE_NON_PROMPT_ELECTRON};//

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      int matchToTruth(const edm::Ptr<reco::GsfElectron> el, //
		       const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);//

      void findFirstNonElectronMother(const reco::Candidate *particle,//
				    int &ancestorPID, int &ancestorStatus);//


      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetToken electronsMiniAODToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;//

      // MVA values and categories (optional)
      edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;//
      edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;//

      // ----------member data ---------------------------

      TTree *mtree;

      //Global info
      Int_t run_;//
      Int_t lumi_;//
      Int_t evtnum_;//

      // all the variables for the output tree
      int numelectron; //number of electrons in the event

      std::vector<float> electron_e;
      std::vector<float> electron_pt;
      std::vector<float> electron_px;
      std::vector<float> electron_py;
      std::vector<float> electron_pz;
      std::vector<float> electron_eta;
      std::vector<float> electron_phi;
      std::vector<float> electron_ch;
      std::vector<float> electron_iso;
      std::vector<bool> electron_veto;//
      std::vector<bool> electron_isLoose;
      std::vector<bool> electron_isMedium;
      std::vector<bool> electron_isTight;
      std::vector<float> electron_dxy;
      std::vector<float> electron_dz;
      std::vector<float> electron_dxyError;
      std::vector<float> electron_dzError;

      std::vector<float> mvaValue_;//
      std::vector<int> mvaCategory_;//
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
ElectronNtuplerVIDwithMVADemo::ElectronNtuplerVIDwithMVADemo(const edm::ParameterSet& iConfig):
 electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
 vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
 mvaValuesMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
 mvaCategoriesMapToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap")))
{

   //miniAOD Tokens
   electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electronsMiniAOD"));

   genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));

   //now do what ever initialization is needed

   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");

  mtree->Branch("run"        ,  &run_     , "run/I");//
  mtree->Branch("lumi"       ,  &lumi_    , "lumi/I");//
  mtree->Branch("evtnum"     ,  &evtnum_  , "evtnum/I");//
 
  mtree->Branch("numberelectron",&numelectron);
  mtree->GetBranch("numberelectron")->SetTitle("number of electrons");
  mtree->Branch("electron_e",&electron_e);
  mtree->GetBranch("electron_e")->SetTitle("electron energy");
  mtree->Branch("electron_pt",&electron_pt);
  mtree->GetBranch("electron_pt")->SetTitle("electron transverse momentum");
  mtree->Branch("electron_px",&electron_px);
  mtree->GetBranch("electron_px")->SetTitle("electron momentum x-component");
  mtree->Branch("electron_py",&electron_py);
  mtree->GetBranch("electron_py")->SetTitle("electron momentum y-component");
  mtree->Branch("electron_pz",&electron_pz);
  mtree->GetBranch("electron_pz")->SetTitle("electron momentum z-component");
  mtree->Branch("electron_eta",&electron_eta);
  mtree->GetBranch("electron_eta")->SetTitle("electron pseudorapidity");
  mtree->Branch("electron_phi",&electron_phi);
  mtree->GetBranch("electron_phi")->SetTitle("electron polar angle");
  mtree->Branch("electron_ch",&electron_ch);
  mtree->GetBranch("electron_ch")->SetTitle("electron charge");
  mtree->Branch("electron_iso",&electron_iso);
  mtree->GetBranch("electron_iso")->SetTitle("electron isolation");
  mtree->Branch("electron_veto",&electron_veto);//
  mtree->GetBranch("electron_veto")->SetTitle("electron veto");//
  mtree->Branch("electron_isLoose",&electron_isLoose);
  mtree->GetBranch("electron_isLoose")->SetTitle("electron tagged loose");
  mtree->Branch("electron_isMedium",&electron_isMedium);
  mtree->GetBranch("electron_isMedium")->SetTitle("electron tagged medium");
  mtree->Branch("electron_isTight",&electron_isTight);
  mtree->GetBranch("electron_isTight")->SetTitle("electron tagged tight");
  mtree->Branch("electron_dxy",&electron_dxy);
  mtree->GetBranch("electron_dxy")->SetTitle("electron transverse plane impact parameter (mm)");
  mtree->Branch("electron_dz",&electron_dz);
  mtree->GetBranch("electron_dz")->SetTitle("electron longitudinal impact parameter (mm)");
  mtree->Branch("electron_dxyError",&electron_dxyError);
  mtree->GetBranch("electron_dxyError")->SetTitle("electron transverse impact parameter uncertainty (mm)");
  mtree->Branch("electron_dzError",&electron_dzError);
  mtree->GetBranch("electron_dzError")->SetTitle("electron longitudinal impact parameter uncertainty (mm)");

  mtree->Branch("mvaVal" ,  &mvaValue_ );//
  mtree->Branch("mvaCat" ,  &mvaCategory_ );//

}


ElectronNtuplerVIDwithMVADemo::~ElectronNtuplerVIDwithMVADemo()
{

   // do anything here that needs to be done at dectruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronNtuplerVIDwithMVADemo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   math::XYZPoint pv(vertices->begin()->position());

   numelectron = 0;
   electron_e.clear();
   electron_pt.clear();
   electron_px.clear();
   electron_py.clear();
   electron_pz.clear();
   electron_eta.clear();
   electron_phi.clear();
   electron_ch.clear();
   electron_iso.clear();
   electron_veto.clear();//
   electron_isLoose.clear();
   electron_isMedium.clear();
   electron_isTight.clear();
   electron_dxy.clear();
   electron_dz.clear();
   electron_dxyError.clear();
   electron_dzError.clear();

    for (const pat::Electron &el : *electrons)
    {
      electron_e.push_back(el.energy());
      electron_pt.push_back(el.pt());
      electron_px.push_back(el.px());
      electron_py.push_back(el.py());
      electron_pz.push_back(el.pz());
      electron_eta.push_back(el.eta());
      electron_phi.push_back(el.phi());
      electron_ch.push_back(el.charge());
      electron_iso.push_back(el.ecalPFClusterIso());
      electron_veto.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-veto"));//
      electron_isLoose.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"));
      electron_isMedium.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"));
      electron_isTight.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight"));
      electron_dxy.push_back(el.gsfTrack()->dxy(pv));
      electron_dz.push_back(el.gsfTrack()->dz(pv));
      electron_dxyError.push_back(el.gsfTrack()->d0Error());
      electron_dzError.push_back(el.gsfTrack()->dzError());
      numelectron++;
    }

  mtree->Fill();
  return;

}


// ------------ method called once each job just before starting event loop  ------------
void
ElectronNtuplerVIDwithMVADemo::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ElectronNtuplerVIDwithMVADemo::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronNtuplerVIDwithMVADemo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronNtuplerVIDwithMVADemo);
