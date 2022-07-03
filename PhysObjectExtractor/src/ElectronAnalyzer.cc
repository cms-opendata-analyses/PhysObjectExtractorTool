// -*- C++ -*-
//
// Package:    Electron/ElectronAnalyzer
// Class:      ElectronAnalyzer
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

#include "DataFormats/Common/interface/ValueMap.h"

//class to extract electron information
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

//Transient track for impact parameter
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

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


//using reco::TrackCollection;

class ElectronAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ElectronAnalyzer(const edm::ParameterSet&);
      ~ElectronAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void myVar(//const pat::Electron&,
                     const reco::Vertex&,
                     const TransientTrackBuilder&, const edm::Event& /*, const edm::EventSetup&*/);
      virtual void analyze(const edm::Event&, const edm::EventSetup& ) override;
      virtual void endJob() override;

      edm::EDGetTokenT<pat::ElectronCollection> electronToken_, electronToken2_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;

      // ID decisions objects
     // edm::EDGetTokenT<edm::ValueMap<bool>> eleLooseIdMapToken_;
     // edm::EDGetTokenT<edm::ValueMap<bool>> eleTightIdMapToken_;;

      // MVA values and categories (optional)
     // edm::EDGetTokenT<edm::ValueMap<float>> mvaValuesMapToken_;
     // edm::EDGetTokenT<edm::ValueMap<int>> mvaCategoriesMapToken_;
      // ----------member data ---------------------------

      // embed various impact parameters with errors
    // embed high level selection

      TTree *mtree;
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
      std::vector<float> electron_ecalIso;
     // std::vector<float> electron_mvaValue;
     // std::vector<int> electron_mvaCategory;
      std::vector<int> electron_ismvaLoose;
      std::vector<int> electron_ismvaTight;
      std::vector<float> electron_ip3d;

      Float_t myMVAVar_ip3d;
      Float_t myMVAVar_ip3dSig;
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
ElectronAnalyzer::ElectronAnalyzer(const edm::ParameterSet& iConfig):
 electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
 vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
// eleLooseIdMapToken_(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
// eleTightIdMapToken_(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleTightIdMap")))
// mvaValuesMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
// mvaCategoriesMapToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap")))
{
   //now do what ever initialization is needed

   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");


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
   mtree->Branch("electron_ecalIso",&electron_ecalIso);
   mtree->GetBranch("electron_ecalIso")->SetTitle("electron Ecal Reconstruction Hit");
   //mtree->Branch("electron_mvaValue",&electron_mvaValue);
   //mtree->GetBranch("electron_mvaValue")->SetTitle("electron mva Value map");
   //mtree->Branch("electron_mvaCategory",&electron_mvaCategory);
   //mtree->GetBranch("electron_mvaCategory")->SetTitle("electron mva category map");
   mtree->Branch("electron_ismvaLoose",&electron_ismvaLoose);
   mtree->GetBranch("electron_ismvaLoose")->SetTitle("electron mva Loose");
   mtree->Branch("electron_ismvaTight",&electron_ismvaTight);
   mtree->GetBranch("electron_ismvaTight")->SetTitle("electron mva Tight");
   mtree->Branch("electron_ip3d",&electron_ip3d);
   mtree->GetBranch("electron_ip3d")->SetTitle("electron impact parameter");
}

//Destructor
ElectronAnalyzer::~ElectronAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
   using namespace edm;

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   math::XYZPoint pv(vertices->begin()->position());
   // Get the electron ID data from the event stream.
   // Note: this implies that the VID ID modules have been run upstream.
   // If you need more info, check with the EGM group.
   //edm::Handle<edm::ValueMap<bool>> loose_id_decisions;
   //edm::Handle<edm::ValueMap<bool>> tight_id_decisions;
   //iEvent.getByToken(eleLooseIdMapToken_,loose_id_decisions);
   //iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);
   // Get MVA values and categories (optional)
   //edm::Handle<edm::ValueMap<float>> mvaValues;
   //edm::Handle<edm::ValueMap<int>> mvaCategories;
  // iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  // iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);

//   edm::Handle<edm::View<reco::GsfElectron>> elec;
//   iEvent.getByToken(electronsMiniAODToken_,elec);
//   bool isPassLoose = (*loose_id_decisions)[el];
//   bool isPassTight  = (*tight_id_decisions)[el];

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
   electron_ecalIso.clear();
  // electron_mvaValue.clear();
  // electron_mvaCategory.clear();
   electron_ismvaLoose.clear();
   electron_ismvaTight.clear();

    for (const pat::Electron &el : *electrons)
    {
     // bool isPassLoose = (*loose_id_decisions)[el];
     // bool isPassTight  = (*tight_id_decisions)[el];

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
     // electron_mvaValue.push_back(el.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp90"));
     // electron_mvaCategory.push_back(el.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp90"));
      electron_ismvaLoose.push_back(el.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp90"));
      electron_ismvaTight.push_back(el.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp80"));

      numelectron++;
    }

  for ( const reco::GsfElectron &ele : *electrons)
  {
      electron_ecalIso.push_back(ele.dr03EcalRecHitSumEt());
      //electron_mvaValue.push_back((*mvaValues)[ele]);
      //electron_mvaCategory.push_back((*mvaCatagories)[ele]);
  }

  mtree->Fill();
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void ElectronAnalyzer::myVar(//const reco::GsfElectron& ele,
			     const reco::Vertex& vertex,
			     const TransientTrackBuilder& transientTrackBuilder,
                             const edm::Event& iEvent/*, const edm::EventSetup& iSetup*/)
{
   using namespace edm;

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken2_, electrons);

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   math::XYZPoint pv(vertices->begin()->position());

   numelectron = 0;
   electron_ip3d.clear();

   //default values for IP3D
   for(const pat::Electron &el : *electrons)
   {
     myMVAVar_ip3d = -999.;
     // myMVAVar_ip3dSig = 0.0;
     if (el.gsfTrack().isNonnull()) {
     const double gsfsign   = ( (-el.gsfTrack()->dxy(pv))   >=0 ) ? 1. : -1.;

     const reco::TransientTrack& tt = transientTrackBuilder.build(el.gsfTrack());
     const std::pair<bool,Measurement1D>& ip3dpv =  IPTools::absoluteImpactParameter3D(tt,vertex);
     if (ip3dpv.first) {
       //double ip3d = gsfsign*ip3dpv.second.value();
	//double ip3derr = ip3dpv.second.error();
       electron_ip3d.push_back(gsfsign*ip3dpv.second.value());
      // electron_ip3d.push_back(myMVAVar_ip3d);
       // myMVAVar_ip3dSig = ip3d/ip3derr;
       }
      }
   numelectron++;
  }
   mtree->Fill();
   return;
}
// ------------ method called once each job just before starting event loop  ------------
void
ElectronAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ElectronAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronAnalyzer);
