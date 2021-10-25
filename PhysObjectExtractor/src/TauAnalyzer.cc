// -*- C++ -*-
//
// Package:    Tau/TauAnalyzer
// Class:      TauAnalyzer
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

//class to extract tau information
#include "DataFormats/PatCandidates/interface/Tau.h"

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



class TauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TauAnalyzer(const edm::ParameterSet&);
      ~TauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<pat::TauCollection> tauToken_;

      // ----------member data ---------------------------
      
    TTree *mtree;
   int numtau; //number of taus in the event
   std::vector<float> tau_e;
   std::vector<float> tau_pt;
   std::vector<float> tau_px;
   std::vector<float> tau_py;
   std::vector<float> tau_pz;
   std::vector<float> tau_eta;
   std::vector<float> tau_phi;
   std::vector<float> tau_ch;
   std::vector<float> tau_mass;
   std::vector<float> tau_decaymode;
   //std::vector<float> tau_reliso_all;
   std::vector<float> tau_iddecaymode;
   std::vector<float> tau_idisoraw;
   std::vector<float> tau_idisovloose;
   std::vector<float> tau_idisoloose;
   std::vector<float> tau_idisomedium;
   std::vector<float> tau_idisotight;
   std::vector<float> tau_idantieletight;
   std::vector<float> tau_idantimutight;  
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
TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig): 
 tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus")))
   
{
   //now do what ever initialization is needed
   
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
  
   mtree->Branch("numbertau",&numtau);
   mtree->GetBranch("numbertau")->SetTitle("number of taus");
   mtree->Branch("tau_e",&tau_e);
   mtree->GetBranch("tau_e")->SetTitle("tau energy");
   mtree->Branch("tau_pt",&tau_pt);
   mtree->GetBranch("tau_pt")->SetTitle("tau transverse momentum");
   mtree->Branch("tau_px",&tau_px);
   mtree->GetBranch("tau_px")->SetTitle("tau momentum x-component");
   mtree->Branch("tau_py",&tau_py);
   mtree->GetBranch("tau_py")->SetTitle("tau momentum y-component");
   mtree->Branch("tau_pz",&tau_pz);
   mtree->GetBranch("tau_pz")->SetTitle("tau momentum z-component");
   mtree->Branch("tau_eta",&tau_eta);
   mtree->GetBranch("tau_eta")->SetTitle("tau pseudorapidity");
   mtree->Branch("tau_phi",&tau_phi);
   mtree->GetBranch("tau_phi")->SetTitle("tau polar angle");
   mtree->Branch("tau_ch",&tau_ch);
   mtree->GetBranch("tau_ch")->SetTitle("tau charge");
   mtree->Branch("tau_mass",&tau_mass);
   mtree->GetBranch("tau_mass")->SetTitle("tau mass");
   mtree->Branch("tau_decaymode",&tau_decaymode);
   mtree->GetBranch("tau_decaymode")->SetTitle("tau decay mode");
   //mtree->Branch("tau_reliso_all",&tau_reliso_all);
   //mtree->GetBranch("tau_reliso_all")->SetTitle("tau reliso all");
   mtree->Branch("tau_iddecaymode",&tau_iddecaymode);
   mtree->GetBranch("tau_iddecaymode")->SetTitle("tau id decay mode");
   mtree->Branch("tau_idisoraw",&tau_idisoraw);
   mtree->GetBranch("tau_idisoraw")->SetTitle("tau id isoraw");
   mtree->Branch("tau_idisovloose",&tau_idisovloose);
   mtree->GetBranch("tau_idisovloose")->SetTitle("tau id isovloose");
   mtree->Branch("tau_idisoloose",&tau_idisoloose);
   mtree->GetBranch("tau_idisoloose")->SetTitle("tau id isoloose");
   mtree->Branch("tau_idisomedium",&tau_idisomedium);
   mtree->GetBranch("tau_idisomedium")->SetTitle("tau id isomedium");
   mtree->Branch("tau_idisotight",&tau_idisotight);
   mtree->GetBranch("tau_idisotight")->SetTitle("tau id isotight");
   mtree->Branch("tau_idantieletight",&tau_idantieletight);
   mtree->GetBranch("tau_idantieletight")->SetTitle("tau id antieletight");
   mtree->Branch("tau_idantimutight",&tau_idantimutight); 
   mtree->GetBranch("tau_idantimutight")->SetTitle("tau id antimutight");
}


TauAnalyzer::~TauAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<pat::TauCollection> taus;
   iEvent.getByToken(tauToken_, taus);

   numtau = 0;
   tau_e.clear();
   tau_pt.clear();
   tau_px.clear();
   tau_py.clear();
   tau_pz.clear();
   tau_eta.clear();
   tau_phi.clear();
   tau_ch.clear();
   tau_mass.clear();
   tau_decaymode.clear();
   //tau_reliso_all.clear();
   tau_iddecaymode.clear();
   tau_idisoraw.clear();
   tau_idisovloose.clear();
   tau_idisoloose.clear();
   tau_idisomedium.clear();
   tau_idisotight.clear();
   tau_idantieletight.clear();
   tau_idantimutight.clear();

    for (const pat::Tau &tau : *taus)
    {
      tau_e.push_back(tau.energy());
	   tau_pt.push_back(tau.pt());
	   tau_px.push_back(tau.px());
	   tau_py.push_back(tau.py());
	   tau_pz.push_back(tau.pz());
	   tau_eta.push_back(tau.eta());
	   tau_phi.push_back(tau.phi());
	   tau_ch.push_back(tau.charge());
	   tau_decaymode.push_back(tau.decayMode());
	   tau_mass.push_back(tau.mass());
      tau_iddecaymode.push_back(tau.tauID("decayModeFinding"));
	   tau_idisoraw.push_back(tau.tauID("byIsolationMVA3newDMwLTraw"));
      tau_idisovloose.push_back(tau.tauID("byVLooseIsolationMVA3newDMwLT"));
      tau_idisoloose.push_back(tau.tauID("byLooseIsolationMVA3newDMwLT")); 
      tau_idisomedium.push_back(tau.tauID("byMediumIsolationMVA3newDMwLT"));
      tau_idisotight.push_back(tau.tauID("byTightIsolationMVA3oldDMwLT"));
      tau_idantieletight.push_back(tau.tauID("againstElectronTightMVA5"));
      tau_idantimutight.push_back(tau.tauID("againstMuonTight3"));

      //tau_reliso_all.push_back((tau.isolationPFChargedHadrCandsPtSum() + tau.isolationPFGammaCandsEtSum()) / tau.pt());
	   numtau++;
    } 

  mtree->Fill();
  return;      
          
 
}


// ------------ method called once each job just before starting event loop  ------------
void
TauAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TauAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
