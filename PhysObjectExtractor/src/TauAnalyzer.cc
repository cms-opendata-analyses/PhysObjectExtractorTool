// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      TauAnalyzer
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

//classes to extract tau information
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"


//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class TauAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TauAnalyzer(const edm::ParameterSet&);
      ~TauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      
      //declare the input tag for PFTauCollection
      edm::InputTag tauInput;
      
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
      std::vector<float> tau_reliso_all;
      std::vector<float> tau_genpartidx;
      std::vector<float> tau_jetidx;
      std::vector<float> tau_iddecaymode;
      std::vector<float> tau_idisoraw;
      std::vector<float> tau_idisovloose;
      std::vector<float> tau_idisoloose;
      std::vector<float> tau_idisomedium;
      std::vector<float> tau_idisotight;
      std::vector<float> tau_idantieleloose;
      std::vector<float> tau_idantielemedium;
      std::vector<float> tau_idantieletight;
      std::vector<float> tau_idantimuloose;
      std::vector<float> tau_idantimumedium;
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

TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
	tauInput = iConfig.getParameter<edm::InputTag>("InputCollection");
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
	mtree->Branch("tau_reliso_all",&tau_reliso_all);
	mtree->GetBranch("tau_reliso_all")->SetTitle("tau reliso all");
	mtree->Branch("tau_genpartidx",&tau_genpartidx);
	mtree->GetBranch("tau_genpartidx")->SetTitle("Index into genParticle list for MC matching to status==2 taus");
	mtree->Branch("tau_jetidx",&tau_jetidx);
	mtree->GetBranch("tau_jetidx")->SetTitle("index of the associated jet");
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
	mtree->Branch("tau_idantieleloose",&tau_idantieleloose);
	mtree->GetBranch("tau_idantieleloose")->SetTitle("tau id antieleloose");
	mtree->Branch("tau_idantielemedium",&tau_idantielemedium);
	mtree->GetBranch("tau_idantielemedium")->SetTitle("tau id antielemedium");
	mtree->Branch("tau_idantieletight",&tau_idantieletight);
	mtree->GetBranch("tau_idantieletight")->SetTitle("tau id antieletight");
	mtree->Branch("tau_idantimuloose",&tau_idantimuloose);
	mtree->GetBranch("tau_idantimuloose")->SetTitle("tau id antimuloose");
	mtree->Branch("tau_idantimumedium",&tau_idantimumedium);
	mtree->GetBranch("tau_idantimumedium")->SetTitle("tau id antimumedium");
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
   using namespace std;

   Handle<reco::PFTauCollection> mytaus;
   iEvent.getByLabel(tauInput, mytaus);
   
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
   tau_reliso_all.clear();
   tau_genpartidx.clear();
   tau_jetidx.clear();
   tau_iddecaymode.clear();
   tau_idisoraw.clear();
   tau_idisovloose.clear();
   tau_idisoloose.clear();
   tau_idisomedium.clear();
   tau_idisotight.clear();
   tau_idantieleloose.clear();
   tau_idantielemedium.clear();
   tau_idantieletight.clear();
   tau_idantimuloose.clear();
   tau_idantimumedium.clear();
   tau_idantimutight.clear();

  if(mytaus.isValid()){
     // get the number of taus in the event
     numtau=mytaus->size();

     const float tau_min_pt = 15;
        for (reco::PFTauCollection::const_iterator itTau=mytaus->begin(); itTau!=mytaus->end(); ++itTau){
           if (itTau->pt() > tau_min_pt) {
    	        tau_e.push_back(itTau->energy());
    	        tau_pt.push_back(itTau->pt());
    	        tau_px.push_back(itTau->px());
    	        tau_py.push_back(itTau->py());
    	        tau_pz.push_back(itTau->pz());
    	        tau_eta.push_back(itTau->eta());
    	        tau_phi.push_back(itTau->phi());
    	        tau_ch.push_back(itTau->charge());
    	        tau_decaymode.push_back(itTau->decayMode());
    	        tau_mass.push_back(itTau->mass());
       
  //Discriminators

  Handle<reco::PFTauDiscriminator> tausLooseIso, tausVLooseIso, tausMediumIso, tausTightIso,
                             tausDecayMode, tausLooseEleRej, tausMediumEleRej,
                             tausTightEleRej, tausLooseMuonRej, tausMediumMuonRej,
                             tausTightMuonRej, tausRawIso, tausLooseIsoMVA, tausMediumIsoMVA, tausTightIsoMVA,
                             tausLooseIso3Hits, tausMediumIso3Hits, tausTightIso3Hits;

  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
          tausDecayMode);

  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr"),
          tausRawIso);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr"),
          tausVLooseIso);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"),
          tausLooseIso);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr"),
          tausMediumIso);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr"),
          tausTightIso);

  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByLooseElectronRejection"),
          tausLooseEleRej);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByMediumElectronRejection"),
          tausMediumEleRej);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByTightElectronRejection"),
          tausTightEleRej);

  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByLooseMuonRejection"),
          tausLooseMuonRej);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByMediumMuonRejection"),
          tausMediumMuonRej);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByTightMuonRejection"),
          tausTightMuonRej);

  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),
          tausLooseIso3Hits);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits"),
          tausMediumIso3Hits);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits"),
          tausTightIso3Hits);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByLooseIsolationMVA"),
          tausLooseIsoMVA);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByMediumIsolationMVA"),
          tausMediumIsoMVA);
  iEvent.getByLabel(InputTag("hpsPFTauDiscriminationByTightIsolationMVA"),
          tausTightIsoMVA);


       	const auto idx = itTau - mytaus->begin();
		tau_iddecaymode.push_back(tausDecayMode->operator[](idx).second);
      	     	tau_idisoraw.push_back(tausRawIso->operator[](idx).second);
       	tau_idisovloose.push_back(tausVLooseIso->operator[](idx).second);
       	tau_idisoloose.push_back(tausLooseIso->operator[](idx).second);
       	tau_idisomedium.push_back(tausMediumIso->operator[](idx).second);
       	tau_idisotight.push_back(tausTightIso->operator[](idx).second);
       	tau_idantieleloose.push_back(tausLooseEleRej->operator[](idx).second);
       	tau_idantielemedium.push_back(tausMediumEleRej->operator[](idx).second);
       	tau_idantieletight.push_back(tausTightEleRej->operator[](idx).second);
       	tau_idantimuloose.push_back(tausLooseMuonRej->operator[](idx).second);
       	tau_idantimumedium.push_back(tausMediumMuonRej->operator[](idx).second);
       	tau_idantimutight.push_back(tausTightMuonRej->operator[](idx).second);
       	tau_reliso_all.push_back((itTau->isolationPFChargedHadrCandsPtSum() + itTau->isolationPFGammaCandsEtSum()) / itTau->pt());
       	tau_jetidx.push_back(-1);
       	tau_genpartidx.push_back(-1);
 
          }   	
        }
  }
  
  mtree->Fill();
  return;
   
}

// ------------ method called once each job just before starting event loop  ------------
void
TauAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
TauAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
TauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
TauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
TauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
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
