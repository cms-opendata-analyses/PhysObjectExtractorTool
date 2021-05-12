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
//declare a function to do the tau analysis
      void analyzeTaus(const edm::Event& iEvent, const edm::Handle<reco::PFTauCollection> &taus);

//declare the input tag for PFTauCollection
      edm::InputTag tauInput;

	  // ----------member data ---------------------------

	int numtau; //number of taus in the event

	TFile *mfile;
	TTree *mtree;

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

   analyzeTaus(iEvent,mytaus);

   mtree->Fill();
   return;
}

void
TauAnalyzer::analyzeTaus(const edm::Event& iEvent, const edm::Handle<reco::PFTauCollection> &taus)
{
   using namespace edm;
   using namespace std;


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

  if(taus.isValid()){
     // get the number of taus in the event
     numtau=taus->size();
     const float tau_min_pt = 15;
        for (reco::PFTauCollection::const_iterator itTau=taus->begin(); itTau!=taus->end(); ++itTau){
          if (itTau->pt() > tau_min_pt) {
    	        tau_e.push_back(itTau->energy());
    	        tau_pt.push_back(itTau->pt());
    	        tau_px.push_back(itTau->px());
    	        tau_py.push_back(itTau->py());
    	        tau_pz.push_back(itTau->pz());
    	        tau_eta.push_back(itTau->eta());
    	        tau_phi.push_back(itTau->phi());
    	        tau_ch.push_back(itTau->charge());
       	        //tau_decaymode.push_back(itTau->decayMode());
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


       	const auto idx = itTau - taus->begin();
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
}

// ------------ method called once each job just before starting event loop  ------------
void
TauAnalyzer::beginJob()
{

mfile = new TFile("TauInfo.root","RECREATE");
mtree = new TTree("mtree","Tau information");

  mtree->Branch("numbertau",&numtau);
  mtree->Branch("tau_e",&tau_e);
  mtree->Branch("tau_pt",&tau_pt);
  mtree->Branch("tau_px",&tau_px);
  mtree->Branch("tau_py",&tau_py);
  mtree->Branch("tau_pz",&tau_pz);
  mtree->Branch("tau_eta",&tau_eta);
  mtree->Branch("tau_phi",&tau_phi);
  mtree->Branch("tau_ch",&tau_ch);
  mtree->Branch("tau_mass",&tau_mass);
  mtree->Branch("tau_decaymode",&tau_decaymode);
  mtree->Branch("tau_reliso_all",&tau_reliso_all);
  mtree->Branch("tau_genpartidx",&tau_genpartidx);
  mtree->Branch("tau_jetidx",&tau_jetidx);
  mtree->Branch("tau_iddecaymode",&tau_iddecaymode);
  mtree->Branch("tau_idisoraw",&tau_idisoraw);
  mtree->Branch("tau_idisovloose",&tau_idisovloose);
  mtree->Branch("tau_idisoloose",&tau_idisoloose);
  mtree->Branch("tau_idisomedium",&tau_idisomedium);
  mtree->Branch("tau_idisotight",&tau_idisotight);
  mtree->Branch("tau_idantieleloose",&tau_idantieleloose);
  mtree->Branch("tau_idantielemedium",&tau_idantielemedium);
  mtree->Branch("tau_idantieletight",&tau_idantieletight);
  mtree->Branch("tau_idantimuloose",&tau_idantimuloose);
  mtree->Branch("tau_idantimumedium",&tau_idantimumedium);
  mtree->Branch("tau_idantimutight",&tau_idantimutight); 

}

// ------------ method called once each job just after ending the event loop  ------------
void
TauAnalyzer::endJob()
{
  mfile->Write();
}

// ------------ method called when starting to processes a run  ------------
void
TauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
TauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void
TauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

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
