// -*- C++ -*-
//
// Package:    JetAnalyzer
// Class:      JetAnalyzer
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//classes to extract PFJet information
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class JetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit JetAnalyzer(const edm::ParameterSet&);
      ~JetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);

      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
//declare a function to do the jet analysis
      void analyzeJets(const edm::Event& iEvent, const edm::Handle<reco::PFJetCollection> &jets);


//se declara el input tag de tipo PFJetCollection
      edm::InputTag jetInput;

	  // ----------member data ---------------------------

	int numjet; //number of jets in the event

	TFile *mfile;
	TTree *mtree;

	  std::vector<float> _e;
  	std::vector<float> _pt;
  	std::vector<float> _px;
  	std::vector<float> _py;
  	std::vector<float> _pz;
  	std::vector<float> _eta;
  	std::vector<float> _phi;
  	std::vector<float> _ch;
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

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
	jetInput = iConfig.getParameter<edm::InputTag>("InputCollection");
}

JetAnalyzer::~JetAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
JetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<reco::PFJetCollection> myjets;
   iEvent.getByLabel(jetInput, myjets);

   analyzeJets(iEvent,myjets);

   mtree->Fill();
   return;
}

void
JetAnalyzer::analyzeJets(const edm::Event& iEvent, const edm::Handle<reco::PFJetCollection> &jets)
{
	  numjet = 0;
	  _e.clear();
	  _pt.clear();
	  _px.clear();
	  _py.clear();
	  _pz.clear();
	  _eta.clear();
	  _phi.clear();
	  _ch.clear();

  if(jets.isValid()){
     // get the number of jets in the event
     numjet=(*jets).size();
        for (reco::PFJetCollection::const_iterator itjet=jets->begin(); itjet!=jets->end(); ++itjet){

        	    _e.push_back(itjet->energy());
        	    _pt.push_back(itjet->pt());
        	    _px.push_back(itjet->px());
        	    _py.push_back(itjet->py());
        	    _pz.push_back(itjet->pz());
        	    _eta.push_back(itjet->eta());
        	    _phi.push_back(itjet->phi());
        	    _ch.push_back(itjet->charge());
        }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void
JetAnalyzer::beginJob()
{

mfile = new TFile("JetInfo.root","RECREATE");
mtree = new TTree("mtree","Jet information");

  mtree->Branch("jet_e",&_e);
  mtree->Branch("jet_pt",&_pt);
  mtree->Branch("jet_px",&_px);
  mtree->Branch("jet_py",&_py);
  mtree->Branch("jet_pz",&_pz);
  mtree->Branch("jet_eta",&_eta);
  mtree->Branch("jet_phi",&_phi);
  mtree->Branch("jet_ch",&_ch);
}

// ------------ method called once each job just after ending the event loop  ------------
void
JetAnalyzer::endJob()
{
  mfile->Write();
}

// ------------ method called when starting to processes a run  ------------
void
JetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
JetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void
JetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
JetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetAnalyzer);
