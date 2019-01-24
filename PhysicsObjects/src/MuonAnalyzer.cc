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

//classes to extract Muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h" 

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"

//class to save the histograms
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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
//declare a function to do the muon analysis
      void analyzeMuons(const edm::Event& iEvent, const edm::Handle<reco::MuonCollection> &muons);


//se declara el input tag de tipo MuonCollection         
      edm::InputTag muonInput;

	  // ----------member data ---------------------------

	int nummuon; //number of muons in the event
	TH1D *muonhisto;
	TH1D *hist_e;
	TH1D *hist_pt;
	TH1D *hist_px;
	TH1D *hist_py;
	TH1D *hist_pz;
	TH1D *hist_eta;
	TH1D *hist_phi;
	TH1D *hist_ch;
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

MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& iConfig)

{
//now do what ever initialization is needed
	edm::Service<TFileService> fs;

// se crean los muon histogramas
	hist_e = fs->make <TH1D>("hist_energy", "Muon energy", 100, 0, 5000);
	hist_pt = fs->make <TH1D>("hist_pt", "Muon pt ", 100,0,5000 );
	hist_px = fs->make <TH1D>("hist_px", "Muon px ", 100, 0, 5000 );
	hist_py = fs->make <TH1D>("hist_py", "Muon py ", 100, 0, 5000 );
	hist_pz = fs->make <TH1D>("hist_pz", "Muon pz ", 100, 0, 5000 );
	hist_eta = fs->make <TH1D>("hist_eta", "Muon eta ", 100, 0, 5000 );
	hist_phi = fs->make <TH1D>("hist_phi", "Muon phi ", 100, 0, 5000 );
	hist_ch =  fs->make <TH1D>("hist_ch", "Muon ch ", 100,0,5000 );
	muonhisto = fs->make <TH1D>("muonhisto", "Muon histo", 100, 0, 5000);

	muonInput = iConfig.getParameter<edm::InputTag>("InputCollection");

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

   analyzeMuons(iEvent,mymuons);

   mtree->Fill();
   return;

}

//************************************************************************

void 
MuonAnalyzer::analyzeMuons(const edm::Event& iEvent, const edm::Handle<reco::MuonCollection> &muons)
{
	  nummuon = 0;
	  _e.clear();
	  _pt.clear();
	  _px.clear();
	  _py.clear();
	  _pz.clear();
	  _eta.clear();
	  _phi.clear();
	  _ch.clear();

  if(muons.isValid()){
     // get the number of muons in the event
     nummuon=(*muons).size();
     muonhisto->Fill(muons->size());
        for (reco::MuonCollection::const_iterator itmuon=muons->begin(); itmuon!=muons->end(); ++itmuon){

	    _e.push_back(itmuon->energy());
	    _pt.push_back(itmuon->pt());
	    _px.push_back(itmuon->px());
	    _py.push_back(itmuon->py());
	    _pz.push_back(itmuon->pz());
	    _eta.push_back(itmuon->eta());
	    _phi.push_back(itmuon->phi());
	    _ch.push_back(itmuon->charge());

	    hist_e->Fill(itmuon->energy());
	    hist_pt->Fill(itmuon->pt());
	    hist_px->Fill(itmuon->px());
	    hist_py->Fill(itmuon->py());
	    hist_pz->Fill(itmuon->pz());
	    hist_eta->Fill(itmuon->eta());
	    hist_phi->Fill(itmuon->phi());
	    hist_ch->Fill(itmuon->charge());

        }
  }
}

//*************************************************************************


// ------------ method called once each job just before starting event loop  ------------
void
MuonAnalyzer::beginJob()
{

mfile = new TFile("MuonInfo.root","RECREATE");
mtree = new TTree("mtree","Muon information");

  mytree->Branch("muon_e",&_e);
  mytree->Branch("muon_pt",&_pt);
  mytree->Branch("muon_px",&_px);
  mytree->Branch("muon_py",&_py);
  mytree->Branch("muon_pz",&_pz);
  mytree->Branch("muon_eta",&_eta);
  mytree->Branch("muon_phi",&_phi);
  mytree->Branch("muon_ch",&_ch);

}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonAnalyzer::endJob()
{
  mfile->Write();
}

// ------------ method called when starting to processes a run  ------------
void
MuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
MuonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void
MuonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
MuonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

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
