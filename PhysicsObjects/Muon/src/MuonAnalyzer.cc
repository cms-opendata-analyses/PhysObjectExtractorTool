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

// se crean los histogramas
	hist_e = fs->make <TH1D>("hist_energy", "Electron energy", 100, 0, 5000);
	hist_pt = fs->make <TH1D>("hist_pt", "Electron pt ", 100,0,5000 );
	hist_px = fs->make <TH1D>("hist_px", "Electron px ", 100, 0, 5000 );
	hist_py = fs->make <TH1D>("hist_py", "Electron py ", 100, 0, 5000 );
	hist_pz = fs->make <TH1D>("hist_pz", "Electron pz ", 100, 0, 5000 );
	hist_eta = fs->make <TH1D>("hist_eta", "Electron eta ", 100, 0, 5000 );
	hist_phi = fs->make <TH1D>("hist_phi", "Electron phi ", 100, 0, 5000 );
	hist_ch =  fs->make <TH1D>("hist_ch", "Electron ch ", 100,0,5000 );
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

  if(objeto.isValid()){
     // get the number of muons in the event
     nummuon=(*muons).size();
     objetohisto->Fill(objeto->size());
        for (reco::MuonCollection::const_iterator itobjeto=objeto->begin(); itobjeto!=objeto->end(); ++itobjeto){

	    _e.push_back(itobjeto->energy());
	    _pt.push_back(itobjeto->pt());
	    _px.push_back(itobjeto->px());
	    _py.push_back(itobjeto->py());
	    _pz.push_back(itobjeto->pz());
	    _eta.push_back(itobjeto->eta());
	    _phi.push_back(itobjeto->phi());
	    _ch.push_back(itobjeto->charge());

	    hist_e->Fill(itobjeto->energy());
	    hist_pt->Fill(itobjeto->pt());
	    hist_px->Fill(itobjeto->px());
	    hist_py->Fill(itobjeto->py());
	    hist_pz->Fill(itobjeto->pz());
	    hist_eta->Fill(itobjeto->eta());
	    hist_phi->Fill(itobjeto->phi());
	    hist_ch->Fill(itobjeto->charge());

        }
  }
}

//*************************************************************************


// ------------ method called once each job just before starting event loop  ------------
void
MuonAnalyzer::beginJob()
{

mfile = new TFile("ObjetoInfo.root","RECREATE");
mtree = new TTree("mtree","Objeto information");

  mytree->Branch("_e",&_e);
  mytree->Branch("_pt",&_pt);
  mytree->Branch("_px",&_px);
  mytree->Branch("_py",&_py);
  mytree->Branch("_pz",&_pz);
  mytree->Branch("_eta",&_eta);
  mytree->Branch("_phi",&_phi);
  mytree->Branch("_ch",&_ch);

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
