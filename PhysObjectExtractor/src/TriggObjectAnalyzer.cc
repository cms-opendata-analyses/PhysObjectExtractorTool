// -*- C++ -*-
//
// Package:    TriggObjectAnalyzer
// Class:      TriggObjectAnalyzer
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


//classes to extract trigger information
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"


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

class TriggObjectAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TriggObjectAnalyzer(const edm::ParameterSet&);
      ~TriggObjectAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);

      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
//declare a function to do the trigger analysis
      void analyzeTriggObject(const edm::Event& iEvent, const edm::Handle<trigger::TriggerEvent> &trigEvent);


               
      //edm::InputTag electronInput;

	  // ----------member data ---------------------------

	int numtrigobj; //number of trigger objects in the event
	TH1D *trigobjhisto;
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

	std::vector<float> obj_e;
  	std::vector<float> obj_pt;
  	std::vector<float> obj_px;
  	std::vector<float> obj_py;
  	std::vector<float> obj_pz;
  	std::vector<float> obj_eta;
  	std::vector<float> obj_phi;
  	std::vector<float> obj_ch;
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

TriggObjectAnalyzer::TriggObjectAnalyzer(const edm::ParameterSet& iConfig)

{
//now do what ever initialization is needed
	edm::Service<TFileService> fs;

// se crean los electron histogramas
	hist_e = fs->make <TH1D>("hist_energy", "obj energy", 100, 0, 5000);
	hist_pt = fs->make <TH1D>("hist_pt", "obj pt ", 100,0,5000 );
	hist_px = fs->make <TH1D>("hist_px", "obj px ", 100, 0, 5000 );
	hist_py = fs->make <TH1D>("hist_py", "obj py ", 100, 0, 5000 );
	hist_pz = fs->make <TH1D>("hist_pz", "obj pz ", 100, 0, 5000 );
	hist_eta = fs->make <TH1D>("hist_eta", "obj eta ", 100, 0, 5000 );
	hist_phi = fs->make <TH1D>("hist_phi", "obj phi ", 100, 0, 5000 );
	hist_ch =  fs->make <TH1D>("hist_ch", "obj ch ", 100,0,5000 );
	trigobjhisto = fs->make <TH1D>("trigobjchisto", "obj histo", 100, 0, 5000);

	//electronInput = iConfig.getParameter<edm::InputTag>("InputCollection");

}


TriggObjectAnalyzer::~TriggObjectAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void 
TriggObjectAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

	InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); //make sure have correct process on MC
	//data process=HLT, MC depends, Spring11 is REDIGI311X
	Handle<trigger::TriggerEvent> mytrigEvent; 
	iEvent.getByLabel(trigEventTag,mytrigEvent);
	analyzeTriggObject(iEvent,mytrigEvent);

    mtree->Fill();
   return;

}


void 
TriggObjectAnalyzer::analyzeTriggObject(const edm::Event& iEvent, const edm::Handle<trigger::TriggerEvent> &trigEvent)
{
	  numtrigobj = 0;
	  obj_e.clear();
	  obj_pt.clear();
	  obj_px.clear();
	  obj_py.clear();
	  obj_pz.clear();
	  obj_eta.clear();
	  obj_phi.clear();
	  obj_ch.clear();

string filterName("hltSingleJet190Regional"); 

//it is important to specify the right HLT process for the filter, not doing this is a common bug
trigger::size_type filterIndex = trigEvent->filterIndex(InputTag(filterName,"",trigEventTag.process())); 
if(filterIndex<trigEvent->sizeFilters()){ 
    const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
	numtrigobj=(*trigkeys).size();
     trigobjhisto->Fill(trigkeys->size());
    //now loop of the trigger objects passing filter
    for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
      const trigger::TriggerObject& obj = trigObjColl[*keyIt];
      //do what you want with the trigger objects, you have
      //eta,phi,pt,mass,p,px,py,pz,et,energy accessors
	    obj_e.push_back(obj->energy());
	    obj_pt.push_back(obj->pt());
	    obj_px.push_back(obj->px());
	    obj_py.push_back(obj->py());
	    obj_pz.push_back(obj->pz());
	    obj_eta.push_back(obj->eta());
	    obj_phi.push_back(obj->phi());
	    obj_ch.push_back(iobj->charge());

	    hist_e->Fill(obj->energy());
	    hist_pt->Fill(obj->pt());
	    hist_px->Fill(obj->px());
	    hist_py->Fill(obj->py());
	    hist_pz->Fill(obj->pz());
	    hist_eta->Fill(obj->eta());
	    hist_phi->Fill(obj->phi());
	    hist_ch->Fill(obj->charge());
    }
}//end filter size check

}

//*************************************************************************


// ------------ method called once each job just before starting event loop  ------------
void
TriggObjectAnalyzer::beginJob()
{

mfile = new TFile("TriggObjectInfo.root","RECREATE");
mtree = new TTree("mtree","TriggObject information");

  mtree->Branch("obj_e",&obj_e);
  mtree->Branch("obj_pt",&obj_pt);
  mtree->Branch("obj_px",&obj_px);
  mtree->Branch("obj_py",&obj_py);
  mtree->Branch("obj_pz",&obj_pz);
  mtree->Branch("obj_eta",&obj_eta);
  mtree->Branch("obj_phi",&obj_phi);
  mtree->Branch("obj_ch",&obj_ch);

}

// ------------ method called once each job just after ending the event loop  ------------
void
TriggObjectAnalyzer::endJob()
{
  mfile->Write();
}

// ------------ method called when starting to processes a run  ------------
void
TriggObjectAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
TriggObjectAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void
TriggObjectAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TriggObjectAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggObjectAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggObjectAnalyzer);
