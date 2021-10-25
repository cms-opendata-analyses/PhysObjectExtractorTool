// -*- C++ -*-
//
// Package:    Jet/JetAnalyzer
// Class:      JetAnalyzer
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

//class to extract jet information
#include "DataFormats/PatCandidates/interface/Jet.h"

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



class JetAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit JetAnalyzer(const edm::ParameterSet&);
      ~JetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

       edm::EDGetTokenT<pat::JetCollection> jetToken_;

      // ----------member data ---------------------------
      
     
  TTree *mtree;

  int numjet; //number of jets in the event
  std::vector<float> jet_e;
  std::vector<float> jet_pt;
  std::vector<float> jet_pt_uncorr;
  std::vector<float> jet_px;
  std::vector<float> jet_py;
  std::vector<float> jet_pz;
  std::vector<float> jet_eta;
  std::vector<float> jet_phi;
  std::vector<float> jet_ch;
  std::vector<float> jet_mass;
  std::vector<double> jet_btag;

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
JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig): 
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
   
{
   //now do what ever initialization is needed
   
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
  
   mtree->Branch("numberjet",&numjet);
   mtree->GetBranch("numberjet")->SetTitle("Number of Jets");
   mtree->Branch("jet_e",&jet_e);
   mtree->GetBranch("jet_e")->SetTitle("Jet Energy");
   mtree->Branch("jet_pt",&jet_pt);
   mtree->GetBranch("jet_pt")->SetTitle("Transverse Jet Momentum");
   mtree->Branch("jet_px",&jet_px);
   mtree->GetBranch("jet_px")->SetTitle("X-Component of Jet Momentum");
   mtree->Branch("jet_py",&jet_py); 
   mtree->GetBranch("jet_py")->SetTitle("Y-Component of Jet Momentum");
   mtree->Branch("jet_pz",&jet_pz);
   mtree->GetBranch("jet_pz")->SetTitle("Z-Component of Jet Momentum");
   mtree->Branch("jet_eta",&jet_eta);
   mtree->GetBranch("jet_eta")->SetTitle("Jet Eta");
   mtree->Branch("jet_phi",&jet_phi);
   mtree->GetBranch("jet_phi")->SetTitle("Jet Phi");
   mtree->Branch("jet_ch",&jet_ch);
   mtree->GetBranch("jet_ch")->SetTitle("Jet Charge");
   mtree->Branch("jet_mass",&jet_mass);
   mtree->GetBranch("jet_mass")->SetTitle("Jet Mass");
   mtree->Branch("jet_btag",&jet_btag);
   mtree->GetBranch("jet_btag")->SetTitle("Jet Btagging Discriminant (CSV)");
   mtree->Branch("jet_pt_uncorr",&jet_pt_uncorr);
   mtree->GetBranch("jet_pt_uncorr")->SetTitle("Uncorrected Transverse Jet Momentum");
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

   Handle<pat::JetCollection> jets;
   iEvent.getByToken(jetToken_, jets);

   numjet = 0;
   jet_e.clear();
   jet_pt.clear();
   jet_px.clear();
   jet_py.clear();
   jet_pz.clear();
   jet_eta.clear();
   jet_phi.clear();
   jet_ch.clear();
   jet_mass.clear();
   jet_btag.clear();
   jet_pt_uncorr.clear();

     for (const pat::Jet &jet : *jets)
    {
    jet_e.push_back(jet.energy());
	   jet_pt.push_back(jet.pt());
	   jet_px.push_back(jet.px());
	   jet_py.push_back(jet.py());
	   jet_pz.push_back(jet.pz());
	   jet_eta.push_back(jet.eta());
	   jet_phi.push_back(jet.phi());
	   jet_ch.push_back(jet.charge());
	   jet_mass.push_back(jet.mass());
     jet_btag.push_back(std::max(0.f,jet.bDiscriminator("combinedSecondaryVertexBJetTags")));
     jet_pt_uncorr.push_back(jet.pt()*jet.jecFactor(0));
	   numjet++;
    } 

  mtree->Fill();
  return;      
          
 
}


// ------------ method called once each job just before starting event loop  ------------
void
JetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
JetAnalyzer::endJob()
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
