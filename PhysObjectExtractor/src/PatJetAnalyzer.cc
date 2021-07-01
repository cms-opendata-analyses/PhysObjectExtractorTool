// -*- C++ -*-
//
// Package:    PatJetAnalyzer
// Class:      PatJetAnalyzer
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
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "math.h"

//classes to extract PFJet information
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

#include "TRandom3.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//
// class declaration
//

class PatJetAnalyzer : public edm::EDAnalyzer {
public:
  explicit PatJetAnalyzer(const edm::ParameterSet&);
  ~PatJetAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  std::vector<float> factorLookup(float eta);
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;  
       
  // ----------member data ---------------------------    
  // jec variables
  std::vector<std::string> jecPayloadNames_;
  std::string              jecUncName_;
  std::string              jerResName_;
  boost::shared_ptr<JetCorrectionUncertainty> jecUnc_;
  boost::shared_ptr<SimpleJetCorrector> jer_;
  bool isData;

  edm::InputTag jetInput;

  int numjet; //number of jets in the event
  TTree *mtree;
  
  std::vector<float> jet_pt;
  std::vector<float> jet_eta;
  std::vector<float> jet_phi;
  std::vector<float> jet_ch;
  std::vector<float> jet_mass;
  std::vector<float> jet_btag;
  std::vector<float> corr_jet_pt;
  std::vector<float> corr_jet_ptUp;
  std::vector<float> corr_jet_ptDown;
  std::vector<float> corr_jet_ptSmearUp;
  std::vector<float> corr_jet_ptSmearDown;
  std::vector<float> corr_jet_mass;
  std::vector<float> corr_jet_e;
  std::vector<float> corr_jet_px;
  std::vector<float> corr_jet_py;
  std::vector<float> corr_jet_pz;
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

PatJetAnalyzer::PatJetAnalyzer(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
  jetInput = iConfig.getParameter<edm::InputTag>("InputCollection");
  isData = iConfig.getParameter<bool>("isData");
  edm::Service<TFileService> fs;
  mtree = fs->make<TTree>("Events", "Events");
  
  jecUncName_ = iConfig.getParameter<edm::FileInPath>("jecUncName").fullPath(); // JEC uncertainties
  jerResName_ = iConfig.getParameter<edm::FileInPath>("jerResName").fullPath(); // JER Resolutions                               

  // Make the FactorizedJetCorrector and Uncertainty
  jecUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecUncName_) );
  JetCorrectorParameters *jerPar = new JetCorrectorParameters(jerResName_);
  jer_ = boost::shared_ptr<SimpleJetCorrector>( new SimpleJetCorrector(*jerPar) );  

  mtree->Branch("numberjet",&numjet);
  mtree->GetBranch("numberjet")->SetTitle("Number of Jets");
  mtree->Branch("jet_pt",&jet_pt);
  mtree->GetBranch("jet_pt")->SetTitle("Uncorrected Transverse Jet Momentum");
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
  mtree->Branch("corr_jet_pt",&corr_jet_pt);
  mtree->GetBranch("corr_jet_pt")->SetTitle("Corrected Transverse Jet Momentum");
  mtree->Branch("corr_jet_ptUp",&corr_jet_ptUp);
  mtree->GetBranch("corr_jet_ptUp")->SetTitle("Corrected Transverse Jet Momentum (JEC Shifted Up)");
  mtree->Branch("corr_jet_ptDown",&corr_jet_ptDown);
  mtree->GetBranch("corr_jet_ptDown")->SetTitle("Corrected Transverse Jet Momentum (JEC Shifted Down)");
  mtree->Branch("corr_jet_ptSmearUp",&corr_jet_ptSmearUp);
  mtree->GetBranch("corr_jet_ptSmearUp")->SetTitle("Corrected Transverse Jet Momentum (JER Shifted Up)");
  mtree->Branch("corr_jet_ptSmearDown",&corr_jet_ptSmearDown);	
  mtree->GetBranch("corr_jet_ptSmearDown")->SetTitle("Corrected Transverse Jet Momentum (JER Shifted Down)");
  mtree->Branch("corr_jet_mass",&corr_jet_mass);  
  mtree->GetBranch("corr_jet_mass")->SetTitle("Corrected Jet Mass");
  mtree->Branch("corr_jet_e",&corr_jet_e);
  mtree->GetBranch("corr_jet_e")->SetTitle("Corrected Jet Energy");
  mtree->Branch("corr_jet_px",&corr_jet_px);
  mtree->GetBranch("corr_jet_px")->SetTitle("Corrected X-Component of Jet Momentum");
  mtree->Branch("corr_jet_py",&corr_jet_py);
  mtree->GetBranch("corr_jet_py")->SetTitle("Corrected Y-Component of Jet Momentum");
  mtree->Branch("corr_jet_pz",&corr_jet_pz);
  mtree->GetBranch("corr_jet_pz")->SetTitle("Corrected Z-Component of Jet Momentum");
}

PatJetAnalyzer::~PatJetAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
PatJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<std::vector<pat::Jet>> myjets;
   iEvent.getByLabel(jetInput, myjets);
   Handle<reco::JetTagCollection> btags;
   iEvent.getByLabel(InputTag("combinedSecondaryVertexBJetTags"), btags);
   Handle<double> rhoHandle;
   iEvent.getByLabel(InputTag("fixedGridRhoAll"), rhoHandle);
   Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
   numjet = 0;
   jet_pt.clear();
   jet_eta.clear();
   jet_phi.clear();
   jet_ch.clear();
   jet_mass.clear();
   jet_btag.clear();
   corr_jet_pt.clear();
   corr_jet_ptUp.clear();
   corr_jet_ptDown.clear();
   corr_jet_ptSmearUp.clear();
   corr_jet_ptSmearDown.clear();
   corr_jet_mass.clear();
   corr_jet_e.clear();
   corr_jet_px.clear();
   corr_jet_py.clear();
   corr_jet_pz.clear();

   if(myjets.isValid()){
     // get the number of jets in the event
     numjet=myjets->size();
     for (std::vector<pat::Jet>::const_iterator itjet=myjets->begin(); itjet!=myjets->end(); ++itjet){
       pat::Jet uncorrJet = itjet->correctedJet(0);     

       double corrUp = 1.0;
       double corrDown = 1.0;
       jecUnc_->setJetEta( itjet->eta() );
       jecUnc_->setJetPt( itjet->pt() );
       corrUp = (1 + fabs(jecUnc_->getUncertainty(1)));
       jecUnc_->setJetEta( itjet->eta() );
       jecUnc_->setJetPt( itjet->pt() );
       corrDown = (1 - fabs(jecUnc_->getUncertainty(-1)));
       
       float ptscale = 1;
       float ptscale_down = 1;
       float ptscale_up = 1;
       float res = 1;
       if(!isData) {
         std::vector<float> factors = factorLookup(fabs(itjet->eta())); // returns in order {factor, factor_down, factor_up}
         std::vector<float> feta;
         std::vector<float> PTNPU;
         feta.push_back( fabs(itjet->eta()) );
         PTNPU.push_back( itjet->pt() );
         PTNPU.push_back( vertices->size() );

         res = jer_->correction(feta, PTNPU);
         float pt = itjet->pt();
	 const reco::GenJet *genJet = itjet->genJet();
	 bool smeared = false;
	 if(genJet){
	   double deltaPt = fabs(genJet->pt() - pt);
           double deltaR = reco::deltaR(genJet->p4(),itjet->p4());
	   if ((deltaR < 0.2) && deltaPt <= 3*pt*res){
	     double gen_pt = genJet->pt();
	     double reco_pt = pt;
             double deltapt = (reco_pt - gen_pt) * factors[0];
             double deltapt_down = (reco_pt - gen_pt) * factors[1];
             double deltapt_up = (reco_pt - gen_pt) * factors[2];
	     ptscale = max(0.0, (reco_pt + deltapt) / reco_pt);
	     ptscale_up = max(0.0, (reco_pt + deltapt_up) / reco_pt);
	     ptscale_down = max(0.0, (reco_pt + deltapt_down) / reco_pt);
	     smeared = true;
	   }
	 } 
         if (!smeared && factors[0]>0) {
            TRandom3 JERrand;

            JERrand.SetSeed(abs(static_cast<int>(itjet->phi()*1e4)));
            ptscale = max(0.0, JERrand.Gaus(itjet->pt(),sqrt(factors[0]*(factors[0]+2))*res*itjet->pt())/itjet->pt());
			
	    JERrand.SetSeed(abs(static_cast<int>(itjet->phi()*1e4)));
            ptscale_down = max(0.0, JERrand.Gaus(itjet->pt(),sqrt(factors[1]*(factors[1]+2))*res*itjet->pt())/itjet->pt());
 
            JERrand.SetSeed(abs(static_cast<int>(itjet->phi()*1e4)));
            ptscale_up = max(0.0, JERrand.Gaus(itjet->pt(),sqrt(factors[2]*(factors[2]+2))*res*itjet->pt())/itjet->pt());
         }
       }
       jet_pt.push_back(uncorrJet.pt());
       jet_eta.push_back(itjet->eta());
       jet_phi.push_back(itjet->phi());
       jet_ch.push_back(itjet->charge());
       jet_mass.push_back(uncorrJet.mass());
       jet_btag.push_back(itjet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
       corr_jet_pt.push_back(ptscale*itjet->pt());
       corr_jet_ptUp.push_back(ptscale*corrUp*itjet->pt());
       corr_jet_ptDown.push_back(ptscale*corrDown*itjet->pt());
       corr_jet_ptSmearUp.push_back(ptscale_up*itjet->pt());
       corr_jet_ptSmearDown.push_back(ptscale_down*itjet->pt()); 
       corr_jet_mass.push_back(itjet->mass());
       corr_jet_e.push_back(itjet->energy());
       corr_jet_px.push_back(itjet->px());
       corr_jet_py.push_back(itjet->py());
       corr_jet_pz.push_back(itjet->pz());
       }
   }
  mtree->Fill();
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void
PatJetAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
PatJetAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
PatJetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
PatJetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
PatJetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
PatJetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PatJetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
std::vector<float>
PatJetAnalyzer::factorLookup(float eta) { //used in jet loop for JER factor value
  if(eta > 3.2) { //input is > 0
    return {1.056, .865, 1.247}; // {factor, factor_down, factor_up}
  }
  else if(eta > 2.8) {
    return {1.395, 1.332, 1.468};
  }
  else if(eta > 2.3) {
    return {1.254, 1.192, 1.316};
  }
  else if(eta > 1.7) {
    return {1.208, 1.162, 1.254};
  }
  else if(eta > 1.1) {
    return {1.121, 1.092, 1.15};
  }
  else if(eta > .5) {
    return {1.099, 1.071, 1.127};
  }
  else {
    return {1.079, 1.053, 1.105};
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatJetAnalyzer);
