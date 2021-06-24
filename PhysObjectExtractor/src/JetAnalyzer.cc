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


//
// class declaration
//

class JetAnalyzer : public edm::EDAnalyzer {
public:
  explicit JetAnalyzer(const edm::ParameterSet&);
  ~JetAnalyzer();
  
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
       
  //declare the input tag for PFJetCollection
  edm::InputTag jetInput;
  
  // ----------member data ---------------------------    
  // jec variables
  std::vector<std::string> jecPayloadNames_;
  std::string              jecL1_;
  std::string              jecL2_;
  std::string              jecL3_;
  std::string              jecRes_;
  std::string              jecUncName_;
  std::string              jerResName_;
  boost::shared_ptr<JetCorrectionUncertainty> jecUnc_;
  boost::shared_ptr<FactorizedJetCorrector> jec_;
  boost::shared_ptr<SimpleJetCorrector> ak5PFCorrector;
  bool isData;

  int numjet; //number of jets in the event
  TTree *mtree;
  std::vector<float> jet_e;
  std::vector<float> jet_pt;
  std::vector<float> jet_px;
  std::vector<float> jet_py;
  std::vector<float> jet_pz;
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
  edm::Service<TFileService> fs;
  mtree = fs->make<TTree>("Events", "Events");
  
  isData = iConfig.getParameter<bool>("isData");
  jecL1_ = iConfig.getParameter<edm::FileInPath>("jecL1Name").fullPath(); // JEC level payloads                     
  jecL2_ = iConfig.getParameter<edm::FileInPath>("jecL2Name").fullPath(); // JEC level payloads                     
  jecL3_ = iConfig.getParameter<edm::FileInPath>("jecL3Name").fullPath(); // JEC level payloads                     
  jecRes_= iConfig.getParameter<edm::FileInPath>("jecResName").fullPath();
  jecUncName_ = iConfig.getParameter<edm::FileInPath>("jecUncName").fullPath(); // JEC uncertainties
  jerResName_ = iConfig.getParameter<edm::FileInPath>("jerResName").fullPath(); // JER Resolutions                               

  //Get the factorized jet corrector parameters.
  jecPayloadNames_.push_back(jecL1_);
  jecPayloadNames_.push_back(jecL2_);
  jecPayloadNames_.push_back(jecL3_);
  if( isData == true ) jecPayloadNames_.push_back(jecRes_);
  std::vector<JetCorrectorParameters> vPar;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNames_.begin(),
	  payloadEnd = jecPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vPar.push_back(pars);
  }

  // Make the FactorizedJetCorrector and Uncertainty                                                                                              
  jec_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  jecUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecUncName_) );
  JetCorrectorParameters *ak5PFPar = new JetCorrectorParameters(jerResName_);
  ak5PFCorrector = boost::shared_ptr<SimpleJetCorrector>( new SimpleJetCorrector(*ak5PFPar) );  

	
  mtree->Branch("numberjet",&numjet);
  mtree->GetBranch("numberjet")->SetTitle("Number of Jets");
  mtree->Branch("jet_e",&jet_e);
  mtree->GetBranch("jet_e")->SetTitle("Uncorrected Jet Energy");
  mtree->Branch("jet_pt",&jet_pt);
  mtree->GetBranch("jet_pt")->SetTitle("Uncorrected Transverse Jet Momentum");
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
   Handle<reco::JetTagCollection> btags;
   iEvent.getByLabel(InputTag("combinedSecondaryVertexBJetTags"), btags);
   Handle<double> rhoHandle;
   iEvent.getByLabel(InputTag("fixedGridRhoAll"), rhoHandle);
   Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
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
   corr_jet_pt.clear();
   corr_jet_ptUp.clear();
   corr_jet_ptDown.clear();

   if(myjets.isValid()){
     // get the number of jets in the event
     numjet=myjets->size();
     for (reco::PFJetCollection::const_iterator itjet=myjets->begin(); itjet!=myjets->end(); ++itjet){
       reco::Candidate::LorentzVector uncorrJet = itjet->p4();
       jec_->setJetEta( uncorrJet.eta() );
       jec_->setJetPt ( uncorrJet.pt() );
       jec_->setJetE  ( uncorrJet.energy() );
       jec_->setJetA  ( itjet->jetArea() );
       jec_->setRho   ( *(rhoHandle.product()) );
       jec_->setNPV   ( vertices->size() );
       double corr = jec_->getCorrection();
       
       double corrUp = 1.0;
       double corrDown = 1.0;
       jecUnc_->setJetEta( uncorrJet.eta() );
       jecUnc_->setJetPt( corr * uncorrJet.pt() );
       corrUp = corr * (1 + fabs(jecUnc_->getUncertainty(1)));
       jecUnc_->setJetEta( uncorrJet.eta() );
       jecUnc_->setJetPt( corr * uncorrJet.pt() );
       corrDown = corr * (1 - fabs(jecUnc_->getUncertainty(-1)));
       
       float ptscale, ptscale_down, ptscale_up;
   
       if(isData) {
         ptscale = 1;
         ptscale_down = 1;
         ptscale_up = 1;
       } 
       else {
         std::vector<float> factors = factorLookup(fabs(itjet->eta())); // returns in order {factor, factor_down, factor_up}
         std::vector<float> feta;
         std::vector<float> PTNPU;
	 feta.push_back( fabs(itjet->eta()) );
         PTNPU.push_back( itjet->pt() );
         PTNPU.push_back( vertices->size() );

         float res = ak5PFCorrector->correction(feta, PTNPU);

         TRandom3 JERrand;

         JERrand.SetSeed(abs(static_cast<int>(itjet->phi()*1e4)));
         ptscale = max(0.0, JERrand.Gaus(itjet->pt(),sqrt(factors[0]*(factors[0]+2))*res*itjet->pt())/itjet->pt());

         JERrand.SetSeed(abs(static_cast<int>(itjet->phi()*1e4)));
         ptscale_down = max(0.0, JERrand.Gaus(itjet->pt(),sqrt(factors[1]*(factors[1]+2))*res*itjet->pt())/itjet->pt());

         JERrand.SetSeed(abs(static_cast<int>(itjet->phi()*1e4)));
         ptscale_up = max(0.0, JERrand.Gaus(itjet->pt(),sqrt(factors[2]*(factors[2]+2))*res*itjet->pt())/itjet->pt());
       }

       jet_e.push_back(itjet->energy());
       jet_pt.push_back(itjet->pt());
       jet_px.push_back(itjet->px());
       jet_py.push_back(itjet->py());
       jet_pz.push_back(itjet->pz());
       jet_eta.push_back(itjet->eta());
       jet_phi.push_back(itjet->phi());
       jet_ch.push_back(itjet->charge());
       jet_mass.push_back(itjet->mass());
       jet_btag.push_back(btags->operator[](itjet - myjets->begin()).second);
       corr_jet_pt.push_back(ptscale*corr*uncorrJet.pt());
       corr_jet_ptUp.push_back(corrUp*uncorrJet.pt());
       corr_jet_ptDown.push_back(corrDown*uncorrJet.pt());
       corr_jet_ptSmearUp.push_back(ptscale_up*corrUp*uncorrJet.pt());
       corr_jet_ptSmearDown.push_back(ptscale_down*corrUp*uncorrJet.pt());
     }
   }
   
  mtree->Fill();
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void
JetAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
JetAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
JetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
JetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
JetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

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
std::vector<float>
JetAnalyzer::factorLookup(float eta) { //used in jet loop for JER factor value
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
DEFINE_FWK_MODULE(JetAnalyzer);
