// -*- C++ -*-
//
// Package:    FatjetAnalyzer
// Class:      FatjetAnalyzer
//

// system include files
#include <memory>
#include <TMath.h>
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

//class to extract jet information
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
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

class FatjetAnalyzer : public edm::EDAnalyzer {
public:
  explicit FatjetAnalyzer(const edm::ParameterSet&);
  ~FatjetAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;  
       
  //declare the input tag for PFJetCollection
  edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  
  // ----------member data ---------------------------    
  // jec variables
  std::vector<std::string> jecPayloadNames_;
  std::string              jecL1_;
  std::string              jecL2_;
  std::string              jecL3_;
  std::string              jecRes_;
  std::string              jetJECUncName_;
  std::string              jetResName_;
  std::string              sfName_;
  boost::shared_ptr<JetCorrectionUncertainty> jecUnc_;
  boost::shared_ptr<FactorizedJetCorrector> jec_;
  boost::shared_ptr<FactorizedJetCorrector> jecL2L3_;
  bool isData;

  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor resolution_sf;

  int numfatjet; //number of jets in the event
  TTree *mtree;
  std::vector<float> fatjet_pt;
  std::vector<float> fatjet_eta;
  std::vector<float> fatjet_phi;
  std::vector<float> fatjet_ch;
  std::vector<float> fatjet_mass;
  std::vector<float> corr_fatjet_pt;
  std::vector<float> corr_fatjet_ptUp;
  std::vector<float> corr_fatjet_ptDown;
  std::vector<float> corr_fatjet_ptSmearUp;
  std::vector<float> corr_fatjet_ptSmearDown;
  std::vector<float> corr_fatjet_mass;
  std::vector<float> corr_fatjet_e;
  std::vector<float> corr_fatjet_px;
  std::vector<float> corr_fatjet_py;
  std::vector<float> corr_fatjet_pz;
  std::vector<float> fatjet_prunedmass;
  std::vector<float> fatjet_softdropmass;
  std::vector<float> fatjet_tau1;
  std::vector<float> fatjet_tau2;
  std::vector<float> fatjet_tau3;
  std::vector<float> fatjet_subjet1_btag;
  std::vector<float> fatjet_subjet2_btag;
  std::vector<float> fatjet_subjet1_hflav;
  std::vector<float> fatjet_subjet2_hflav;

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

FatjetAnalyzer::FatjetAnalyzer(const edm::ParameterSet& iConfig):
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  jecL1_(iConfig.getParameter<edm::FileInPath>("jecL1Name").fullPath()), // JEC level payloads for test                  
  jecL2_(iConfig.getParameter<edm::FileInPath>("jecL2Name").fullPath()), // JEC level payloads for test            
  jecL3_(iConfig.getParameter<edm::FileInPath>("jecL3Name").fullPath()), // JEC level payloads for test            
  jecRes_(iConfig.getParameter<edm::FileInPath>("jecResName").fullPath()),
  jetJECUncName_(iConfig.getParameter<edm::FileInPath>("jetJECUncName").fullPath()), // JEC uncertainties
  jetResName_(iConfig.getParameter<edm::FileInPath>("jerResName").fullPath()), // JER Resolutions
  sfName_(iConfig.getParameter<edm::FileInPath>("jerSFName").fullPath()), // JER Resolutions
  isData(iConfig.getParameter<bool>("isData"))
{
//now do what ever initialization is needed
  edm::Service<TFileService> fs;
  mtree = fs->make<TTree>("Events", "Events");

  //Get the factorized jet corrector parameters.
  jecPayloadNames_.push_back(jecL1_);
  jecPayloadNames_.push_back(jecL2_);
  jecPayloadNames_.push_back(jecL3_);
  if( isData == true ) jecPayloadNames_.push_back(jecRes_);
  std::vector<JetCorrectorParameters> vPar;
  std::vector<JetCorrectorParameters> vParL2L3;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNames_.begin(),
	  payloadEnd = jecPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vPar.push_back(pars);
    if ((*ipayload).find("L1FastJet") == std::string::npos) vParL2L3.push_back(pars);
  }

  // Make the FactorizedJetCorrector and Uncertainty                                                                                              
  jec_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  jecL2L3_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vParL2L3) );
  jecUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jetJECUncName_) );

  resolution = JME::JetResolution(jetResName_);
  resolution_sf = JME::JetResolutionScaleFactor(sfName_);

  edm::InputTag rhotag("fixedGridRhoFastjetAll");
  rhoToken_ = consumes<double>(rhotag);
  edm::InputTag pvtag("offlineSlimmedPrimaryVertices");
  pvToken_ = consumes<reco::VertexCollection>(pvtag);
	
  mtree->Branch("numberfatjet",&numfatjet);
  mtree->GetBranch("numberfatjet")->SetTitle("Number of Fatjets");
  mtree->Branch("fatjet_pt",&fatjet_pt);
  mtree->GetBranch("fatjet_pt")->SetTitle("Uncorrected Transverse Fatjet Momentum");
  mtree->Branch("fatjet_eta",&fatjet_eta);
  mtree->GetBranch("fatjet_eta")->SetTitle("Fatjet Eta");
  mtree->Branch("fatjet_phi",&fatjet_phi);
  mtree->GetBranch("fatjet_phi")->SetTitle("Fatjet Phi");
  mtree->Branch("fatjet_ch",&fatjet_ch);
  mtree->GetBranch("fatjet_ch")->SetTitle("Fatjet Charge");
  mtree->Branch("fatjet_mass",&fatjet_mass);
  mtree->GetBranch("fatjet_mass")->SetTitle("Fatjet Mass");
  mtree->Branch("corr_fatjet_pt",&corr_fatjet_pt);
  mtree->GetBranch("corr_fatjet_pt")->SetTitle("Corrected Transverse Fatjet Momentum");
  mtree->Branch("corr_fatjet_ptUp",&corr_fatjet_ptUp);
  mtree->GetBranch("corr_fatjet_ptUp")->SetTitle("Corrected Transverse Fatjet Momentum (JEC Shifted Up)");
  mtree->Branch("corr_fatjet_ptDown",&corr_fatjet_ptDown);
  mtree->GetBranch("corr_fatjet_ptDown")->SetTitle("Corrected Transverse Fatjet Momentum (JEC Shifted Down)");
  mtree->Branch("corr_fatjet_ptSmearUp",&corr_fatjet_ptSmearUp);
  mtree->GetBranch("corr_fatjet_ptSmearUp")->SetTitle("Corrected Transverse Fatjet Momentum (JER Shifted Up)");
  mtree->Branch("corr_fatjet_ptSmearDown",&corr_fatjet_ptSmearDown);	
  mtree->GetBranch("corr_fatjet_ptSmearDown")->SetTitle("Corrected Transverse Fatjet Momentum (JER Shifted Down)");
  mtree->Branch("corr_fatjet_mass",&corr_fatjet_mass);  
  mtree->GetBranch("corr_fatjet_mass")->SetTitle("Corrected Fatjet Mass");
  mtree->Branch("corr_fatjet_e",&corr_fatjet_e);
  mtree->GetBranch("corr_fatjet_e")->SetTitle("Corrected Fatjet Energy");
  mtree->Branch("corr_fatjet_px",&corr_fatjet_px);
  mtree->GetBranch("corr_fatjet_px")->SetTitle("Corrected X-Component of Fatjet Momentum");
  mtree->Branch("corr_fatjet_py",&corr_fatjet_py);
  mtree->GetBranch("corr_fatjet_py")->SetTitle("Corrected Y-Component of Fatjet Momentum");
  mtree->Branch("corr_fatjet_pz",&corr_fatjet_pz);
  mtree->GetBranch("corr_fatjet_pz")->SetTitle("Corrected Z-Component of Fatjet Momentum");
  mtree->Branch("fatjet_prunedmass",&fatjet_prunedmass);
  mtree->GetBranch("fatjet_prunedmass")->SetTitle("L2+L3-corrected pruned mass of Fatjet");
  mtree->Branch("fatjet_softdropmass",&fatjet_softdropmass);
  mtree->GetBranch("fatjet_softdropmass")->SetTitle("L2+L3-corrected softdrop mass of Fatjet");
  mtree->Branch("fatjet_tau1",&fatjet_tau1);
  mtree->GetBranch("fatjet_tau1")->SetTitle("N-subjettiness tau_1 of Fatjet");
  mtree->Branch("fatjet_tau2",&fatjet_tau2);
  mtree->GetBranch("fatjet_tau2")->SetTitle("N-subjettiness tau_2 of Fatjet");
  mtree->Branch("fatjet_tau3",&fatjet_tau3);
  mtree->GetBranch("fatjet_tau3")->SetTitle("N-subjettiness tau_3 of Fatjet");
  mtree->Branch("fatjet_subjet1_btag",&fatjet_subjet1_btag);
  mtree->GetBranch("fatjet_subjet1_btag")->SetTitle("Leading softdrop subjet 1 b discriminant");
  mtree->Branch("fatjet_subjet2_btag",&fatjet_subjet2_btag);
  mtree->GetBranch("fatjet_subjet2_btag")->SetTitle("Leading softdrop subjet 2 b discriminant");
  mtree->Branch("fatjet_subjet1_hflav",&fatjet_subjet1_hflav);
  mtree->GetBranch("fatjet_subjet1_hflav")->SetTitle("Leading softdrop subjet 1 hadron flavour");
  mtree->Branch("fatjet_subjet2_hflav",&fatjet_subjet2_hflav);
  mtree->GetBranch("fatjet_subjet2_hflav")->SetTitle("Leading softdrop subjet 2 hadron flavour");
}

FatjetAnalyzer::~FatjetAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
FatjetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;

  Handle<pat::JetCollection> fatjets;
  iEvent.getByToken(fatjetToken_, fatjets);

  Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(pvToken_, vertices);

  numfatjet = 0;
  fatjet_pt.clear();
  fatjet_eta.clear();
  fatjet_phi.clear();
  fatjet_ch.clear();
  fatjet_mass.clear();
  corr_fatjet_pt.clear();
  corr_fatjet_ptUp.clear();
  corr_fatjet_ptDown.clear();
  corr_fatjet_ptSmearUp.clear();
  corr_fatjet_ptSmearDown.clear();
  corr_fatjet_mass.clear();
  corr_fatjet_e.clear();
  corr_fatjet_px.clear();
  corr_fatjet_py.clear();
  corr_fatjet_pz.clear();
  fatjet_prunedmass.clear();
  fatjet_softdropmass.clear();
  fatjet_tau1.clear();
  fatjet_tau2.clear();
  fatjet_tau3.clear();
  fatjet_subjet1_btag.clear();
  fatjet_subjet2_btag.clear();
  fatjet_subjet1_hflav.clear();
  fatjet_subjet2_hflav.clear();

  double corrpt;
  double corrUp, corrDown, corr;
  float ptscale, ptscale_down, ptscale_up;
  int min_pt = 200;
  
  if(fatjets.isValid()){

    for (const pat::Jet &fatjet : *fatjets){
      pat::Jet uncorrFatjet = fatjet.correctedJet(0);

      corr = 1;
      jec_->setJetEta( uncorrFatjet.eta() );
      jec_->setJetPt ( uncorrFatjet.pt() );
      jec_->setJetE  ( uncorrFatjet.energy() );
      jec_->setJetA  ( fatjet.jetArea() );
      jec_->setRho   ( *(rhoHandle.product()) );
      jec_->setNPV   ( vertices->size() );
      corr = jec_->getCorrection();

      corrpt = corr*fatjet.pt();
      corrUp = 1.0;
      corrDown = 1.0;

      if( fabs(fatjet.eta()) < 5) jecUnc_->setJetEta( fatjet.eta() );
      else jecUnc_->setJetEta( 4.99 );
      jecUnc_->setJetPt( corrpt );
      corrUp = (1 + fabs(jecUnc_->getUncertainty(1)));

      if( fabs(fatjet.eta()) < 5) jecUnc_->setJetEta( fatjet.eta() );
      else jecUnc_->setJetEta( 4.99 );
      jecUnc_->setJetPt( corrpt );
      corrDown = (1 - fabs(jecUnc_->getUncertainty(-1)));        

      ptscale = 1;
      ptscale_down = 1;
      ptscale_up = 1;

      if(!isData) {

	JME::JetParameters JERparameters = {{JME::Binning::JetPt, corrpt}, {JME::Binning::JetEta, fatjet.eta()}, {JME::Binning::Rho, *(rhoHandle.product())}};
	float res = resolution.getResolution(JERparameters); 
	float sf = resolution_sf.getScaleFactor(JERparameters);
	float sf_up = resolution_sf.getScaleFactor(JERparameters, Variation::UP);
	float sf_down = resolution_sf.getScaleFactor(JERparameters, Variation::DOWN);

	const reco::GenJet *genFatjet = fatjet.genJet();
	bool smeared = false;
	if(genFatjet){
	  double deltaPt = fabs(genFatjet->pt() - corrpt);
	  double deltaR = reco::deltaR(genFatjet->p4(),fatjet.p4());
	  if ((deltaR < 0.4) && deltaPt <= 3*corrpt*res){
	    ptscale = max(0.0, 1 + (sf - 1.0)*(corrpt - genFatjet->pt())/corrpt);
	    ptscale_down = max(0.0, 1 + (sf_down - 1.0)*(corrpt - genFatjet->pt())/corrpt);
	    ptscale_up = max(0.0, 1 + (sf_up - 1.0)*(corrpt - genFatjet->pt())/corrpt);
	    smeared = true;
	  }
	} 
	if (!smeared) {
	  TRandom3 JERrand;
	  
	  JERrand.SetSeed(abs(static_cast<int>(fatjet.phi()*1e4)));
	  ptscale = max(0.0, 1.0 + JERrand.Gaus(0, res)*sqrt(max(0.0, sf*sf - 1.0)));
	  
	  JERrand.SetSeed(abs(static_cast<int>(fatjet.phi()*1e4)));
	  ptscale_down = max(0.0, 1.0 + JERrand.Gaus(0, res)*sqrt(max(0.0, sf_down*sf_down - 1.0)));
	  
	  JERrand.SetSeed(abs(static_cast<int>(fatjet.phi()*1e4)));
	  ptscale_up = max(0.0, 1.0 + JERrand.Gaus(0, res)*sqrt(max(0.0, sf_up*sf_up - 1.0)));
	}

      }
      
      if( ptscale*corrpt <= min_pt) continue;

      pat::Jet recorrFatjet = uncorrFatjet;
      recorrFatjet.scaleEnergy(ptscale*corr);

      fatjet_pt.push_back(uncorrFatjet.pt());
      fatjet_eta.push_back(uncorrFatjet.eta());
      fatjet_phi.push_back(uncorrFatjet.phi());
      fatjet_ch.push_back(uncorrFatjet.charge());
      fatjet_mass.push_back(uncorrFatjet.mass());

      corr_fatjet_pt.push_back(recorrFatjet.pt());
      corr_fatjet_ptUp.push_back(corrUp*recorrFatjet.pt());
      corr_fatjet_ptDown.push_back(corrDown*recorrFatjet.pt());
      corr_fatjet_ptSmearUp.push_back(ptscale_up*recorrFatjet.pt()/ptscale);
      corr_fatjet_ptSmearDown.push_back(ptscale_down*recorrFatjet.pt()/ptscale); 
      corr_fatjet_mass.push_back(recorrFatjet.mass());
      corr_fatjet_e.push_back(recorrFatjet.energy());
      corr_fatjet_px.push_back(recorrFatjet.px());
      corr_fatjet_py.push_back(recorrFatjet.py());
      corr_fatjet_pz.push_back(recorrFatjet.pz());

      double corrL2L3 = 1;
      jecL2L3_->setJetEta( uncorrFatjet.eta() );
      jecL2L3_->setJetPt ( uncorrFatjet.pt() );
      jecL2L3_->setJetE  ( uncorrFatjet.energy() );
      jecL2L3_->setJetA  ( fatjet.jetArea() );
      jecL2L3_->setRho   ( *(rhoHandle.product()) );
      jecL2L3_->setNPV   ( vertices->size() );
      corrL2L3 = jecL2L3_->getCorrection();

      fatjet_prunedmass.push_back(corrL2L3*(double)recorrFatjet.userFloat("ak8PFJetsCHSPrunedMass"));
      fatjet_softdropmass.push_back(corrL2L3*(double)recorrFatjet.userFloat("ak8PFJetsCHSSoftDropMass"));
      fatjet_tau1.push_back((double)recorrFatjet.userFloat("NjettinessAK8:tau1"));
      fatjet_tau2.push_back((double)recorrFatjet.userFloat("NjettinessAK8:tau2"));
      fatjet_tau3.push_back((double)recorrFatjet.userFloat("NjettinessAK8:tau3"));

      auto const & sdSubjets = recorrFatjet.subjets("SoftDrop");
      int nSDSubJets = sdSubjets.size();
      if(nSDSubJets > 0){
	pat::Jet subjet1 = sdSubjets.at(0);
	fatjet_subjet1_btag.push_back(subjet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	fatjet_subjet1_hflav.push_back(subjet1.hadronFlavour());
      }else{
	fatjet_subjet1_btag.push_back(-999);
	fatjet_subjet1_hflav.push_back(-999);
      }
      if(nSDSubJets > 1){
	pat::Jet subjet2 = sdSubjets.at(1);
	fatjet_subjet2_btag.push_back(subjet2.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	fatjet_subjet2_hflav.push_back(subjet2.hadronFlavour());
      }else{
	fatjet_subjet2_btag.push_back(-999);
	fatjet_subjet2_hflav.push_back(-999);
      }

      ++numfatjet;
    }
  }
  
  mtree->Fill();
  return;
  
}

// ------------ method called once each job just before starting event loop  ------------
void
FatjetAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
FatjetAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
FatjetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
FatjetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
FatjetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
FatjetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FatjetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FatjetAnalyzer);
