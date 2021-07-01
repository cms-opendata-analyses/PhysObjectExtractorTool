// -*- C++ -*-
//
// Package:    FWLJMET/WeightAnalyzerBEff
// Class:      WeightAnalyzerBEff
//
/**\class WeightAnalyzerBEff WeightAnalyzerBEff.cc FWLJMET/WeightAnalyzerBEff/plugins/WeightAnalyzerBEff.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rizki Syarif
//         Created:  Fri, 10 May 2019 18:16:01 GMT
//
//


// system include files
#include <memory>
//#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

//#include "LHAPDF/LHAPDF.h"
//#include "LHAPDF/GridPDF.h"
//#include "LHAPDF/Info.h"
//#include "LHAPDF/Config.h"
//#include "LHAPDF/PDFInfo.h"
//#include "LHAPDF/PDFSet.h"
//#include "LHAPDF/Factories.h"



//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class WeightAnalyzerBEff : public edm::EDAnalyzer {
public:
 explicit WeightAnalyzerBEff(const edm::ParameterSet &);
  ~WeightAnalyzerBEff();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  TH1D* BEff_Dptbins_b;   
  TH1D* BEff_Dptbins_c;   
  TH1D* BEff_Dptbins_udsg;
  TH1D* BEffTight_Nptbins_b;  
  TH1D* BEffTight_Nptbins_c;   
  TH1D* BEffTight_Nptbins_udsg;
  TH1D* BEffMed_Nptbins_b;  
  TH1D* BEffMed_Nptbins_c;   
  TH1D* BEffMed_Nptbins_udsg;
  TH1D* BEffLoose_Nptbins_b;  
  TH1D* BEffLoose_Nptbins_c;   
  TH1D* BEffLoose_Nptbins_udsg;
  
  // ----------member data ---------------------------
//  edm::EDGetTokenT<pat::JetCollection> JETtokenAK4;  
//  edm::EDGetTokenT<GenEventInfoProduct> GEIPtoken;

const edm::InputTag GEIPtoken;


  const edm::InputTag jetsTagAK4;
  const double  discriminatorValueT;
  const double  discriminatorValueM;
  const double  discriminatorValueL;

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
WeightAnalyzerBEff::WeightAnalyzerBEff(const edm::ParameterSet& iConfig):
  jetsTagAK4(iConfig.getParameter<edm::InputTag>("jetTag")),
  discriminatorValueT(iConfig.getParameter<double>("DiscriminatorValueTight")),
  discriminatorValueM(iConfig.getParameter<double>("DiscriminatorValueMedium")),
  discriminatorValueL(iConfig.getParameter<double>("DiscriminatorValueLoose"))
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;  

 // JETtokenAK4 = consumes<pat::JetCollection>(jetsTagAK4);

	//getlabel change here?

  edm::InputTag GEIPtag("generator");
  //GEIPtoken = consumes<GenEventInfoProduct>(GEIPtag);

  //Make histograms to save counts
  double ptbinsB[11] = {0, 25, 50, 75, 100, 125, 150, 200, 400, 800, 1200};
  BEff_Dptbins_b    = fs->make<TH1D>("BEff_Dptbins_b   ","",10,ptbinsB); BEff_Dptbins_b->Sumw2();
  BEff_Dptbins_c    = fs->make<TH1D>("BEff_Dptbins_c   ","",10,ptbinsB); BEff_Dptbins_c->Sumw2();
  BEff_Dptbins_udsg = fs->make<TH1D>("BEff_Dptbins_udsg","",10,ptbinsB); BEff_Dptbins_udsg->Sumw2();
  BEffTight_Nptbins_b      = fs->make<TH1D>("BEffTight_Nptbins_b     ","",10,ptbinsB); BEffTight_Nptbins_b->Sumw2();
  BEffTight_Nptbins_c      = fs->make<TH1D>("BEffTight_Nptbins_c     ","",10,ptbinsB); BEffTight_Nptbins_c->Sumw2();
  BEffTight_Nptbins_udsg   = fs->make<TH1D>("BEffTight_Nptbins_udsg  ","",10,ptbinsB); BEffTight_Nptbins_udsg->Sumw2();
  BEffMed_Nptbins_b      = fs->make<TH1D>("BEffMed_Nptbins_b     ","",10,ptbinsB); BEffMed_Nptbins_b->Sumw2();
  BEffMed_Nptbins_c      = fs->make<TH1D>("BEffMed_Nptbins_c     ","",10,ptbinsB); BEffMed_Nptbins_c->Sumw2();
  BEffMed_Nptbins_udsg   = fs->make<TH1D>("BEffMed_Nptbins_udsg  ","",10,ptbinsB); BEffMed_Nptbins_udsg->Sumw2();
  BEffLoose_Nptbins_b      = fs->make<TH1D>("BEffLoose_Nptbins_b     ","",10,ptbinsB); BEffLoose_Nptbins_b->Sumw2();
  BEffLoose_Nptbins_c      = fs->make<TH1D>("BEffLoose_Nptbins_c     ","",10,ptbinsB); BEffLoose_Nptbins_c->Sumw2();
  BEffLoose_Nptbins_udsg   = fs->make<TH1D>("BEffLoose_Nptbins_udsg  ","",10,ptbinsB); BEffLoose_Nptbins_udsg->Sumw2();
 //15 was original instead of 9
}


WeightAnalyzerBEff::~WeightAnalyzerBEff()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
WeightAnalyzerBEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  
  // Basic +/- event weight
  Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByLabel("generator", genEvtInfo);  

  double weight = 1.0;
  if(genEvtInfo->weight() < 0) weight = -1;


  Handle<std::vector<pat::Jet>> jetsAK4;
  iEvent.getByLabel(jetsTagAK4, jetsAK4);

  for(std::vector<pat::Jet>::const_iterator it = jetsAK4->begin(); it != jetsAK4->end(); ++it){
    double disc = it->bDiscriminator("combinedSecondaryVertexBJetTags");
    int hadronFlavor = it->partonFlavour();

      if( abs(hadronFlavor)==5 ){
      BEff_Dptbins_b->Fill(it->pt(),weight);      
      if( disc >= discriminatorValueT) BEffTight_Nptbins_b->Fill(it->pt(),weight);
      if( disc >= discriminatorValueM) BEffMed_Nptbins_b->Fill(it->pt(),weight);
      if( disc >= discriminatorValueL) BEffLoose_Nptbins_b->Fill(it->pt(),weight);
      
    }
    else if( abs(hadronFlavor)==4 ){
      BEff_Dptbins_c->Fill(it->pt(),weight);
      
      if( disc >= discriminatorValueT) BEffTight_Nptbins_c->Fill(it->pt(),weight);
      if( disc >= discriminatorValueM) BEffMed_Nptbins_c->Fill(it->pt(),weight); if( disc >= discriminatorValueL) BEffLoose_Nptbins_c->Fill(it->pt(),weight);
    }
    else{
      BEff_Dptbins_udsg->Fill(it->pt(),weight);
      
      if( disc >= discriminatorValueT) BEffTight_Nptbins_udsg->Fill(it->pt(),weight);
      if( disc >= discriminatorValueM) BEffMed_Nptbins_udsg->Fill(it->pt(),weight);
      if( disc >= discriminatorValueL) BEffLoose_Nptbins_udsg->Fill(it->pt(),weight);
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void
WeightAnalyzerBEff::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
WeightAnalyzerBEff::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
WeightAnalyzerBEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WeightAnalyzerBEff);
