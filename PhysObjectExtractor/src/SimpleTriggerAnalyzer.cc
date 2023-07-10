// -*- C++ -*-
//
// Package:    SimpleTriggerAnalyzer
// Class:      SimpleTriggerAnalyzer
//

// system include files
#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <boost/foreach.hpp>
#include <iterator>
#include <map>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include <cassert>

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class SimpleTriggerAnalyzer : public edm::stream::EDAnalyzer< >  {
   public:
      explicit SimpleTriggerAnalyzer(const edm::ParameterSet&);
      ~SimpleTriggerAnalyzer();
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual bool analyzeSimplePrescales(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);
      virtual void initPattern(const edm::TriggerResults & result,
                        const edm::EventSetup& iSetup,
                        const edm::TriggerNames & triggerNames);
      //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      //virtual void beginJob();
      //virtual void endJob();


   private:
      // ----------member data ---------------------------
      // from HLTEventAnalyzerAOD.h
      /// module config parameters
      std::string   processName_;
      edm::InputTag triggerResultsTag_;
      const edm::EDGetTokenT<edm::TriggerResults>   triggerResultsToken_;
      /// HLT trigger names
      edm::ParameterSetID triggerNamesID_;
      // additional class data memebers
      // these are actually the containers where we will store
      // the trigger HLTPattrmation
      edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
      //to get prescales
      HLTConfigProvider hltConfig_;
      HLTPrescaleProvider hltPrescaleProvider_;

      //inspired by https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/HLTrigger/HLTfilters/interface/HLTHighLevel.h
      // input patterns that will be expanded into trigger names
      std::vector<std::string>  HLTPatterns_;

      /// list of required HLT triggers by HLT name
      std::vector<std::string>  HLTPathsByName_;

  TTree *mtree;
  //std::map<std::string, int> trigmap;
  //specific triggers bit names for a given analysis

  //These triggers need to match containers below for consistency
  //They are input by hand for the workshop exercise
  bool trig_Ele22_eta2p1_WPLoose_Gsf;
  bool trig_IsoMu20;
  bool trig_IsoTkMu20;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

SimpleTriggerAnalyzer::SimpleTriggerAnalyzer(const edm::ParameterSet& ps):
processName_(ps.getParameter<std::string>("processName")),
triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
triggerResultsToken_(consumes<edm::TriggerResults>(triggerResultsTag_)),
triggerNamesID_(),
hltPrescaleProvider_(ps, consumesCollector(), *this)
//HLTPatterns_(ps.getParameter<std::vector<std::string> >("triggerPatterns")),
//HLTPathsByName_()
{
   //now do what ever initialization is needed
   using namespace std;
   using namespace edm;

   //the order is important to match with names in the branches
   HLTPatterns_.push_back("HLT_Ele22_eta2p1_WPLoose_Gsf_v*");
   HLTPatterns_.push_back("HLT_IsoMu20_v*");
   HLTPatterns_.push_back("HLT_IsoTkMu20_v*");

   //for storing HLTPatt
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
   //mtree->Branch("triggermap", &trigmap);
   //second  stores the multiplication acceptbit*L1ps*HLTps
   //so, if negative, it means that L1ps couldn't be found.
   //look below in the code to understand the specifics
   //mtree->GetBranch("triggermap")->SetTitle("first:name of trigger, second: acceptbit*L1ps*HLTps");
   mtree->Branch("trig_Ele22_eta2p1_WPLoose_Gsf", &trig_Ele22_eta2p1_WPLoose_Gsf);
   mtree->GetBranch("trig_Ele22_eta2p1_WPLoose_Gsf")->SetTitle("Ele22_eta2p1_WPLoose_Gsf trig acceptance bit");
   mtree->Branch("trig_IsoMu20", &trig_IsoMu20);
   mtree->GetBranch("trig_IsoMu20")->SetTitle("IsoMu20 trig acceptance bit");
   mtree->Branch("trig_IsoTkMu20", &trig_IsoTkMu20);
   mtree->GetBranch("trig_IsoTkMu20")->SetTitle("IsoTkMu20 trig acceptance bit");
}



SimpleTriggerAnalyzer::~SimpleTriggerAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
// ------------ method called when ending the processing of a run  ------------
void SimpleTriggerAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a run  ------------
void SimpleTriggerAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
//--------------------------------------------------------------------------
{
    using namespace std;
    using namespace edm;

    bool changed(true);
    hltPrescaleProvider_.init(iRun,iSetup,processName_,changed);
    if (changed){
     cout<<"HLTConfig has changed for this Run. . . "<<endl;
   }
} //------------------- beginRun()


// ------------ method called for each event  ------------
void
SimpleTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   iEvent.getByToken(triggerResultsToken_,triggerResultsHandle_);
   if (!triggerResultsHandle_.isValid()) {
       LogVerbatim("SimpleTriggerAnalyzer") << "SimpleTriggerAnalyzer::analyze: Error in getting TriggerResults product from Event!" << endl;
      return;
   }
   
   HLTConfigProvider const&  hltConfig = hltPrescaleProvider_.hltConfigProvider();

   // sanity check
   assert(triggerResultsHandle_->size()==hltConfig.size());

   //Inspired in https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/HLTrigger/HLTfilters/src/HLTHighLevel.cc
   // init the TriggerNames with the TriggerResults
   const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResultsHandle_);
   bool config_changed = false;
   if (triggerNamesID_ != triggerNames.parameterSetID()) {
      triggerNamesID_ = triggerNames.parameterSetID();
      config_changed = true;
   }

   // (re)run the initialization of the container with the trigger patterns
   // - this is the first event
   // - or the HLT table has changed
   if (config_changed) {
      initPattern(*triggerResultsHandle_, iSetup, triggerNames);
   }

   trig_Ele22_eta2p1_WPLoose_Gsf = analyzeSimplePrescales(iEvent,iSetup,HLTPathsByName_.at(0));
   trig_IsoMu20 = analyzeSimplePrescales(iEvent,iSetup,HLTPathsByName_.at(1));
   trig_IsoTkMu20= analyzeSimplePrescales(iEvent,iSetup,HLTPathsByName_.at(2));

  mtree->Fill();
  return;


}

//---------------------------Actual trigger analysis-------------
bool SimpleTriggerAnalyzer::analyzeSimplePrescales(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName)
//-----------------------------------------------------------------
{

  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  int acc_bit = 0;
  //int L1_ps = -1;
  //int HLT_ps = 0;

  HLTConfigProvider const&  hltConfig = hltPrescaleProvider_.hltConfigProvider();

  //cout<<"Currently analyzing trigger "<<triggerName<<endl;

  //Check the current configuration to see how many total triggers there are
  const unsigned int n(hltConfig.size());
  //Get the trigger index for the current trigger
  const unsigned int triggerIndex(hltConfig.triggerIndex(triggerName));
  //check that the trigger in the event and in the configuration agree
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));
  // abort on invalid trigger name
  if (triggerIndex>=n) {
    cout << "Trigger path"<< triggerName << " not found!" << endl;
    throw exception();
  }
  
   acc_bit = triggerResultsHandle_->accept(triggerIndex);
  
  return acc_bit;
}


//--------------------------analyzeSimplePrescales() ------------------------

void SimpleTriggerAnalyzer::initPattern(const edm::TriggerResults & result,
                        const edm::EventSetup& iSetup,
                        const edm::TriggerNames & triggerNames)
//--------------------------------------------------------------------------
{
    unsigned int n;

    // clean up old data
    HLTPathsByName_.clear();

    if (HLTPatterns_.empty()) {
        // for empty input vector, default to all HLT trigger paths
        n = result.size();
        HLTPathsByName_.resize(n);
        for (unsigned int i = 0; i < n; ++i) {
            HLTPathsByName_[i] = triggerNames.triggerName(i);
        }
    } else {
        // otherwise, expand wildcards in trigger names...
        BOOST_FOREACH(const std::string & pattern, HLTPatterns_) {
            if (edm::is_glob(pattern)) {
                // found a glob pattern, expand it
                std::vector< std::vector<std::string>::const_iterator > matches = edm::regexMatch(triggerNames.triggerNames(), pattern);
                if (matches.empty()) {
                    // pattern does not match any trigger paths
                    std::cout<<"No patterns found.  Please check quality file... "<<std::endl;
                    exit(0);
                } else {
                    // store the matching patterns
                    BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches)
                        HLTPathsByName_.push_back(*match);
                }
            } else {
                // found a trigger name, just copy it
                HLTPathsByName_.push_back(pattern);
            }
        }

    }
}



// ------------ method called once each job just before starting event loop  ------------
//void
//SimpleTriggerAnalyzer::beginJob()
//{
//}

// ------------ method called once each job just after ending the event loop  ------------
//void
//SimpleTriggerAnalyzer::endJob()
//{
//}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//void
//SimpleTriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
//  edm::ParameterSetDescription desc;
//  desc.setUnknown();
//  descriptions.addDefault(desc);

//}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleTriggerAnalyzer);
