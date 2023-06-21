// -*- C++ -*-
//
// Package:    Trigger/TriggerAnalyzer
// Class:      TriggerAnalyzer
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



class TriggerAnalyzer : public edm::stream::EDAnalyzer< >  {
   public:
      explicit TriggerAnalyzer(const edm::ParameterSet&);
      ~TriggerAnalyzer();
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void analyzeSimplePrescales(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);
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
     std::map<std::string, int> trigmap;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
// constructors and destructor
// Notice that here, using the parameter set tool,
// you need to point the code to
// the right branch (in the EDM root files) where the trigger HLTPattrmation
// is stored.
// Also, at configuration time
// you will need to point to the appropiate triggers
// you want to look at. Alternatively (is also shown below), you
// could select the trigger names dynamically; for example getting them
// from the HLTConfigProvider.
// To start out, you need to define a processName, which is the name of
// the CMSSW computing process that originally wrote the products in the root
// file. Originally, this is always "HLT", by default.
// In triggerName, you can
// use wildcards, which will be described later.
// As for the InputTags, these shall match the name of the ROOT branches
// where the HLTPattrmation is stored.  This was essentially fixed and will
// most likely be the same always.
//This should match your configuration python file
//
TriggerAnalyzer::TriggerAnalyzer(const edm::ParameterSet& ps):
processName_(ps.getParameter<std::string>("processName")),
triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
triggerResultsToken_(consumes<edm::TriggerResults>(triggerResultsTag_)),
triggerNamesID_(),
hltPrescaleProvider_(ps, consumesCollector(), *this),
HLTPatterns_(ps.getParameter<std::vector<std::string> >("triggerPatterns")),
HLTPathsByName_()
{
   //now do what ever initialization is needed
   using namespace std;
   using namespace edm;

   //for storing HLTPatt
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
   mtree->Branch("triggermap", &trigmap);
   //second  stores the multiplication acceptbit*L1ps*HLTps
   //so, if negative, it means that L1ps couldn't be found.
   //look below in the code to understand the specifics
   mtree->GetBranch("triggermap")->SetTitle("first:name of trigger, second: acceptbit*L1ps*HLTps");
}


TriggerAnalyzer::~TriggerAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
// ------------ method called when ending the processing of a run  ------------
void TriggerAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a run  ------------
void TriggerAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
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
TriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   // Get event products:
   // In the following, the code is trying to access the HLTPattrmation
   // from the ROOT files and point the containers (that we created),
   // namely triggerResultsHandle_,
   // to the correct "address", given at configuration time
   // and assigned to triggerResultsTag_
    // After that, a simple sanity check is done.

   iEvent.getByToken(triggerResultsToken_,triggerResultsHandle_);
   if (!triggerResultsHandle_.isValid()) {
       LogVerbatim("TriggerAnalyzer") << "TriggerAnalyzer::analyze: Error in getting TriggerResults product from Event!" << endl;
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

   unsigned int n = HLTPathsByName_.size();

   //clear trigmap
   trigmap.clear();

   //Loop over all triggers in the pattern
   for (unsigned int i=0; i!=n; ++i) {
       analyzeSimplePrescales(iEvent,iSetup,HLTPathsByName_[i]);
   }

   //Printout for testing
   //map<string, int>::iterator it;
   //for (it = trigmap.begin(); it != trigmap.end(); it++)
   //{
   //    cout << "********Trigger Name: "<<it->first<<"  accbit*L1ps*HLTps: "<<it->second<< endl;
  //}


  mtree->Fill();
  return;


}

//---------------------------Actual trigger analysis-------------
void TriggerAnalyzer::analyzeSimplePrescales(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName)
//-----------------------------------------------------------------
{

  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  int acc_bit = 0;
  int L1_ps = -1;
  int HLT_ps = 0;

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


  //  L1 and HLT prescale values via (L1) EventSetup
  // Current (default) prescale set index - to be taken from L1GtUtil via Event.
  // Try to get the L1 and HLT prescale values that were actually used
  // for this event.
  // This example needs the conditions stored in the Global Tag,
  // which is some sort of snapshot of the by-then current detector
  // conditions.  They need to be extracted from the database
  // and for that the "Frontier Conditions" lines need to
  // be added in the python configuration file along with the
  // name for the global tag.
  const std::pair<int,int> prescales(hltPrescaleProvider_.prescaleValues(iEvent,iSetup,triggerName));
  L1_ps = prescales.first;
  HLT_ps = prescales.second;

  // Find out if the trigger was active, accepted, or in error.
  // We could also find out whether the trigger was active (wasrun),
  // if it accepted the event (accept) or if it gave an error (error).
  // Results from TriggerResults product
  int wr = triggerResultsHandle_->wasrun(triggerIndex);
  int acc = triggerResultsHandle_->accept(triggerIndex);
  int err = triggerResultsHandle_->error(triggerIndex);
  //cout<<"wr = "<<wr<<endl;
  //cout<<"acc = "<<acc<<endl;
  //cout<<"err = "<<err<<endl;

  //just interested on whether the trigger was fired
  if (wr == 1 && acc == 1 && err == 0){
     acc_bit = 1; //Accepted State
  }

  //if accepted the prescale will be different from zero but positive
  //if accepted but L1 prescale can't be read, it will be negative
  trigmap.insert(pair<string,int>(triggerName,acc_bit*L1_ps*HLT_ps));

  return;
}

//--------------------------analyzeSimplePrescales() ------------------------

void TriggerAnalyzer::initPattern(const edm::TriggerResults & result,
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
//TriggerAnalyzer::beginJob()
//{
//}

// ------------ method called once each job just after ending the event loop  ------------
//void
//TriggerAnalyzer::endJob()
//{
//}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//void
//TriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
//  edm::ParameterSetDescription desc;
//  desc.setUnknown();
//  descriptions.addDefault(desc);

//}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzer);
