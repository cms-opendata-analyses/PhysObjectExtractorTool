// -*- C++ -*-
//
// Package:    SimpleMuTauFilter
// Class:      SimpleMuTauFilter
// 
/**\class SimpleMuTauFilter SimpleMuTauFilter.cc PhysObjectExtractorTool/SimpleMuTauFilter/src/SimpleMuTauFilter.cc

 Description: [one line class summary]

This is a simple filter example to filter on at least a muon and at least
a tau of certain characteristics, which are mostly configurable.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Sat Jul 17 22:23:23 CEST 2021
// $Id$
//
//


// system include files
#include <memory>
#include <math.h> 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//classes to extract Muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//classes to extract Electron information
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//
// class declaration
//

class Simple4LeptonFilter : public edm::EDFilter {
   public:
      explicit Simple4LeptonFilter(const edm::ParameterSet&);
      ~Simple4LeptonFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      
      virtual float effectiveArea0p3cone(float eta);

      // ----------member data ---------------------------
  edm::InputTag muonInput;
  edm::InputTag electronInput;  

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
Simple4LeptonFilter::Simple4LeptonFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  muonInput = iConfig.getParameter<edm::InputTag>("InputCollectionMuons");
  electronInput = iConfig.getParameter<edm::InputTag>("InputCollectionElectrons");

}


Simple4LeptonFilter::~Simple4LeptonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

float Simple4LeptonFilter::effectiveArea0p3cone(float eta){
  if(fabs(eta) < 1.0) return 0.13;
  else if(fabs(eta) < 1.479) return 0.14;
  else if(fabs(eta) < 2.0) return 0.07;
  else if(fabs(eta) < 2.2) return 0.09;
  else if(fabs(eta) < 2.3) return 0.11;
  else if(fabs(eta) < 2.4) return 0.11;
  else return 0.14;
}

// ------------ method called on each new Event  ------------
bool Simple4LeptonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;
  using namespace std;


 //Filter on at least one good muon
  Handle<reco::MuonCollection> mymuons;
  iEvent.getByLabel(muonInput, mymuons);
  
  Handle<reco::GsfElectronCollection> myelectrons;
  iEvent.getByLabel(electronInput, myelectrons);

  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
  math::XYZPoint pv(vertices->begin()->position());
  
  Handle<double> rhoHandle;
  iEvent.getByLabel(InputTag("fixedGridRhoAll"), rhoHandle);

  bool isGood4Lepton = false;
  
  int nummuons = 0;
  int goodmuonpt = 0;
  int goodmuoneta = 0;
  int goodmuon_pfreliso04all = 0;
  int goodmuon_sip3d = 0;
  int goodmuon_dxy = 0;
  int goodmuon_dz = 0;
  float muon_dxy = -999;
  float muon_dxyErr = -999;
  float muon_dz = -999;
  float muon_dzErr = -999;
  float muon_ip3d = sqrt(muon_dxy*muon_dxy + muon_dz*muon_dz);
  float muon_sip3d = muon_ip3d/sqrt(muon_dxyErr*muon_dxyErr + muon_dzErr*muon_dzErr);
  float muon_pfreliso04all = -999;
  float muon_ch = 0;
  int muon_postch = 0;
  int muon_negch = 0;
  float muon_phi1 = 0;
  float muon_phi2 = 0;
  float muon_eta1 = 0;
  float muon_eta2 = 0;
  float muon_dr = 0;
  float muon_pt1 = 0;
  float muon_pt2 = 0;
  bool goodmuonpair_pt = false;
  int i = 0;
  
  int numElec = 0;
  int goodElecpt = 0;
  int goodEleceta = 0;
  int goodElec_pfreliso03all = 0;
  int goodElec_sip3d = 0;
  int goodElec_dxy = 0;
  int goodElec_dz = 0;
  float Elec_dxy = -999;
  float Elec_dxyErr = -999;
  float Elec_dz = -999;
  float Elec_dzErr = -999;
  float Elec_ip3d = sqrt(Elec_dxy*Elec_dxy + Elec_dz*Elec_dz);
  float Elec_sip3d = Elec_ip3d/sqrt(Elec_dxyErr*Elec_dxyErr + Elec_dzErr*Elec_dzErr);
  float Elec_pfreliso03all = -999;
  float Elec_ch = 0;
  int Elec_postch = 0;
  int Elec_negch = 0;
  float Elec_phi1 = 0;
  float Elec_phi2 = 0;
  float Elec_eta1 = 0;
  float Elec_eta2 = 0;
  float Elec_dr = 0;
  float Elec_pt1 = 0;
  float Elec_pt2 = 0;
  bool goodElecpair_pt = false;
  int j = 0;
  
  bool goodDR = true; //Good DeltaR
  
  if(mymuons.isValid()){
  
    nummuons = mymuons->size();
    
    for(reco::MuonCollection::const_iterator itmuon=mymuons->begin(); itmuon!=mymuons->end(); ++itmuon){
    
       muon_ch = itmuon->charge();
       
       if(nummuons == 2){
         if(i == 0){
           muon_eta1 = itmuon -> eta();
           muon_phi1 = itmuon -> phi();
           muon_pt1 = itmuon -> pt();
           i++;
         }
         
         if(i == 1){
           muon_eta2 = itmuon -> eta();
           muon_phi2 = itmuon -> phi();
           muon_pt2 = itmuon -> pt();
         }
           
       }
       
       if(muon_ch > 0){
         muon_postch++;
       }
       
       if(muon_ch < 0){
         muon_negch++;
       }
      
       if(itmuon->pt() > 5){
         goodmuonpt++;  //Good Muon Transversal Momentum
       }
       
       if(abs(itmuon->eta()) < 2.4){
         goodmuoneta++; //Good Muon Eta
       }
       
       if(itmuon->isPFMuon() && itmuon->isPFIsolationValid()) {
         auto iso04 = itmuon->pfIsolationR04();
         muon_pfreliso04all = (iso04.sumChargedHadronPt + iso04.sumNeutralHadronEt + iso04.sumPhotonEt)/itmuon->pt();
         } 
       
       if(abs(muon_pfreliso04all)<0.40){
         goodmuon_pfreliso04all++; //Good Muon Isolation
         }
       
       auto trk = itmuon->globalTrack();
       
       if(trk.isNonnull()) {
         muon_dxy = trk->dxy(pv);
         muon_dz = trk->dz(pv);
         muon_dxyErr = trk->d0Error();
         muon_dzErr = trk->dzError();
         muon_ip3d = sqrt(muon_dxy*muon_dxy + muon_dz*muon_dz);
         muon_sip3d = muon_ip3d/sqrt(muon_dxyErr*muon_dxyErr + muon_dzErr*muon_dzErr);
         } 
       
       if(muon_sip3d < 4){
         goodmuon_sip3d++;
         }
       
       if(abs(muon_dxy)<0.5){
         goodmuon_dxy++;
         }
         
       if(abs(muon_dz)<1.0){
         goodmuon_dz++;
         }
         
       }  
    
  }
  
  if(nummuons == 2){

      muon_dr = sqrt(pow(muon_eta1 - muon_eta2, 2) + pow(muon_phi1 - muon_phi2, 2));
    
      if(muon_pt1 > 20 && muon_pt2 >10){
        goodmuonpair_pt = true;
      }
      
      if(muon_pt2 > 20 && muon_pt1 >10){
        goodmuonpair_pt = true;
      } 
    
    }
  
  if(myelectrons.isValid()){
  
    numElec = myelectrons->size();
    
    for(reco::GsfElectronCollection::const_iterator itElec=myelectrons->begin(); itElec!=myelectrons->end(); ++itElec){
    
       Elec_ch = itElec->charge();
       
       if(numElec == 2){
         if(j == 0){
           Elec_eta1 = itElec -> eta();
           Elec_phi1 = itElec -> phi();
           Elec_pt1 = itElec -> pt();
           j++;
         }
         
         if(j == 1){
           Elec_eta2 = itElec -> eta();
           Elec_phi2 = itElec -> phi();
           Elec_pt2 = itElec -> pt();
         }
           
       }
       
       if(Elec_ch > 0){
         Elec_postch++;
       }
       
       if(Elec_ch < 0){
         Elec_negch++;
       }
       
       if(itElec->pt() > 7){
         goodElecpt++;  //Good Electron Transversal Momentum
       }
       
       if(abs(itElec->eta()) < 2.5){
         goodEleceta++; //Good Electron Eta
       }
       
       if(itElec->passingPflowPreselection()){
         double rho = 0;
         if(rhoHandle.isValid()) rho = *(rhoHandle.product());
         double Aeff = effectiveArea0p3cone(itElec->eta());
         auto iso03 = itElec->pfIsolationVariables();
         Elec_pfreliso03all = (iso03.chargedHadronIso + std::max(0.0,iso03.neutralHadronIso + iso03.photonIso - rho*Aeff))/itElec->pt();
         } 
           
       if(abs(Elec_pfreliso03all)<0.40){
         goodElec_pfreliso03all++; //Good Electron Isolation
         }
       
       auto trk1 = itElec->gsfTrack();
       
       if(trk1.isNonnull()) {
         Elec_dxy = trk1->dxy(pv);
         Elec_dz = trk1->dz(pv);
         Elec_dxyErr = trk1->d0Error();
         Elec_dzErr = trk1->dzError();
         Elec_ip3d = sqrt(Elec_dxy*Elec_dxy + Elec_dz*Elec_dz);
         Elec_sip3d = Elec_ip3d/sqrt(Elec_dxyErr*Elec_dxyErr + Elec_dzErr*Elec_dzErr);
         }
       
       if(Elec_sip3d < 4){
         goodElec_sip3d++;
         }
       
       if(abs(Elec_dxy)<0.5){
         goodElec_dxy++;
         }
         
       if(abs(Elec_dz)<1.0){
         goodElec_dz++;
         }
         
       }  
    
  }
  
  if(numElec == 2){
    
    Elec_dr = sqrt(pow(Elec_eta1 - Elec_eta2, 2) + pow(Elec_phi1 - Elec_phi2, 2));
         
    if(Elec_pt1 > 20 && Elec_pt2 >10){
      goodElecpair_pt = true;
    }
         
    if(Elec_pt2 > 20 && Elec_pt1 >10){
       goodElecpair_pt = true;
    } 
       
  }
        
  
  if(muon_dr < 0.02 || Elec_dr < 0.02){ 
    goodDR = false;
    }

  //4 Muons in the Event
  if(nummuons >= 4){ //Good Number of Muons      
      if(goodmuon_pfreliso04all == nummuons){ //All Muons Must Have a Good Isolation
        if(goodmuonpt == nummuons && goodmuoneta == nummuons){ //Good Muon Kinematics
          if(goodmuon_sip3d == nummuons && goodmuon_dxy == nummuons && goodmuon_dz == nummuons){ //Muon Track close to primary vertex with small uncertainty   
            if(nummuons == 4 && muon_postch == 2 && muon_negch == 2){ //At least Two positive and two negative Muons
              return isGood4Lepton = true;
            }
          }
        }
      }    
    }
    
  //4 Electrons in the Event
  if(numElec >= 4){ //Good Number of Electrons     
      if(goodElec_pfreliso03all == numElec){ //All Electrons Must Have a Good Isolation
        if(goodElecpt == numElec && goodEleceta == numElec){ //Good Electron Kinematics
          if(goodElec_sip3d == numElec && goodElec_dxy == numElec && goodElec_dz == numElec){ //Electron Track close to primary vertex with small uncertainty           
            if(numElec == 4 && Elec_postch == 2 && Elec_negch == 2){ //At least Two positive and two negative Electrons
              return isGood4Lepton = true;
            }
          }
        }
      }    
    }
    
  //2 Electrons and 2 Muons in the Event
  if(numElec >= 2 && nummuons >= 2){ //Good Number of Leptons   
    if(goodEleceta == numElec && goodmuoneta == nummuons){ //Eta Cuts
      if(goodElec_pfreliso03all == numElec && goodmuon_pfreliso04all == nummuons){ //All Leptons Must Have a Good Isolation
        if(goodElec_sip3d == numElec && goodElec_dxy == numElec && goodElec_dz == numElec){ //Electron Track close to primary vertex with small uncertainty 
          if(goodmuon_sip3d == nummuons && goodmuon_dxy == nummuons && goodmuon_dz == nummuons){ //Muon Track close to primary vertex with small uncertainty  
            if(Elec_postch == 1 && Elec_negch == 1 && muon_postch == 1 && muon_negch == 1 ){ //Two opposite charged electron and muon pairs
              if(goodElecpair_pt == true && goodmuonpair_pt == true){ //Pt cuts
                if(goodDR == true){ //DeltaR Cuts
                  return isGood4Lepton = true;
                  }
                }
              }
            }          
          }
        }
      }       
    }
  
  using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   return isGood4Lepton;
}

// ------------ method called once each job just before starting event loop  ------------
void 
Simple4LeptonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Simple4LeptonFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
Simple4LeptonFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
Simple4LeptonFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
Simple4LeptonFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
Simple4LeptonFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Simple4LeptonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(Simple4LeptonFilter);
