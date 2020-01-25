import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/SingleElectron/AOD/12Oct2013-v1/10000/1045436C-1240-E311-851B-003048D2BF1C.root'
    )
)

#needed to get the actual prescale values used from the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')
process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'

#Here, you can enter the desired input tag, corresponding to each container, In addition, you can add more containers.
process.demo = cms.EDAnalyzer('PhysicsObjectsInfo'
                              ,ElectronInputCollection = cms.InputTag("gsfElectrons")
                              ,JetInputCollection = cms.InputTag("ak5PFJets")
                              ,MetInputCollection = cms.InputTag("pfMet")
                              ,MuonInputCollection = cms.InputTag("muons")
                              ,PhotonInputCollection = cms.InputTag("photons")
                              ,TauInputCollection = cms.InputTag("hpsPFTauProducer")
                              ,filterName = cms.string("hltSingleJet190Regional")
                             )
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable

process.p = cms.Path(process.demo)
