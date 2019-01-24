import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(25) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/SingleElectron/AOD/12Oct2013-v1/10000/1045436C-1240-E311-851B-003048D2BF1C.root'
    )
)

process.demo = cms.EDAnalyzer('ElectronAnalyzer'
                              ,
                              InputCollection = cms.InputTag("gsfElectrons")
)
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable
process.TFileService = cms.Service("TFileService",
              fileName = cms.string('histo.root')
                                   )


process.p = cms.Path(process.demo)
