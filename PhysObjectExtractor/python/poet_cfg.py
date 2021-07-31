import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("POET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("POET")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'root://cms-xrd-global.cern.ch//store/data/Run2018A/DoubleMuon/MINIAOD/17Sep2018-v2/00000/0068887D-518C-2B47-8C29-E699E295A2CC.root'
                )
                            )

process.myelectrons = cms.EDAnalyzer('ElectronAnalyzer',electrons = cms.InputTag("slimmedElectrons"), 
                               vertices=cms.InputTag("offlineSlimmedPrimaryVertices"))
                              

process.TFileService = cms.Service("TFileService", fileName=cms.string("myoutput.root"))

process.p = cms.Path(process.myelectrons)
