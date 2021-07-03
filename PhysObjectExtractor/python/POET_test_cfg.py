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

#select -1 for all events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/000A9D3F-CE4C-E311-84F8-001E673969D2.root'
    )
)

process.myevents = cms.EDAnalyzer('EventAnalyzer')	                             
process.mygenparticle= cms.EDAnalyzer('GenParticleAnalyzer',
			#collect particles with specific pdgid:status
			#if 0:0, collect them all	
			input_particle = cms.vstring("0:0")
			)

process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("myoutput.root"))


process.p = cms.Path(process.myevents+process.mygenparticle)