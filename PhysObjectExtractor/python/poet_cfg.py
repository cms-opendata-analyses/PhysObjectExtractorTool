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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring( 
            # 'root://eospublic.cern.ch//eos/opendata/cms/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/00039A2E-D7A7-E511-98EE-3417EBE64696.root'   
            'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0005EA25-8CB8-E511-A910-00266CF85DA0.root'   
                )
                            )

#globaltag for 2015 collision data
#process.load("Configuration.StandardSequences.Services_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/76X_dataRun2_16Dec2015_v0.db')
#process.GlobalTag.globaltag = '76X_dataRun2_16Dec2015_v0'

process.myelectrons = cms.EDAnalyzer('ElectronAnalyzer', electrons = cms.InputTag("slimmedElectrons"), 
                               vertices=cms.InputTag("offlineSlimmedPrimaryVertices"))
                              
process.mymuons = cms.EDAnalyzer('MuonAnalyzer', muons = cms.InputTag("slimmedMuons"), 
                               vertices=cms.InputTag("offlineSlimmedPrimaryVertices"))

process.mytaus = cms.EDAnalyzer('TauAnalyzer', taus=cms.InputTag("slimmedTaus"))

process.myphotons = cms.EDAnalyzer('PhotonAnalyzer', photons=cms.InputTag("slimmedPhotons"))

process.mymets = cms.EDAnalyzer('MetAnalyzer', mets=cms.InputTag("slimmedMETs"))

process.mytriggers = cms.EDAnalyzer('TriggObjectAnalyzer', objects = cms.InputTag("selectedPatTrigger"))

process.mypvertex = cms.EDAnalyzer('VertexAnalyzer', vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
beams=cms.InputTag("offlineBeamSpot"))

process.mygenparticle = cms.EDAnalyzer('GenParticleAnalyzer', pruned=cms.InputTag("prunedGenParticles"),
                                       #---- Collect particles with specific "pdgid:status"
				                       #---- if 0:0, collect them all	
				                       input_particle = cms.vstring("1:11","1:13","1:22","2:15"))

process.myjets = cms.EDAnalyzer('JetAnalyzer', jets = cms.InputTag("slimmedJets"))

process.TFileService = cms.Service("TFileService", fileName=cms.string("myoutput.root"))

process.p = cms.Path(process.myelectrons+process.mymuons+process.mytaus+process.myphotons+process.mymets+process.mytriggers+process.mypvertex+process.mygenparticle+process.myjets)

