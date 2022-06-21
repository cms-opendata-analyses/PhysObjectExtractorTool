import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("POET")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

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

#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)




#
# Configure an example module for user analysis with electrons
#

process.ntupler = cms.EDAnalyzer('ElectronAnalyzer',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #
                                 # ... none ...
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 #
                                 # ID decisions (common to all formats)
                                 eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
                                 eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
                                 #
                                 # ValueMaps with MVA results
                                 #
                                 mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
                                 mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories")
                                )

                              
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

process.p = cms.Path(process.egmGsfElectronIDSequence * process.ntupler+process.mymuons+process.mytaus+process.myphotons+process.mymets+process.mytriggers+process.mypvertex+process.mygenparticle+process.myjets)

