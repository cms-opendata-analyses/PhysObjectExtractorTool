import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import sys

if len(sys.argv) > 2:
    isData = eval(sys.argv[2])
else:
    isData = False
isMC = True
if isData: isMC = False

process = cms.Process("POET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("POET")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(         
        #'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0005EA25-8CB8-E511-A910-00266CF85DA0.root'   
        'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000/0007DBD0-2ED2-E511-AD0D-20CF3019DEF5.root'
        )
)
if isData:
    process.source.fileNames = cms.untracked.vstring(
        'root://eospublic.cern.ch//eos/opendata/cms/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/00039A2E-D7A7-E511-98EE-3417EBE64696.root'
        )
    #---- Apply the data quality JSON file filter. This example is for 2015 data
    #---- It needs to be done after the process.source definition
    #---- Make sure the location of the file agrees with your setup
    goodJSON = "data/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt"
    myLumis = LumiList.LumiList(filename=goodJSON).getCMSSWString().split(",")
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)


#---- These two lines are needed if you require access to the conditions database. E.g., to get jet energy corrections, trigger prescales, etc.
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#---- Uncomment and arrange a line like this if you are getting access to the conditions database through CVMFS snapshot files (requires installing CVMFS client)
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/76X_dataRun2_16Dec2015_v0.db')
#---- The global tag must correspond to the needed epoch (comment out if no conditions needed)
#if isData: process.GlobalTag.globaltag = '76X_dataRun2_16Dec2015_v0'
#else: process.GlobalTag.globaltag = "76X_mcRun2_asymptotic_RunIIFall15DR76_v1"
#---- If the container has local DB files available, uncomment lines like the ones below
#---- instead of the corresponding lines above
#if isData: process.GlobalTag.connect = cms.string('sqlite_file:/opt/cms-opendata-conddb/FT53_V21A_AN6_FULL_data_stripped.db')
#else:  process.GlobalTag.connect = cms.string('sqlite_file:/opt/cms-opendata-conddb/START53_V27_MC_stripped.db')
#if isData: process.GlobalTag.globaltag = 'FT53_V21A_AN6_FULL::All'
#else: process.GlobalTag.globaltag = "START53_V27::All"


#globaltag for 2015 collision data
#process.load("Configuration.StandardSequences.Services_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/76X_dataRun2_16Dec2015_v0.db')
#process.GlobalTag.globaltag = '76X_dataRun2_16Dec2015_v0'

process.myelectrons = cms.EDAnalyzer('ElectronAnalyzer', electrons = cms.InputTag("slimmedElectrons"), 
                               vertices=cms.InputTag("offlineSlimmedPrimaryVertices"))
                              
process.mymuons = cms.EDAnalyzer('MuonAnalyzer', muons = cms.InputTag("slimmedMuons"), vertices=cms.InputTag("offlineSlimmedPrimaryVertices"))

process.mytaus = cms.EDAnalyzer('TauAnalyzer', taus=cms.InputTag("slimmedTaus"))

process.myphotons = cms.EDAnalyzer('PhotonAnalyzer', photons=cms.InputTag("slimmedPhotons"))

process.mytriggers = cms.EDAnalyzer('TriggObjectAnalyzer', objects = cms.InputTag("selectedPatTrigger"))

process.mypvertex = cms.EDAnalyzer('VertexAnalyzer', vertices=cms.InputTag("offlineSlimmedPrimaryVertices"), beams=cms.InputTag("offlineBeamSpot"))

process.mygenparticle = cms.EDAnalyzer('GenParticleAnalyzer', pruned=cms.InputTag("prunedGenParticles"),
                                       #---- Collect particles with specific "pdgid:status"
                           #---- if 0:0, collect them all 
                           input_particle = cms.vstring("1:11","1:13","1:22","2:15"))

JecString = 'MC'
if isData: JecString = 'DATA'
process.myjets = cms.EDAnalyzer('JetAnalyzer', 
				jets = cms.InputTag("slimmedJets"),
				isData = cms.bool(isData),
                                jecL1Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_'+JecString+'_L1FastJet_AK4PFchs.txt'), 
                                jecL2Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_'+JecString+'_L2Relative_AK4PFchs.txt'), 
                                jecL3Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_'+JecString+'_L3Absolute_AK4PFchs.txt'), 
                                jecResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt'), 
				jetJECUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt'), 
                                jerResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt'),
                                jerSFName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_SF_AK4PFchs.txt'),
				)
process.myfatjets = cms.EDAnalyzer('FatjetAnalyzer', 
				fatjets = cms.InputTag("slimmedJetsAK8"),
				isData = cms.bool(isData),
                                jecL1Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_'+JecString+'_L1FastJet_AK8PFchs.txt'), 
                                jecL2Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_'+JecString+'_L2Relative_AK8PFchs.txt'), 
                                jecL3Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_'+JecString+'_L3Absolute_AK8PFchs.txt'), 
                                jecResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_DATA_L2L3Residual_AK8PFchs.txt'), 
				jetJECUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_Uncertainty_AK8PFchs.txt'), 
                                jerResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_PtResolution_AK8PFchs.txt'),
                                jerSFName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_SF_AK4PFchs.txt'), # AK8 == AK4
				)

process.mymets = cms.EDAnalyzer('MetAnalyzer',mets=cms.InputTag("slimmedMETs"))

process.TFileService = cms.Service("TFileService", fileName=cms.string("myoutput.root"))

process.p = cms.Path(process.myelectrons+process.mymuons+process.mytaus+process.myphotons+process.mymets+process.mytriggers+process.mypvertex+process.mygenparticle+process.myjets+process.myfatjets)
