import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
from PhysicsTools.PatAlgos.patTemplate_cfg import *

isData = False
doPat = True

process = cms.Process("POET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("POET")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/SingleElectron/AOD/12Oct2013-v1/10000/1045436C-1240-E311-851B-003048D2BF1C.root'
      #'file:/playground/002F62E1-B53D-E311-A49F-003048F1B950.root'
       'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/000A9D3F-CE4C-E311-84F8-001E673969D2.root'
    )
)

#These two lines are needed if you require access to the conditions database. E.g., to get jet energy corrections, trigger prescales, etc.
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
if isData: process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'
else: process.GlobalTag.globaltag = "START53_V27::All"

#Uncomment this line if you are getting access to the conditions database through CVMFS snapshot files (requires installing CVMFS client)
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')

#More information about InputCollections at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable
process.myevents = cms.EDAnalyzer('EventAnalyzer')	                             
process.myelectrons = cms.EDAnalyzer('ElectronAnalyzer',
				     InputCollection = cms.InputTag("gsfElectrons")
				    )
process.mymuons = cms.EDAnalyzer('MuonAnalyzer',
				 InputCollection = cms.InputTag("muons")
				 )
process.myphotons = cms.EDAnalyzer('PhotonAnalyzer',
                                   InputCollection = cms.InputTag("photons")
                             )
#Path Strings: These correspond to the Global Tag. Run jec_cfg.py first to get .txt files
JecString = 'START53_V27_'
if isData: JecString = 'FT53_V21A_AN6_'

if doPat:
 # Load PAT config
 process.load("PhysicsTools.PatAlgos.patSequences_cff")
 process.load('Configuration.StandardSequences.Reconstruction_cff')
 process.load('RecoJets.Configuration.RecoPFJets_cff')
 process.load('RecoJets.Configuration.RecoJets_cff')
 process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
 process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

 from PhysicsTools.PatAlgos.tools.pfTools import *
 from PhysicsTools.PatAlgos.tools.coreTools import *
 from PhysicsTools.PatAlgos.tools.metTools import *	
 from PhysicsTools.PatAlgos.tools.jetTools import *
 from PhysicsTools.PatAlgos.tools.coreTools import *
 from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

 # Set up the new jet collection
 process.ak5PFJets.doAreaFastjet = True
 addPfMET(process, 'PF')
 
 addJetCollection(process,cms.InputTag('ak5PFJets'),
 		 'AK5', 'PFCorr',
		 doJTA        = True,
		 doBTagging   = True,
		 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative','L3Absolute'])),
		 doType1MET   = True,
		 doL1Cleaning = True,
		 doL1Counters = False,
		 doJetID      = True,
		 jetIdLabel   = "ak5",
		 )
 process.myjets= cms.EDAnalyzer('PatJetAnalyzer',
				   InputCollection = cms.InputTag("selectedPatJetsAK5PFCorr"),
                                   isData = cms.bool(isData),
                                   jecUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'Uncertainty_AK5PF.txt'), 
                                   jerResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/JetResolutionInputAK5PF.txt')         
                               )
else:
    process.myjets= cms.EDAnalyzer('JetAnalyzer',
                                   InputCollection = cms.InputTag("ak5PFJets"),
                                   isData = cms.bool(isData),
                                   jecL1Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L1FastJet_AK5PF.txt'), 
                                   jecL2Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L2Relative_AK5PF.txt'),     #Don't forget to run jec_cfg.py
                                   jecL3Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L3Absolute_AK5PF.txt'),     #to get these .txt files :)
                                   jecResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L2L3Residual_AK5PF.txt'),
                                   jecUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'Uncertainty_AK5PF.txt'),
                                   jerResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/JetResolutionInputAK5PF.txt')
                               )
process.mymets= cms.EDAnalyzer('MetAnalyzer',
                               InputCollection = cms.InputTag("pfMet")
                              )
process.mytaus = cms.EDAnalyzer('TauAnalyzer',
                                InputCollection = cms.InputTag("hpsPFTauProducer")
                               )
process.mytrigEvent = cms.EDAnalyzer('TriggObjectAnalyzer',
                                     filterName = cms.string("hltSingleJet190Regional"),
                             )

process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("myoutput.root"))

if doPat:
	process.p = cms.Path(process.patDefaultSequence+process.myevents+process.myelectrons+process.mymuons+process.myphotons+process.myjets+process.mymets+process.mytaus+process.mytrigEvent)
else: process.p = cms.Path(process.myevents+process.myelectrons+process.mymuons+process.myphotons+process.myjets+process.mymets+process.mytaus+process.mytrigEvent)
