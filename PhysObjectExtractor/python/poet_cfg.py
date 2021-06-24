import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

isData = True

process = cms.Process("POET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("POET")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/SingleElectron/AOD/12Oct2013-v1/10000/1045436C-1240-E311-851B-003048D2BF1C.root'
#	 'file:/playground/002F62E1-B53D-E311-A49F-003048F1B950.root'
  
    )
)

#These two lines are needed if you require access to the conditions database. E.g., to get jet energy corrections, trigger prescales, etc.
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'

#Uncomment this line if you are getting access to the conditions database through CVMFS snapshot files (requires installing CVMFS client)
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')

#Here, you can enter the desired input tag, corresponding to each container, In addition, you can add more containers.
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
if isData:
    process.myjets= cms.EDAnalyzer('JetAnalyzer',
                                   InputCollection = cms.InputTag("ak5PFJets"),
                                   isData = cms.bool(True),
                                   jecL1Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/FT53_V21A_AN6_L1FastJet_AK5PF.txt'),
                                   jecL2Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/FT53_V21A_AN6_L2Relative_AK5PF.txt'),
                                   jecL3Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/FT53_V21A_AN6_L3Absolute_AK5PF.txt'),
                                   jecResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/FT53_V21A_AN6_L2L3Residual_AK5PF.txt'),
                                   jecUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/FT53_V21A_AN6_Uncertainty_AK5PF.txt'),
                               )
else:
    process.myjets= cms.EDAnalyzer('JetAnalyzer',
                                   InputCollection = cms.InputTag("ak5PFJets"),
                                   isData = cms.bool(False),
                                   jecL1Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/START53_V27_L1FastJet_AK5PF.txt'),
                                   jecL2Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/START53_V27_L2Relative_AK5PF.txt'),
                                   jecL3Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/START53_V27_L3Absolute_AK5PF.txt'),
                                   jecResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/START53_V27_L2L3Residual_AK5PF.txt'),
                                   jecUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/START53_V27_Uncertainty_AK5PF.txt'),
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

process.mypvertex = cms.EDAnalyzer('VertexAnalyzer')
process.mytracks= cms.EDAnalyzer('TrackAnalyzer')
process.mygenparticle= cms.EDAnalyzer('GenParticleAnalyzer',input_particle = cms.vstring("electron","muon","photon","tau"))

process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("myoutput.root"))

process.p = cms.Path(process.myevents+process.myelectrons+process.mymuons+process.myphotons+process.myjets+process.mymets+process.mytaus+process.mytrigEvent+process.mypvertex+process.mytracks+process.mygenparticle)

