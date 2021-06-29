import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#Work with data (if False, assumed MC simulations)
isData = False
#Get jet corrections using PAT (Physics Analysis Tool) infrastructure
doPat = True

process = cms.Process("POET")

#Configure the framework messaging system
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("POET")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

#Select the maximum number of events to process (if -1, run over all events)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

#Load needed configuration
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

#Define the source files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'root://eospublic.cern.ch//eos/opendata/cms/Run2012B/DoubleMuParked/AOD/22Jan2013-v1/10000/1EC938EF-ABEC-E211-94E0-90E6BA442F24.root'
      #'file:/playground/002F62E1-B53D-E311-A49F-003048F1B950.root'
       'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/000A9D3F-CE4C-E311-84F8-001E673969D2.root'
    )
)

#Alternatively, to run on larger scale, one could use index files as obtained from the Cern Open Data Portal
#files = FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_10000_file_index.txt")
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_20000_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_20001_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_20002_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_210000_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_30000_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_310000_file_index.txt"))
#process.source = cms.Source(
#    "PoolSource", fileNames=cms.untracked.vstring(*files))



#These two lines are needed if you require access to the conditions database. E.g., to get jet energy corrections, trigger prescales, etc.
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
#Uncomment and arrange a line like this if you are getting access to the conditions database through CVMFS snapshot files (requires installing CVMFS client)
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')
#The global tag must correspond to the needed epoch (comment out if no conditions needed)
if isData: process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'
else: process.GlobalTag.globaltag = "START53_V27::All"

if isData:
	# Apply JSON file with lumi mask for data quality purposes (needs to be done after the process.source definition)
	goodJSON = "data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt"
	myLumis = LumiList.LumiList(filename=goodJSON).getCMSSWString().split(",")
	process.source.lumisToProcess = CfgTypes.untracked(
	    	CfgTypes.VLuminosityBlockRange())
	process.source.lumisToProcess.extend(myLumis)


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

 if isData:
   runOnData(process, ['All'], "", None, [])

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

process.mypvertex = cms.EDAnalyzer('VertexAnalyzer')
process.mytracks= cms.EDAnalyzer('TrackAnalyzer')
process.mygenparticle= cms.EDAnalyzer('GenParticleAnalyzer',
			#collect particles with specific pdgid:status
			#if 0:0, collect them all	
			input_particle = cms.vstring("1:11","1:13","1:22","2:15")
			)

process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("myoutput.root"))

if doPat:
	process.p = cms.Path(process.patDefaultSequence+process.myevents+process.myelectrons+process.mymuons+process.myphotons+process.myjets+process.mymets+process.mytaus+process.mytrigEvent+process.mypvertex+process.mytracks+process.mygenparticle)
else: process.p = cms.Path(process.myevents+process.myelectrons+process.mymuons+process.myphotons+process.myjets+process.mymets+process.mytaus+process.mytrigEvent+process.mypvertex+process.mytracks+process.mygenparticle)
