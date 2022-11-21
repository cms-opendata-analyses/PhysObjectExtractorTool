import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

import os 
import sys

relBase = os.environ['CMSSW_BASE']

#---- sys.argv takes the parameters given as input cmsRun PhysObjectExtractor/python/poet_cfg.py <isData (default=False)> <doPat (default=False)> 
#----  e.g: cmsRun PhysObjectExtractor/python/poet_cfg.py True True
#---- NB the first two parameters are always "cmsRun" and the config file name
#---- Work with data (if False, assumed MC simulations)
#---- This needs to be in agreement with the input files/datasets below.
if len(sys.argv) > 2:
    isData = eval(sys.argv[2])
else:
    isData = False

#---- Flag for using the Physics Analysis Toolkit for jets and MET
if len(sys.argv) > 3:
    doPat = eval(sys.argv[3])
else:
    doPat = False


process = cms.Process("POET")

#---- Configure the framework messaging system
#---- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("POET")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

#---- Select the maximum number of events to process (if -1, run over all events)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

#---- Load needed configuration
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

#---- Define the test source files to be read using the xrootd protocol (root://), or local files (file:)
#---- Several files can be comma-separated
#---- A local file, for testing, can be downloaded using, e.g., the cern open data client (https://cernopendata-client.readthedocs.io/en/latest/):
#---- python cernopendata-client download-files --recid 6004 --filter-range 1-1
#---- For running over larger number of files, comment out this section and use/uncomment the FileUtils infrastructure below
if isData: 
	sourceFile='root://eospublic.cern.ch//eos/opendata/cms/Run2011A/DoubleMu/AOD/12Oct2013-v1/10000/000D143E-9535-E311-B88B-002618943934.root'
else: 
	sourceFile='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/ZZTo2e2mu_mll4_7TeV-powheg-pythia6/AODSIM/PU_S13_START53_LV6-v1/00000/0A571FF3-6392-E411-AB3D-0025904B12FC.root'
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/playground/1EC938EF-ABEC-E211-94E0-90E6BA442F24.root'
	sourceFile
    )
)

#---- Alternatively, to run on larger scale, one could use index files as obtained from the Cern Open Data Portal
#---- and pass them into the PoolSource.  The example is for 2012 data
#files = FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_10000_file_index.txt")
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_20000_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_20001_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_20002_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_210000_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_30000_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_310000_file_index.txt"))
#process.source = cms.Source(
#    "PoolSource", fileNames=cms.untracked.vstring(*files))



#---- These two lines are needed if you require access to the conditions database. E.g., to get jet energy corrections, trigger prescales, etc.
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
#---- Uncomment the two followingh lines if you are getting access to the conditions database through CVMFS snapshot files (requires installing CVMFS client) 
#if isData: process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb/FT_53_LV5_AN1.db')
#else:  process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb/START53_LV6A1.db')
#---- If the container has local DB files available, uncomment the two lines like the ones below
#---- instead of the corresponding lines above
#if isData: process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_data_stripped.db')
#else:  process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb/START53_LV6A1_MC_stripped.db')
#---- The global tag must correspond to the needed epoch (comment out if no conditions needed)
if isData: process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
else: process.GlobalTag.globaltag = "START53_LV6A1::All"


#---- Apply the data quality JSON file filter. This example is for 2011 data
#---- It needs to be done after the process.source definition
#---- Make sure the location of the file agrees with your setup
if isData:
	goodJSON = "data/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt"
	myLumis = LumiList.LumiList(filename=goodJSON).getCMSSWString().split(",")
	process.source.lumisToProcess = CfgTypes.untracked(
	    	CfgTypes.VLuminosityBlockRange())
	process.source.lumisToProcess.extend(myLumis)


#---- Configure the PhysObjectExtractor modules!

#---- More information about InputCollections at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable
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

#---- Jet correction paths -- these correspond to the Global Tag. **Run jec_cfg.py first to get .txt files!!**
JecString = 'START53_LV6A1_'
if isData: JecString = 'FT_53_LV5_AN1_'

#---- Jets are simpler to work with in "Physics Analysis Toolkit" format. See more at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPAT
if doPat:
	#---- Load PAT configs and build some light sequences to process jets and MET
	process.load('PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff')
	process.load('PhysicsTools.PatAlgos.producersLayer1.metProducer_cff')
	process.load('PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi')
	process.patCandidates = cms.Sequence(process.makePatJets+process.makePatMETs)
	process.selectedPatCandidates = cms.Sequence(process.selectedPatJets)
	process.patDefaultSequence = cms.Sequence(process.patCandidates * process.selectedPatCandidates)
	process.load('RecoJets.Configuration.RecoPFJets_cff')
	from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection, runBTagging
	from PhysicsTools.PatAlgos.tools.coreTools import runOnData

	#---- Choose which jet correction levels to apply
	jetcorrlabels = ['L1FastJet','L2Relative','L3Absolute']
	if isData:
		#---- For data we need to remove generator-level matching processes
		runOnData(process, ['Jets','METs'], "", None, [])
		jetcorrlabels.append('L2L3Residual')

	#---- Configure the addJetCollection tool
	#---- This process will make corrected jets with b-tagging included, and will make Type1-corrected MET
	process.ak5PFJets.doAreaFastjet = True
	addJetCollection(process,cms.InputTag('ak5PFJets'),
			 'AK5', 'PFCorr',
			 doJTA        = True,
			 doBTagging   = True, 
			 jetCorrLabel = ('AK5PF', cms.vstring(jetcorrlabels)),
			 doType1MET   = True,
			 doL1Cleaning = False,
			 doL1Counters = False,
			 doJetID      = True,
			 jetIdLabel   = "ak5",
			 ) 
 
	#---- Configure the POET jet analyzer
	#---- Don't forget to run jec_cfg.py to get these .txt files!
	process.myjets= cms.EDAnalyzer('PatJetAnalyzer',
				       InputCollection = cms.InputTag("selectedPatJetsAK5PFCorr"),
				       isData = cms.bool(isData),
				       jecUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'Uncertainty_AK5PF.txt'), 
				       jerResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/JetResolutionInputAK5PF.txt')         
				       )
else:
	if not isData:
		#---- Get non-PAT access to the jet flavour information
		from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
		process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone()
		from PhysicsTools.JetMCAlgos.AK5PFJetsMCFlavourInfos_cfi import ak5JetFlavourInfos
		process.jetFlavourInfosAK5PFJets = ak5JetFlavourInfos.clone()

	#---- Configure the POET jet analyzer
	#---- Don't forget to run jec_cfg.py to get these .txt files!
	process.myjets= cms.EDAnalyzer('JetAnalyzer',
				       InputCollection = cms.InputTag("ak5PFJets"),
				       isData = cms.bool(isData),
				       jecL1Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L1FastJet_AK5PF.txt'), 
				       jecL2Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L2Relative_AK5PF.txt'),
				       jecL3Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L3Absolute_AK5PF.txt'),
				       jecResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L2L3Residual_AK5PF.txt'),
				       jecUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'Uncertainty_AK5PF.txt'),
				       jerResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/JetResolutionInputAK5PF.txt')
				       )

process.mymets= cms.EDAnalyzer('MetAnalyzer',
                               InputCollection = cms.InputTag("pfMet"),
			       doPat = cms.bool(doPat),
			       )
if doPat: process.mymets.InputCollectionPat = cms.InputTag("pfType1CorrectedMet")

process.mytaus = cms.EDAnalyzer('TauAnalyzer',
                                InputCollection = cms.InputTag("hpsPFTauProducer")
				)

process.mytrigEvent = cms.EDAnalyzer('TriggObjectAnalyzer',
                                     filterName = cms.string("hltL2DoubleMu23NoVertexL2PreFiltered"),
                             )

process.mypvertex = cms.EDAnalyzer('VertexAnalyzer')

process.mytracks= cms.EDAnalyzer('TrackAnalyzer')

process.mygenparticle= cms.EDAnalyzer('GenParticleAnalyzer',
				      #---- Collect particles with specific "pdgid:status"
                                      #---- Check PDG ID in the PDG.
				      #---- if 0:0, collect them all	
				      input_particle = cms.vstring("1:11","1:13","1:22","2:15")
				      )

process.mytriggers = cms.EDAnalyzer('TriggerAnalyzer',
                              processName = cms.string("HLT"),
                              #---- These are example triggers for 2011 DoubleMu dataset
                              #---- Wildcards * and ? are accepted (with usual meanings)
                               #---- If left empty, all triggers will run              
#                              triggerPatterns = cms.vstring("HLT_L2DoubleMu23_NoVertex_v*","HLT_Mu13_Mu8_v*", "HLT_DoubleMu45_v*", "HLT_Mu8_Jet40_v*", "HLT_TripleMu5_v*"), 
                              triggerPatterns = cms.vstring("HLT_L2DoubleMu23_NoVertex_v*","HLT_Mu13_Mu8_v*"),
                              triggerResults = cms.InputTag("TriggerResults","","HLT"),
                              triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")                             
                              )


#---- Example of a CMSSW filter that can be used to cut on a given set of triggers
#---- This filter, however, does know about prescales
#---- A previous trigger study would be needed to cut hard on a given trigger or set of triggers
#---- The filter can be added to the path below if needed but is not applied by default
#process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
#process.hltHighLevel.HLTPaths = cms.vstring('HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*')




#---- Example of a very basic home-made filter to select only events of interest
#---- The filter can be added to the running path below if needed but is not applied by default
#process.mutaufilter = cms.EDFilter('SimpleMuTauFilter',
#                                  InputCollectionMuons = cms.InputTag("muons"),
#                                   InputCollectionTaus = cms.InputTag("hpsPFTauProducer"),
#                                   mu_minpt = cms.double(17),
#                                   mu_etacut = cms.double(2.1),
#                                   tau_minpt = cms.double(20),
#                                   tau_etacut = cms.double(2.3)
#                                   )




#---- Configure the output ROOT filename
process.TFileService = cms.Service(
	"TFileService", fileName=cms.string("myoutput.root"))

#---- Finally run everything! 
#---- Separation by * implies that processing order is important.
#---- separation by + implies that any order will work
#---- One can put in or take out the needed processes
if doPat:
	process.p = cms.Path(process.patDefaultSequence+process.myevents+process.myelectrons+process.mymuons+process.myphotons+process.myjets+process.mymets+process.mytaus+process.mytrigEvent+process.mypvertex+process.mytracks+process.mygenparticle+process.mytriggers)
else: 
	if isData: process.p = cms.Path(process.myevents+process.myelectrons+process.mymuons+process.myphotons+process.myjets+process.mymets+process.mytaus+process.mytrigEvent+process.mypvertex+process.mytracks+process.mygenparticle+process.mytriggers)
	else: process.p = cms.Path(process.selectedHadronsAndPartons * process.jetFlavourInfosAK5PFJets * process.myevents+process.myelectrons+process.mymuons+process.myphotons+process.myjets+process.mymets+process.mytaus+process.mytrigEvent+process.mypvertex+process.mytracks+process.mygenparticle+process.mytriggers)
