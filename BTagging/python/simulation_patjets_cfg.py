import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from PhysicsTools.PatAlgos.patTemplate_cfg import *
import os 

relBase = os.environ['CMSSW_BASE']

process = cms.Process("AOD2NanoAOD")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("AOD2NanoAOD")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# Set the maximum number of events to be processed (-1 processes all events)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

##### ------ This set of file setup is for LOTS of files, not a short test ------
# Define files of dataset
#files = FileUtils.loadListFromFile("data/CMS_MonteCarlo2012_Summer12_DR53X_TTbar_8TeV-Madspin_aMCatNLO-herwig_AODSIM_PU_S10_START53_V19-v2_00000_file_index.txt")
#files.extend(FileUtils.loadListFromFile("data/CMS_MonteCarlo2012_Summer12_DR53X_TTbar_8TeV-Madspin_aMCatNLO-herwig_AODSIM_PU_S10_START53_V19-v2_20000_file_index.txt"))

#process.source = cms.Source(
#   "PoolSource", fileNames=cms.untracked.vstring(*files))

##### ------- This is a test file
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
#'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/000A9D3F-CE4C-E311-84F8-001E673969D2.root',
#'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/000EA23D-084D-E311-901C-0025B3E06508.root',
#'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/002854E5-E94C-E311-BA6F-0025B3E05C60.root',
#'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/006E132C-FD4C-E311-B905-002590A88806.root',
#'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/02110278-274D-E311-B574-002590A3CA14.root'))
#'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/0244F8C4-1A4D-E311-B14B-001E67397B07.root',
#'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/02EED3E1-184D-E311-8B96-001E6739665D.root',
#'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/041E194D-CD4C-E311-86F4-0025B3E05E1C.root',
'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/0448AB4F-364D-E311-B1E7-001E67397CBA.root',
'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/04FCA1D5-E74C-E311-92CE-002590A887F0.root'))

# Set global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START53_V27::All"

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

# Number of events to be skipped (0 by default)
process.source.skipEvents = cms.untracked.uint32(0)

#-------------------------------------------------------------------------
# AK5 jets
#-------------------------------------------------------------------------
#from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
#process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone()

#from PhysicsTools.JetMCAlgos.AK5PFJetsMCFlavourInfos_cfi import ak5JetFlavourInfos
#process.jetFlavourInfosAK5PFJets = ak5JetFlavourInfos.clone()


################################
## MC Weights Analyzer
################################
process.mcweightanalyzer = cms.EDAnalyzer(
    "WeightAnalyzerBEff",
    jetTag = cms.InputTag("selectedPatJetsAK5PFCorr"),
    DiscriminatorValueTight = cms.double(0.898),
    DiscriminatorValueMedium = cms.double(0.679),
    DiscriminatorValueLoose = cms.double(0.244),
    )
     
process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("flavortagefficiencies.root"))

process.p = cms.Path(process.patDefaultSequence * process.mcweightanalyzer)

