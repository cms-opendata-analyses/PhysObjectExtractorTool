import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
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
        'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/TTbar_8TeV-Madspin_aMCatNLO-herwig/AODSIM/PU_S10_START53_V19-v2/00000/04FCA1D5-E74C-E311-92CE-002590A887F0.root'))

# Set global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START53_V27::All"

# Load PAT configs and build some light sequences
process.load('PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff')
process.load('PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi')
process.patDefaultSequence = cms.Sequence(process.makePatJets * process.selectedPatJets)
process.load('RecoJets.Configuration.RecoPFJets_cff')
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection, runBTagging
jetcorrlabels = ['L1FastJet','L2Relative','L3Absolute']

# Set up the new jet collection
process.ak5PFJets.doAreaFastjet = True
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PFCorr',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(jetcorrlabels)),
                 doType1MET   = False,
                 doL1Cleaning = False,
                 doL1Counters = False,
                 doJetID      = True,
                 jetIdLabel   = "ak5",
                 ) 

# Number of events to be skipped (0 by default)
process.source.skipEvents = cms.untracked.uint32(0)


################################
## MC Weights Analyzer
################################
process.mcweightanalyzer = cms.EDAnalyzer(
    "WeightAnalyzerBEff",
    jetTag = cms.InputTag("selectedPatJetsAK5PFCorr"),
    discriminator = cms.string("combinedSecondaryVertexBJetTags"),
    DiscriminatorValueTight = cms.double(0.898),
    DiscriminatorValueMedium = cms.double(0.679),
    DiscriminatorValueLoose = cms.double(0.244),
    )
     
process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("flavortagefficiencies.root"))

process.p = cms.Path(process.patDefaultSequence * process.mcweightanalyzer)

