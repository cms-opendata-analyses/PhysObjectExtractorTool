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
        'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/0005D1FB-4BCF-E311-9FE4-002590A8312A.root',
        'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/002B0269-0DC7-E311-8AC7-002481E150DA.root',
        'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/002D845F-65C7-E311-96FE-001E6739689C.root',
        'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/002FB2D4-F6C7-E311-8D85-002590A3C95E.root',
        'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/00332A5F-6DC8-E311-9955-002481E154CE.root',
        'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/00345E78-7CC5-E311-A9C7-001E6739730A.root',
        ))

# Set global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db') # use this in the VM
# process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1_MC_stripped.db') # use this in the container
process.GlobalTag.globaltag = "START53_LV6A1::All"

# Load PAT configs and build some light sequences
process.load('PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff')
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets
patJets.discriminatorSources = cms.VInputTag(
    cms.InputTag("combinedSecondaryVertexBJetTags"),
    cms.InputTag("combinedSecondaryVertexMVABJetTags"),
    cms.InputTag("jetBProbabilityBJetTags"),
    cms.InputTag("jetProbabilityBJetTags"),
    cms.InputTag("simpleSecondaryVertexHighEffBJetTags"),
    cms.InputTag("simpleSecondaryVertexHighPurBJetTags"),
    cms.InputTag("trackCountingHighEffBJetTags"),
    cms.InputTag("trackCountingHighPurBJetTags"),
    )
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

