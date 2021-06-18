import FWCore.ParameterSet.Config as cms

isData = False
#if len(sys.argv) > 1: isData = bool(eval(sys.argv[1]))
#print 'Writing JEC text files. isData = ',isData

# CMS process initialization
process = cms.Process('jecprocess')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# connect to global tag
if isData:
#    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL.db')
    process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'
else:
#    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_V27.db')
    process.GlobalTag.globaltag = 'START53_V27::All'


# setup JetCorrectorDBReader 
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.source = cms.Source('EmptySource')
process.ak5 = cms.EDAnalyzer('JetCorrectorDBReader', 
                             payloadName=cms.untracked.string('AK5PF'),
                             printScreen=cms.untracked.bool(False),
                             createTextFile=cms.untracked.bool(True))

if isData:
    process.ak5.globalTag = cms.untracked.string('FT53_V21A_AN6')
else:
    process.ak5.globalTag = cms.untracked.string('START53_V27')

process.p = cms.Path(process.ak5)
