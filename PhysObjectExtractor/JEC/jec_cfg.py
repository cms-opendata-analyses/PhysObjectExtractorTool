import FWCore.ParameterSet.Config as cms

import sys

if len(sys.argv) > 2: isData = bool(eval(sys.argv[2]))
print 'Writing JEC text files. isData = ',isData

# CMS process initialization
process = cms.Process('jecprocess')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# connect to global tag
# use "connect" lines only if you are getting access to the conditions database through CVMFS snapshot files (requires installing CVMFS client or using VM)
# if isData:
##    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1.db') # use this in the VM
#    process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
# else:
##    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1') # use this in the VM
#    process.GlobalTag.globaltag = 'START53_LV6A1::All'
#---- If the container has local DB files available, uncomment the four lines below
#---- instead of the corresponding lines above
if isData:
    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_data_stripped.db')
    process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
else:
    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1_MC_stripped.db')
    process.GlobalTag.globaltag = 'START53_LV6A1::All'


# setup JetCorrectorDBReader 
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.source = cms.Source('EmptySource')
process.ak5 = cms.EDAnalyzer('JetCorrectorDBReader', 
                             payloadName=cms.untracked.string('AK5PF'),
                             printScreen=cms.untracked.bool(False),
                             createTextFile=cms.untracked.bool(True))

if isData:
    process.ak5.globalTag = cms.untracked.string('FT_53_LV5_AN1')
else:
    process.ak5.globalTag = cms.untracked.string('START53_LV6A1')

process.p = cms.Path(process.ak5)
