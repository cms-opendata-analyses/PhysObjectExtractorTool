#!/bin/sh -l
# Parameters: $1 number of events (default: -1 (all events)), $2 configuration file (default: physicsobjectsinfo_cfg.py)
sudo chown $USER /mnt/vol

# For the plain github action with docker, the area would be available in /mnt/vol
cp -r /mnt/vol PhysObjectExtractorTool
cd PhysObjectExtractorTool
cd PhysObjectExtractor
scram b

if [ -z "$1" ]; then nev=-1; else nev=$1; fi
if [ -z "$2" ]; then config=physicsobjectsinfo_cfg.py; else config=$2; fi
eventline=$(grep maxEvents python/$config)
sed -i "s/$eventline/process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32($nev) ))/g" python/$config

sed -i "s/process.GlobalTag.connect/#process.GlobalTag.connect/g" python/physicsobjectsinfo_cfg.py
cat python/physicsobjectsinfo_cfg.py
cmsRun python/physicsobjectsinfo_cfg.py

cp *.root /mnt/vol/
echo ls -l /mnt/vol
ls -l /mnt/vol
