#!/bin/sh -l
# Parameters: $1 number of events (default: x), $2 configuration file (default: physicsobjectsinfo_cfg.py)
sudo chown $USER /mnt/vol

mkdir workspace
cd workspace
# For the plain github action with docker, the area would be available in /mnt/vol
cp -r /mnt/vol PhysObjectExtractorTool
cd PhysObjectExtractorTool
cd PhysObjectExtractor
scram b

sed -i "s/process.GlobalTag.connect/#process.GlobalTag.connect/g" python/physicsobjectsinfo_cfg.py
cat python/physicsobjectsinfo_cfg.py
cmsRun python/physicsobjectsinfo_cfg.py

cp *.root /mnt/vol/
echo ls -l /mnt/vol
ls -l /mnt/vol
