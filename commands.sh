#!/bin/sh -l
# exit on erro
set -e
sudo chown $USER /mnt/vol

# For the plain github action with docker, the area would be available in /mnt/vol
cp -r /mnt/vol PhysObjectExtractorTool
cd PhysObjectExtractorTool
cd PhysObjectExtractor
scram b

# sed -i "s/process.GlobalTag.connect/#process.GlobalTag.connect/g" python/poet_cfg.py
cmsRun python/poet_cfg.py

cp *.root /mnt/vol/
echo ls -l /mnt/vol
ls -l /mnt/vol
