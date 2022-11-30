# parameters: $1 runas, $2 config, $3 isData
# if running outside github actions, give any other 1st paramater than github i.e. commands.sh mywork (or omit the paramaters)
# defaults:
if [ -z "$1" ]; then runas=other; else runas=$1; fi
if [ -z "$2" ]; then config=python/poet_cfg.py; else config=$2; fi
if [ -z "$3" ]; then isData=False; else isData=$3; fi

set -e

# For the plain github action with docker, the area would be available in /mnt/vol
if [ $runas = github ]
then
  sudo chown $USER /mnt/vol
  cp -r /mnt/vol PhysObjectExtractorTool
else
  git clone -b 2015MiniAOD https://github.com/cms-opendata-analyses/PhysObjectExtractorTool.git
fi

cd PhysObjectExtractorTool
cd PhysObjectExtractor
scram b

echo Going to run $config with isData $isData 
cmsRun $config $isData 

if [ $runas = github ]
then
  cp *.root /mnt/vol/
  echo ls -l /mnt/vol
  ls -l /mnt/vol
fi  
