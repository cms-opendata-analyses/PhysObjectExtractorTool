#!/bin/bash

cd ${_CONDOR_SCRATCH_DIR}
cat <<'EndOfFile' > execute.sh
#!/bin/bash

# Exit on error
set -e

echo "### Begin of job"

ID=$1
echo "ID:" $ID

PROCESS=$2
echo "Process:" $PROCESS

FILE=$3
echo "File:" $FILE

ISDATA=True
DOPAT=True

EOS_HOME=/eos/user/e/ecarrera
#EOS_HOME=/eos/user/FIRST_LETTER/USERNAME
echo "EOS home:" $EOS_HOME

OUTPUT_DIR=${EOS_HOME}/dpoa/workshop2021/
#OUTPUT_DIR=${EOS_HOME}/opendata_files/
echo "Output directory:" $OUTPUT_DIR

CMSSW_BASE=/afs/cern.ch/work/e/ecarrera/dpoa/workshop2021/CMSSW_5_3_32
#CMSSW_BASE=/afs/cern.ch/work/FIRST_LETTER/USERNAME/CMSSW_5_3_32
echo "CMSSW base:" $CMSSW_BASE

 CONFIG=${CMSSW_BASE}/src/PhysObjectExtractorTool/PhysObjectExtractor/python/poet_cfg.py

echo "CMSSW config:" $CONFIG

echo "Hostname:" `hostname`

echo "How am I?" `id`

echo "Where am I?" `pwd`

echo "What is my system?" `uname -a`

echo "### Start working"

# Trigger auto mount of EOS
ls -la $EOS_HOME

# Make output directory
mkdir -p ${OUTPUT_DIR}/${PROCESS}

# Setup CMSSW
THIS_DIR=$PWD
cd $CMSSW_BASE
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $THIS_DIR

# Copy config file
mkdir -p configs/
CONFIG_COPY=configs/cfg_${ID}.py
cp $CONFIG $CONFIG_COPY


#Decide if it is data or simulations
if [[ ${FILE} == *"Run2012"* ]]; then
     ISDATA=True
else
     ISDATA=False
fi

echo "ISDATA:" $ISDATA

echo "DOPAT:" $DOPAT


# Modify CMSSW config to run only a single file
sed -i -e "s,^files =,files = ['"${FILE}"'] #,g" $CONFIG_COPY
sed -i -e 's,^files.extend,#files.extend,g' $CONFIG_COPY

# Modify CMSSW config to read lumi mask from EOS
sed -i -e 's,data/Cert,'${CMSSW_BASE}'/src/PhysObjectExtractorTool/PhysObjectExtractor/data/Cert,g' $CONFIG_COPY

# Modify config to write output directly to EOS
sed -i -e 's,myoutput.root,'${PROCESS}_${ID}.root',g' $CONFIG_COPY

# Print config
cat $CONFIG_COPY

# Run CMSSW config
cmsRun $CONFIG_COPY $ISDATA $DOPAT

# Copy output file
xrdcp -f ${PROCESS}_${ID}.root root://eosuser.cern.ch/${OUTPUT_DIR}/${PROCESS}/${PROCESS}_${ID}.root
rm ${PROCESS}_${ID}.root

echo "### End of job"

EndOfFile

# Trigger auto mount of EOS
ls -la $EOS_HOME

# Make file executable
chmod +x execute.sh

export SINGULARITY_CACHEDIR="/tmp/$(whoami)/singularity"
singularity run -B /afs -B /eos -B /cvmfs -B /usr/libexec/condor -B /pool --no-home docker://cmssw/slc6:latest $(echo $(pwd)/execute.sh $1 $2 $3)
