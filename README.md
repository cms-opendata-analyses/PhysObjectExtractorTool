# Objects Information Extractor Tool
## Description
This page contains instructions and examples of a way to extract information from MiniAOD file

## Usage instrucctions

1. Install CERN [virtual machine](http://opendata.cern.ch/docs/cms-virtual-machine-2011) from the CMS open data website.
2. Set up your enviroment

* Create a project area:
```
cmsrel CMSSW_5_3_32
```
* Change to the CMSSW_5_3_32/src/ directory:

```
cd CMSSW_5_3_32/src/
```
* Then, run the following command to create the CMS runtime variables:

```
cmsenv
```
3. Obtain the code from git:
```
git clone git://github.com/cms-legacydata-analyses/PhysObjectExtractorTool.git
cd PhysObjectExtractorTool
```
4. Compile everything:
```
cd PhysObjectExtractor
scram b
```
5. Make a soft link to the python configuration file
```
ln -s python/poet_cfg.py .
```
6. Make symbolic links to the conditions database: (trigger analizer)
```
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db FT_53_LV5_AN1_RUNA.db
```
7. Run the CMSSW configuration file:
```
cmsRun poet_cfg.py
```

##### As a result you will get a myoutput.root file with simple variables.








