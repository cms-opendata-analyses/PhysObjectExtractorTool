# Objects Information Extractor Tool
## Description
This page contains instructions and examples of a way to extract information from a CMS root file type EDM.
The root files contain collections of physical objects arranged within this structure, these objects can be: electrons, muons, photons, jets, tracks, etc.
Scientific research requires specific information to obtain good results. This information can be extracted from an EDM file by creating an EDAnalyzer.

The instructions to write your own EDAnalyzer are [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule).

In the code there are several comments explaining the logic that follows, also, some twiki links were added to expand the knowledge on the subject.

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








