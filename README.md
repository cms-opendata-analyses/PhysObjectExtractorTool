# Objects Information Extractor Tool
## Description
This page contains instructions and examples of a way to extract information from MiniAOD file

## Usage instrucctions

* Create a project area:
```
cmsrel CMSSW_7_6_7
```
* Change to the CMSSW_7_6_7/src/ directory:

```
cd CMSSW_7_6_7/src/
```
* Then, run the following command to create the CMS runtime variables:

```
cmsenv
```
3. Obtain the code from git:
```
git clone  https://github.com/apetkovi1/PhysObjectExtractorTool.git
cd PhysObjectExtractorTool
git checkout MiniAOD
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
6. Run the CMSSW configuration file:
```
cmsRun poet_cfg.py
```

##### As a result you will get a myoutput.root file with simple variables.








