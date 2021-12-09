# Objects Information Extractor Tool
## Description
This page contains instructions and examples of a way to extract information from 2015 MiniAOD file

## Usage instrucctions

You can run the code in the container if you have docker available:

* Start the container
```
docker run -it --name my_od -P -p 5901:5901 cmsopendata/cmssw_7_6_7-slc6_amd64_gcc493 /bin/bash
```

* Get this code:
```
git clone https://github.com/cms-legacydata-analyses/PhysObjectExtractorTool.git
cd PhysObjectExtractorTool
git checkout 2015MiniAOD
```

* Compile everything:
```
scram b
```
* Run the job:
```
cmsRun PhysObjectExtractorTool/PhysObjectExtractor/python/poet_cfg.py
```
* As a result you will get myoutput.root file with simple variables. To check the contents of the .root file graphically, install a vnc viewer (e.g TigerVNC Viewer) and run:
```
start_vnc
//When prompted enter a password
root myoutpoot.root
TBrowser t
```









