# Objects Information Extractor Tool
## Description
This page contains instructions and examples of a way to extract information from 2015 MiniAOD file

## Usage instrucctions

* Get the slc6 CMSSW_7_6_7 image:
```
docker pull gitlab-registry.cern.ch/cms-cloud/cmssw-docker/cmssw_7_6_7-slc6_amd64_gcc493:2021-08-11-75df77cf
```
* Add vnc to the image:
```
cd cmssw-docker
git checkout vnc
cd vnc
vi Docker file 
//replace the IMAGEBASE argument with gitlab-registry.cern.ch/cms-cloud/cmssw-docker/cmssw_7_6_7-slc6_amd64_gcc493:2021-08-11-75df77cf
docker build .
```
* Get the MiniAOD POET:
```
cd ../../
mkdir MiniAOD
cd MiniAOD
git clone https://github.com/apetkovi1/PhysObjectExtractorTool.git
cd PhysObjectExtractorTool
git checkout MiniAOD
```
* Create a container based on slc6 CMSSW_7_6_7 image with vnc installed and mount MiniAOD local volume:
```
cd ../../
docker run -it --name slc6_MiniAOD -P -p 5901:5901 -v ${HOME}/MiniAOD:/home/cmsusr/MiniAOD  [put your image ID here] /bin/bash
//run all the following commands from the container
ln -s ~/MiniAOD/   

```
* Create a CMS runtime variables and compile everything:
```
cmsenv
scram b
```
* Run the configuration file:
```
cmsRun MiniAOD/PhysObjectExtractorTool/PhysObjectExtractor/python/poet_cfg.py
```
* As a result you will get myoutput.root file with simple variables. To check the contents of the .root file graphically, install a vnc viewer (e.g TigerVNC Viewer) and run:
```
start_vnc
//When prompted enter a password
root myoutpoot.root
TBrowser tb
```









