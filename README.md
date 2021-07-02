# Objects Information Extractor Tool
## Description
This page contains instructions and examples of a way to extract information from a CMS root file type EDM.
The root files contain collections of physical objects arranged within this structure, these objects can be: electrons, muons, photons, jets, tracks, etc.
Scientific research requires specific information to obtain good results. This information can be extracted from an EDM file by creating an EDAnalyzer.

The instructions to write your own EDAnalyzer are [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule).

In the code there are several comments explaining the logic that follows, also, some twiki links were added to expand the knowledge on the subject.

## Usage instrucctions

1. Install CERN [virtual machine](http://opendata.cern.ch/docs/cms-virtual-machine-2011) or [Docker container](https://opendata.cern.ch/docs/cms-guide-docker) from the CMS open data website.

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
6. **Only if your are using the VM** and not the Docker container, make symbolic links to the conditions database.  In fact if you are not extracting information about corrections (e.g., jet corrections) or the trigger, you could actually comment out all the lines that have to do with the *global tag* (conditions database).
```
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db FT_53_LV5_AN1_RUNA.db
```
7. Run the CMSSW configuration file:
```
cmsRun poet_cfg.py
```


**As a result you will get a `myoutput.root` file with simple variables.**  This file has been divided by subdirectories, corresponding to each type of physics object (or kind of information).  If needed, in the `PhysObjectExtractor/test` directory, a `ROOT` analysis example code can be found.  It runs over the families of ROOT trees in the the `myoutput.root` file.



## References:

Muons

* [Muon IDs and Isolation](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId)
* [Muon Analysis](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis)

Taus

* [Tau IDs](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPFTauTagging#Legacy_Tau_ID_Run_I)
* [Tau Discriminators](https://twiki.cern.ch/twiki/bin/view/CMSPublic/NutShellRecipeFor5312AndNewer)

Jets

* [JEC+JER main documentation](https://arxiv.org/abs/1607.03663.pdf): Includes jet pt resolution functions and scale factors (Fig 41) applied for JER
* [Application of Hybrid Smearing Method](https://github.com/cms-legacydata-analyses/PhysObjectExtractorTool/blob/master/PhysObjectExtractor/src/JetAnalyzer.cc#L243-L270)
* [JER Scale Factor Instructions](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideJetResolution)
* [JER Scale Factor Lookup Table](https://github.com/adrager/cmssw/blob/JetResolution53/CondFormats/JetMETObjects/data/JetResolutionInputAK5PF.txt)
* [Read in JEC+JER Code](https://github.com/cms-legacydata-analyses/PhysObjectExtractorTool/blob/master/PhysObjectExtractor/src/JetAnalyzer.cc#L119-L143)
* [Apply JEC+JER Code](https://github.com/cms-legacydata-analyses/PhysObjectExtractorTool/blob/master/PhysObjectExtractor/src/JetAnalyzer.cc#L226-L270)
* [Produce JEC Text Files Code](https://github.com/cms-legacydata-analyses/PhysObjectExtractorTool/blob/master/PhysObjectExtractor/JEC/jec_cfg.py)
* [Read JEC+JER in Config](https://github.com/cms-legacydata-analyses/PhysObjectExtractorTool/blob/master/PhysObjectExtractor/python/poet_cfg.py#L125-L141)

* [BTagging CSV.csv File/Table](https://twiki.cern.ch/twiki/bin/view/CMSPublic/BtagRecommendation2011OpenData#Methods_to_Apply)
* [BTagging Help Understand CSV.csv](https://twiki.cern.ch/twiki/bin/view/CMSPublic/BTagCalibration) 
* [BTagging Weight Calculating Methods](https://twiki.cern.ch/twiki/bin/view/CMSPublic/BtagRecommendation2011OpenData#Methods_to_Apply_b_Tagging_Effic)
* [Application of Scale Factor Lookup Tables](https://github.com/cms-legacydata-analyses/PhysObjectExtractorTool/blob/master/PhysObjectExtractor/src/JetAnalyzer.cc#L335-L357)

Photons
* [Photon ID Paper](https://cms-physics.web.cern.ch/cms-physics/public/EGM-10-006-pas.pdf)
* [Photon ID Code](https://github.com/cms-legacydata-analyses/PhysObjectExtractorTool/blob/master/PhysObjectExtractor/src/PhotonAnalyzer.cc#L166-L242)

Electrons
* [Electron ID Twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/EgammaPublicData)
* [Electron ID Code](https://github.com/cms-legacydata-analyses/PhysObjectExtractorTool/blob/master/PhysObjectExtractor/src/ElectronAnalyzer.cc#L168-L212)

