# B-tagging efficiency calculator

The scripts in this directory can be used to compute the efficiency for tagging b-flavour hadrons in jets, using any simulation sample. These efficiencies are needed in order to apply data/simulation scale factors to correct the b-tagging performance in simulation. The main documentation for OpenData b-tagging is on the [BTaggingRecommendation Twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/BtagRecommendation2011OpenData).

## Important files

 * `plugins/WeightAnalyzerBEff.cc`: this EDAnalyzer looks at each jet, determines its flavour (b, c, or light/gluon), and tests whether it passes pre-defined cuts on a certain b-tagging algorithm's discriminator distribution. The results are stored in histograms as a function of jet transverse momentum. 
 * `python/befficiency_patjets_cfg.py`: this configuration file is used to run the analyzer. You can configure several items:
   * `fileNames`: this is a comma-separated list of input file names that will be processed by the analyzer.  
   * `jetTag`: this parameter contains the name of the collection of jets you wish to analyze.
   * `discriminator`: this string is the name of the b-tagging algorithm you will study.
   * `DiscriminatorValueTight/Medium/Loose`: these float values represent the pre-defined cut points for "tight" (high purity, low background) through "loose" (high efficiency, higher background) settings of the b-tagging algorithm. See the reference TWiki linked above.
 *  `plotBEff.C`: this utility macro will produce plots and text dumps of the efficiencies from the histograms in the ROOT file

## Usage instrucctions

1. Install CERN [virtual machine](http://opendata.cern.ch/docs/cms-virtual-machine-2011) or [Docker container](https://opendata.cern.ch/docs/cms-guide-docker) from the CMS open data website.

2. Set up your enviroment (done for you in docker!)

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
3. If you have not done so already, obtain the code from git:
```
git clone git://github.com/cms-legacydata-analyses/PhysObjectExtractorTool.git
cd PhysObjectExtractorTool
```
4. Compile everything, including code in the `BTagging` folder:
```
scram b
```
5. Enter the BTagging directory
```
cd BTagging
```
6. **Only if your are using the VM** and not the Docker container, make symbolic links to the conditions database.  In fact if you are not extracting information about corrections (e.g., jet corrections) or the trigger, you could actually comment out all the lines that have to do with the *global tag* (conditions database).
```
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_V27 START53_V27
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_V27.db START53_V27.db
```
7. Run the CMSSW configuration file:
```
cmsRun python/befficiency_patjets_cfg.py
```
**As a result you will get a `flavortagefficiencies.root` file with histograms stored inside a folder called `mcweightanalyzer`. 

8. To run the plotting macro, execute: 
```
root -l -b -q plotBEff.C >& myefficiencies.log
```
This will save plots of the efficiencies as well as provide you with test lists of numerical values. 

9. Enter your efficiencies in the jet analyzers:
```
cd ../PhysObjectExtractor/src/
vi PatJetAnalyzer.cc # or JetAnalyzer.cc
```
Enter your efficiencies into the functions named `getBtagEfficiency`, `getCtagEfficiency`, and `getLFtagEfficiency`. Then recompile:
```
scram b
```



