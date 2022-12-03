# Physics Objects Extractor (PhysObjectExtractor) for 2012 (AOD) data

## Description

The `PhysObjectExtractor` package is the heart of the POET repository.  It contains a collection of [EDAnalyzers](https://cms-opendata-guide.web.cern.ch/cmssw/cmsswanalyzers/) that extract information of different physics objects into a [ROOT](https://cms-opendata-guide.web.cern.ch/tools/root/) file called `myoutput.root`.  They have been written separately for clarity and can be executed modularly using a single configuration file called `poet_cfg.py`.

The package is meant to resemble a [CMSSW](https://cms-opendata-guide.web.cern.ch/cmssw/cmsswoverview/) package.  It installs and runs in a similar way.

The `data` directory is reserved for files with information about the location of datafiles, data quality files, etc.

The `JEC` directory hosts text files with energy corrections for jets. These were obtained executing the `jec_cfg.py` config file that resides within.

The `python` directory hosts the configuration file `poet_cfg.py` and possibly special versions of it.

The `src` directory hosts the C++ source code of all the different `EDAnalyzers` that can be configured using `python/poet_cfg.py`.  In the code, there are several comments explaining the logic.  Some links to appropiate references were added to expand the knowledge on the subject.

The `test` directory contains `ROOT` analysis templates/examples that can be used to analyze the resulting `myoutput.root` files if needed.

The `condor` directory hosts scripts that are useful if an HTCondor cluster can be used to run containerized jobs.  These scripts work from the `PhysObjectExtractorTool/PhysObjectExtractor` level and have to be modified according to particular needs.  




## Usage instrucctions

1. Install CERN [virtual machine](http://opendata.cern.ch/docs/cms-virtual-machine-2011) or [Docker container](https://opendata.cern.ch/docs/cms-guide-docker) from the CMS open data website.

2. In the container, the environment is already set up and you can skip this step (and **do not** run `cmsenv`). In the virtual machine, set up your enviroment

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

3. Obtain the code from git (in the VM, do this in the "Outer shell" in the CMSSW_5_3_32/src/ directory):
```
git clone -b 2012 https://github.com/cms-opendata-analyses/PhysObjectExtractorTool.git
```
4. Compile everything (in the VM, compile and run in the "CMS shell"):
```
cd PhysObjectExtractorTool/PhysObjectExtractor
scram b
```
5. Make a soft link to the python configuration file
```
ln -s python/poet_cfg.py .
```
6. **Only if your are using the VM** and not the Docker container, make symbolic links to the [conditions databases](https://opendata.cern.ch/docs/cms-guide-for-condition-database) for data and MC. 
```
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL FT53_V21A_AN6
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL.db FT53_V21A_AN6_FULL.db
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL FT53_V21A_AN6_FULL
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_V27 START53_V27
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_V27.db START53_V27.db
```
In addition (only if you are using the VM), in the configuration file `poet_cfg.py`, remove `_data_stripped` and `MC_stripped` in the database names in the `process.GlobalTag.connect` lines, and remove `_FULL` from the GlobalTag name for data. Note that if you are not extracting information about corrections (e.g., jet corrections) or the trigger, you could actually comment out all the lines that have to do with the *global tag* (conditions database).

7. Run the CMSSW configuration file:
```
cmsRun poet_cfg.py <isData> <doPat>
```

`<isData>` (to run on Data or not) is an optional boolean argument (default is `False`, i.e., runs over MC simulations)

`<doPat>` (to run Physics Analysis Tool routines) is an optional boolean argument (default is `False`)


The result will be a `myoutput.root` file with simple variables.  The information in this file is organized in `ROOT` directories corresponding to each type of physics object (or kind of information).  In order to further analyze the data in these file, the analysis templates in the `test` folder may be used.
