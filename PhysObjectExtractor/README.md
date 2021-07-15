# Physics Objects Extractor (PhysObjectExtractor)

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
cmsRun poet_cfg.py <isData> <doPat>
```

`<isData>` (to run on Data or not) is an optional boolean argument (default is `False`, i.e., runs over MC simulations)

`<doPat>` (to run Physics Analysis Tool routines) is an optional boolean argument (default is `False`)


The result will be a `myoutput.root` file with simple variables.  The information in this file is organized in `ROOT` directories corresponding to each type of physics object (or kind of information).  In order to further analyze the data in these file, the analysis templates in the `test` folder may be used.



