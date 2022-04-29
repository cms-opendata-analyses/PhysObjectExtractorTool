# Physics Objects Extractor (PhysObjectExtractor) for 2015MiniAOD data

## Description
The `PhysObjectExtractor` package is the heart of the POET repository.  It contains a collection of [EDAnalyzers](https://cms-opendata-guide.web.cern.ch/cmssw/cmsswanalyzers/) that extract information of different physics objects into a [ROOT](https://cms-opendata-guide.web.cern.ch/tools/root/) file called `myoutput.root`.  They have been written separately for clarity and can be executed modularly using a single configuration file called `poet_cfg.py`.

The package is meant to resemble a [CMSSW](https://cms-opendata-guide.web.cern.ch/cmssw/cmsswoverview/) package.  It installs and runs in a similar way.

The `data` directory is reserved for files with information about the location of datafiles, data quality files, etc.

The `python` directory hosts the configuration file `poet_cfg.py` and possibly special versions of it.

The `src` directory hosts the C++ source code of all the different `EDAnalyzers` that can be configured using `python/poet_cfg.py`.  In the code, there are several comments explaining the logic.  Some links to appropiate references were added to expand the knowledge on the subject.

The `test` directory contains `ROOT` analysis templates/examples that can be used to analyze the resulting `myoutput.root` files if needed.

## Usage instrucctions

1. Set up a [Docker container](https://opendata.cern.ch/docs/cms-guide-docker) and start it (you could also work with a [virtual machine](https://opendata.cern.ch/docs/cms-virtual-machine-2015)):
```
docker run -it --name my_od -P -p 5901:5901 cmsopendata/cmssw_7_6_7-slc6_amd64_gcc493 /bin/bash
```

1. Obtain the code from git (this repository):
```
git clone -b 2015MiniAOD https://github.com/cms-opendata-analyses/PhysObjectExtractorTool.git
cd PhysObjectExtractorTool
```

1. Compile everything:
```
cd PhysObjectExtractor
scram b
```

1. For convenience, make a soft link to the python configuration file
```
ln -s python/poet_cfg.py .
```

1. Run the job:
```
cmsRun poet_cfg.py
```

1. As a result you will get myoutput.root file with simple variables.  If using the Windows OS, to check the contents of the .root file graphically, install a vnc viewer (e.g TigerVNC Viewer) and run:
```
start_vnc
//When prompted enter a password
root myoutpoot.root
TBrowser t
```









