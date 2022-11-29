# Physics Objects Extractor (PhysObjectExtractor) for 2015MiniAOD data

## Description
The `PhysObjectExtractor` package is the heart of the POET repository.  It contains a collection of [EDAnalyzers](https://cms-opendata-guide.web.cern.ch/cmssw/cmsswanalyzers/) that extract information from different physics objects into a [ROOT](https://cms-opendata-guide.web.cern.ch/tools/root/) file called `myoutput.root`.  They have been written separately for clarity and can be executed modularly using a single configuration file called `poet_cfg.py`.

The package is meant to resemble a [CMSSW](https://cms-opendata-guide.web.cern.ch/cmssw/cmsswoverview/) package.  It installs and runs in a similar way.

The `data` directory is reserved for files with information about the location of datafiles, data quality files, etc.

The `python` directory hosts the configuration file `poet_cfg.py` and possibly special versions of it.

The `src` directory hosts the C++ source code of all the different `EDAnalyzers` (and possibly some example filters) that can be configured using `python/poet_cfg.py`.  In the code, there are several comments explaining the logic.  Some links to appropiate references were added to expand the knowledge on the subject.

The `test` will contain analysis examples

## Usage instrucctions

1. Set up a [Docker container](https://opendata.cern.ch/docs/cms-guide-docker) and start it interactively (you could also work with a [virtual machine](https://opendata.cern.ch/docs/cms-virtual-machine-2015)):
  ```
  docker run -it --name my_od -P -p 5901:5901  -p 6080:6080 cmsopendata/cmssw_7_6_7-slc6_amd64_gcc493 /bin/bash
  ```

2. Obtain the code from git (this repository):
  ```
  git clone -b 2015MiniAOD https://github.com/cms-opendata-analyses/PhysObjectExtractorTool.git
  cd PhysObjectExtractorTool
  ```

3. Compile everything:
  ```
  cd PhysObjectExtractor
  scram b
  ```

4. Run the CMSSW job with the configuration file:
  ```
  cmsRun python/poet_cfg.py <isData>
  ```

  `<isData>` (to run on Data or not) is an optional boolean argument (default is False, i.e., runs over MC simulations)

5. As a result you will get myoutput.root file with simple variables. To check the contents of the .root file graphically, use ROOT. To open a graphics window, start the VNC application with:
  ```
  start_vnc
  ```

  Open the http link given in the message in a browser window of your local computer and connect using the password `cms.cern`. Alternatively, install a VNC viewer on your local computer and connect through it.

  Open the file with ROOT and start a graphics browser:
  ```
  root myoutpoot.root
  TBrowser t
  ```

6. Quit ROOT by typing `.q` in the terminal or selecting `Quit ROOT` in the `Broswer`menu of the graphical window. 

  Stop the VNC application before exiting the container:
  ```
  stop_vnc
  ```










