## Higgs to Tau Tau EventLoopAnalysis example

This is an example of a simplified analysis (mostly used for outreach purposes) where we use the EventLoopAnalysisTemplate in this package.  The [original analysis](http://opendata.web.cern.ch/record/12350) uses the ROOT's [RDataFrame framework](https://root.cern/doc/master/classROOT_1_1RDataFrame.html).  The version here presents a recast of the original example for pedagogical purposes.

The files used by this example are listed below:

```
root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/Run2012B_TauPlusX.root
root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/Run2012C_TauPlusX.root
root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/DYJetsToLL.root
root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/TTbar.root
root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/W1JetsToLNu.root
root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/W2JetsToLNu.root
root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/W3JetsToLNu.root
root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/GluGluToHToTauTau.root
root://eospublic.cern.ch//eos/opendata/cms/upload/od-workshop/ws2021/VBF_HToTauTau.root
```

These ntuples were obtained with the `PhysObjectExtractor`(POET) package using the version of the configuration file `poet_cfg.py` in this directory, which is present for reference.

The `EventLoopAnalysis` example code can be compiled with

`g++ -std=c++11 -g -O3 -Wall -Wextra -o EventLoopAnalysis EventLoopAnalysis.cxx $(root-config --cflags --libs) -lGenVector`

This analysis code loops over the files in the list above (signal, background and data) and produces a `histograms.root` file with a collection of example histograms.  Then, the `plot.py` can be used to make nice, colored plots of those histograms.


## Usage

1. Compile with (requires ROOT 5.32 or newer):
`g++ -std=c++11 -g -O3 -Wall -Wextra -o EventLoopAnalysis EventLoopAnalysis.cxx $(root-config --cflags --libs) -lGenVector`


1. Run with:
`./EventLoopAnalysis`
to obtain a `histograms.root` file.

2. Make some example plots with
`python plot.py`
