"""Test Poet by histograming the extracted data.

The histograms and their expected mean value are defined in a 
yaml file. During the test creation the data file is read and 
histograms are filled.

The implementation follows the example at  https://docs.pytest.org/en/6.2.x/example/nonpython.html

Usage:

    pytest --input=myoutput.root --output=myhistogram.root --junit-xml=test_poet.junit.xml

Author:
    Dietrich Liko <Dietrich.Liko@oeaw.ac.at>
"""
import yaml
import ROOT

import pytest

# init ROOT

ROOT.gROOT.SetBatch()
ROOT.PyConfig.IgnoreCommandLineOptions = True
# ROOT.EnableImplicitMT()

def pytest_addoption(parser):
    """Add options to the pytest cmdline."""
    parser.addoption(
        "--input", action="store", default="myoutput.root", help="Input file"
    )
    parser.addoption(
        "--output", action="store", default="myhistograms.root", help="Output file"
    )
    parser.addoption(
        "--data", action="store_true", default=False, help="Test for real data"
    )

def pytest_configure(config):
    """Read the options and store them in global variables."""
    global input_file
    global output_file
    global data_flag
    input_file = config.getoption("--input")
    output_file =  config.getoption("--output")
    data_flag =  config.getoption("--data")


def pytest_collect_file(parent, path):
    """Read all yaml file with the histogram definitions."""
    if path.ext == ".yaml" and path.basename.startswith("test"):
        return YamlFile.from_parent(parent, fspath=path)

class YamlFile(pytest.File):
    """Process a yaml file.
    
    The datafile is opened and histograms are defined. Then
    the datafile is read and histograms are filled. 

    Finally tests are created to verify the histogram content 
    and returned to the framework.
    """
    def collect(self):
        """Process yaml file."""
        data = yaml.safe_load(self.fspath.open())

        # open datafile
        trees = data['trees']
        input = ROOT.TFile(f"file:{input_file}", "READ")
        events = input.Get(f"{trees[0]}/Events")
        friends = []
        for tree in trees[1:]:
            friends.append(input.Get(f"{tree}/Events"))
            events.AddFriend(friends[-1])

        # define RDataFrame and Histograms
        df = ROOT.RDataFrame(events)
        event_counter = df.Count()

        histos = {}
        for name, h1def in data['histo1D'].items():
            print('Data', data)
            if data_flag and 'GenPart' in name:
                continue
            title = h1def['title'] if 'title' in h1def else name
            bins = h1def['bins']
            histos[name] = df.Histo1D((name, title, *bins), name)

        # fill histograms
        # ROOT.RDF.RunGraphs(histos.values())

        # write histograms
        output = ROOT.TFile(output_file, "RECREATE")
        for h1d in histos.values():
            h1d.Write()
        output.Close()

        # create a test for the number of events read
        if 'events' in data:
            yield EventCounterItem.from_parent(self, name="Events", event_counter=event_counter, value=data['events'])

        # create a test for the mean value of the various histograms
        for name, h1def in data['histo1D'].items():
            if data_flag and 'GenPart' in name:
                continue
            yield Histo1DItem.from_parent(self, name=name, histo=histos[name], range=h1def['test_mean'])

class EventCounterItem(pytest.Item):
    """Test for events read."""

    def __init__(self, name, parent, event_counter, value):
        super().__init__(name, parent)
        self.event_counter = event_counter
        self.value = value

    def runtest(self):

        assert self.event_counter.GetValue() == self.value

class Histo1DItem(pytest.Item):
    """Test for histogram mean value."""
    def __init__(self, name, parent, histo, range):
        super().__init__(name, parent)
        self.histo = histo
        self.mean_min, self.mean_max = range

    def runtest(self):

        assert self.mean_min < self.histo.GetMean()
        assert self.histo.GetMean() < self.mean_max

