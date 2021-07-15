## Poet Analysis Template

This is a template to perform analysis over the `myoutput.root` type of files that are obtained with the POET tool, which is a simple AOD skimmer.
  
There are currently two examples. One based on the traditional event loop over ROOT trees called `EventLoopAnalysisTemplate.cxx` (built from the TTree MakeClass() method in ROOT version 5.32), and the other one with a more sophisticated columnar analysis approach called `RDFAnalysisTemplate.cxx` (which needs ROOT version > 6.22).  

The code can be compiled with

`g++ -g -O3 -Wall -Wextra -o EventLoopAnalysis EventLoopAnalysisTemplate.cxx $(root-config --cflags --libs)`

or

`g++ -g -O3 -Wall -Wextra -Wpedantic -o RDFAnalysis RDFAnalysisTemplate.cxx $(root-config --cflags --libs)`


Since the `myoutput.root` files may contain different directories for different 
kind of physics objects, the `AddFriend` method (`TFriendElement`) from the `TTree` ROOT class is used.  This allows for keeping things simple and modular.

## Usage

1. Compile with: 
`g++ -g -O3 -Wall -Wextra -o EventLoopAnalysis EventLoopAnalysisTemplate.cxx $(root-config --cflags --libs)`

or

`g++ -g -O3 -Wall -Wextra -Wpedantic -o RDFAnalysis RDFAnalysisTemplate.cxx $(root-config --cflags --libs)`

1. Run with:
./EventLoopAnalysis 
or
./RDFAnalysis

Equivalent histogram files will be created, `myhistograms.root` (from the event loop analysis) and the `myRDFhistograms.root` (from the RDF analysis).
