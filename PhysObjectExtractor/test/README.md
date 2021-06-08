## Poet Analysis Template

This is a template to perform analysis over the `myoutput.root` type of files that are obtained with the POET tool, which is a simple AOD skimmer.
  
The code file has actually two template functions.  One based on the traditional event loop over ROOT trees and the other one with a more sophisticated columnar analysis approach.  

This template needs ROOT version greater than 6.18.

The code can be compiled with

`/g++ -g -O3 -Wall -Wextra -Wpedantic -o PoetAnalysis PoetAnalysisTemplate.cxx $(root-config --cflags --libs)`


Since the `myoutput.root` files may contain different directories for different 
kind of physics objects, the `AddFriend` method (`TFriendElement`) from the `TTree` ROOT class is used.  This allows for keeping things simple and modular.

##Usage

1. Compile with: 
`/g++ -g -O3 -Wall -Wextra -Wpedantic -o PoetAnalysis PoetAnalysisTemplate.cxx $(root-config --cflags --libs)`

1. Run like:
./PoetAnalysis

Two equivalent histogram files will be created, `myhistograms.root` (from the event loop analysis) and the `myRDFhistograms.root` (from the RDF analysis).
