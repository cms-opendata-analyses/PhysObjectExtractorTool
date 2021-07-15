#!/usr/bin/env python

import ROOT
from ROOT import TFile, TDirectory, TTree
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
import os
import sys
import re

def parse_arguments():
    if not len(sys.argv) == 2:
        raise Exception("Run with './merge_jobs.py path/to/input/directory'.")
    return sys.argv[1]

def main(input_dir):
    # Sanitize input dir string
    if input_dir[-1] == "/":
        input_dir = input_dir[:-1]

    # Sanitize input directory
    print("Input directory: %s"%(input_dir))
    if not os.path.exists(input_dir):
        raise Exception("Input directory does not exist: %s"%(input_dir))
    if not os.path.isdir(input_dir):
        raise Exception("Input is no directory: %s"%(input_dir))

    # Extract process from path
    process = os.path.basename(input_dir)
    print("Process: %s"%(process))

    # Get expected number of files
    if not os.path.exists("data/") or not os.path.isdir("data/"):
        raise Exception("Directory \"data\" does not exist.")
    count_expected = 0
    count_line = 0
    filelist_combined = {}
    for f in os.listdir("data/"):
        if process in f:
            filelist = open(os.path.join("data", f)).readlines()
            count_expected += len(filelist)
            for line in filelist:
                filelist_combined[count_line] = line.rstrip()
                count_line += 1
    print("Expect %u files in input directory."%(count_expected))

    # Go through files and find missing ones
    files = {}
    for f in os.listdir(input_dir):
        if f[0] == ".": continue
        if not process+"_" in f:
            raise Exception("File %s does not match job file."%(f))
        n = re.search("%s_(.*).root"%(process), f).group(1)
        files[int(n)] = os.path.join(input_dir, f)

    missing_file = False
    argument_list = []
    for i in range(count_expected):
        if not i in files:
            argument_list.append("%u %s %s"%(i, process, filelist_combined[i]))
            print("Miss file with ID %u."%(i))
            missing_file = True
    print("Found %u files of %u expected files in input directory."%(len(files), count_expected))

    """
    # Try to open files and see whether they are corrupted
    count_zombies = 0
    for i in files:
        tfile = ROOT.TFile(files[i])
        if tfile.IsZombie():
            argument_list.append("%u %s %s"%(i, process, filelist_combined[i]))
            print("Found zombie file with ID %u."%(i))
            missing_file = True
            count_zombies += 1
        tfile.Close()
    print("Found %u zombie files of %u files in input directory."%(count_zombies, len(files)))
    """

    if missing_file:
        path_list = "arguments.txt"
        out_list = open(path_list, "w")
        for a in argument_list:
            out_list.write(a+"\n")
        #raise Exception("Found missing files, wrote arguments list to %s."%(path_list))


    #https://root-forum.cern.ch/t/moving-ttrees-into-tdirectories/25386/8
    #Extract directories and trees structure
    oldf = TFile.Open(files.values()[0],"READ")
    treename = "Events"
    mychains = []
    mydirs = []
    for k in oldf.GetListOfKeys():
        mydirs.append(k.ReadObj().GetName())
        thechain = ROOT.TChain(k.ReadObj().GetName()+"/"+treename)
        mychains.append(thechain)

    # Merge trees in directories
    mycounter = 1
    for f in files.values():
        if((mycounter%1.)==0):
            print ("Merging file # "+str(mycounter)+" "+f)
        for ch in mychains:
            ch.Add(f)
        mycounter+=1
    
    #clone chains and write to a new file
    output_path = os.path.join(input_dir.replace(process, ""), process+".root")
    newf = TFile.Open(output_path,'RECREATE')

    #create directories in new file:
    for d in mydirs:
        newf.cd()
        newf.mkdir(d)

    mycounter = 1
    #Clone chain in tree for corresponding dir in new file:
    for d,ch in zip(mydirs,mychains):
        if (mycounter%1.==0):
            print ("Cloning and writiing directory "+d)
        newf.cd()
        newf.cd(d)
        newtree = ch.CloneTree()
        newtree.Write()

    newf.Close()
    
    print("Wrote merged file to %s."%(output_path))

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
