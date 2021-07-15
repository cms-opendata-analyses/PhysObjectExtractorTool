# Condor Scripts

## Description

These scripts can be used to run `PhysObjectExtractor` in an HTCondor cluster if it has the capability to run containerized applications.

General instructions to use them are as follows:

* Copy the files in this `condor` folder to the `PhysObjectExtractorTool/PhysObjectExtractor`

* Uncomment the `FileUtils` way of handling files in the `PoolSource` in the `poet_cfg.py` file.

* Test the job locally and then change number of events to `-1` (all events)

* Everything is executed with `./submit_jobs.sh`.  This script calls `create_job.py` to create an HTCondor submission file, `job.jdl` based on the `job.sh` template.  There are a few variables that need to be adjusted for you own situation.

* The `merge_jobs.py` checks whether all the files were processed and then merge the root files into a single file.  Merging can fail if the files are too big.