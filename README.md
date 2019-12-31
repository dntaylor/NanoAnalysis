# NanoAnalysis
Analysis of CMS data using the NanoAOD data format and the coffea framework

## Setup

To setup the environment
```bash
bash setup.sh
```
The generated virtualenv will allow the jobs to be run 
on a distributed system.

To update the file lists
```bash
python get_datasets.py
python make_filesets.py
```
This will create a new directory structure with filesets
generated from DAS.

## Run analysis

Can run locally or using condor.
Distribution via parsl is the only supported condor mode right now.
```bash
./run_processor.py -j 200 --condor --parsl hzzProcessor 2018 filesets/2018/all.json
```

Then convert output to flat histograms for plotting with ROOT
```bash
./flatten.py hists/*/*/*.coffea
```

## Development

All processors (source and compiled) should be stored in the 
[processors](processors) directory.
Any corrections that need to be loaded should be added to the 
[corrections builder](corrections/build_corrections.py).
