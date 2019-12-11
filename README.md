# NanoAnalysis
Analysis of CMS data using the NanoAOD data format and the coffea framework

## Setup

To setup the environment
```bash
bash setup.sh
```

To update the file lists
```bash
python get_datasets.py
```

## Run analysis

Locally running only so far
```bash
./run_processor.py HZZ --test
```

The convert output to flat histograms for plotting with ROOT
```bash
python flatten.py
```
