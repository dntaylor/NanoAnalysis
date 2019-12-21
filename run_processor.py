#!/usr/bin/env python
import os
import json
import argparse
import glob
from tqdm import tqdm

import uproot
import numpy as np
from coffea import hist, processor
from coffea.util import load, save
from dask.distributed import Client, LocalCluster

from hzzProcessor import HZZProcessor

processor_map = {
    'HZZ': HZZProcessor,
}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run coffea file')
    parser.add_argument('processor', type=str, default=None, help='The analysis, either precompiled or a string in processor_map')
    parser.add_argument('--year', choices=['2016','2017','2018'], default='2018', help='Data taking year (default: %(default)s)')
    parser.add_argument('--output', default='hists.coffea', help='Output histogram filename (default: %(default)s)')
    parser.add_argument('-j', '--workers', type=int, default=1, help='Number of workers to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
    parser.add_argument('--dask', action='store_true', help='Use dask to distribute')
    parser.add_argument('--local', action='store_true', help='Use dask, but with local cluster')
    parser.add_argument('--test', action='store_true', help='Only process a few files')
    args = parser.parse_args()

    if args.dask:
        if args.local:
            cluster = LocalCluster(args.workers,local_dir='/scratch/dntaylor/dask-worker-space')
            client = Client(cluster)
        else:
            client = Client(os.environ['DASK_SCHEDULER'])

    redirector = 'root://cmsxrootd.fnal.gov/'
    #if args.test: redirector = 'root://cmsxrootd.hep.wisc.edu/'

    if args.test:
        fileset = {'DY': ['dy_2018.root'], 'HZZ': ['hzz_2018.root'], 'DoubleMuon': ['DoubleMuon_2018.root'],}
    else:
        with open('data.json') as f:
            fulldatafileset = json.load(f)
        with open('mc.json') as f:
            fullmcfileset = json.load(f)

        fileset = {}
        for s in fulldatafileset[args.year]:
            fileset[s] = []
            for d in fulldatafileset[args.year][s]['datasets']:
                if args.test and 'Run2018A' not in d: continue
                fileset[s] += [redirector+x for x in fulldatafileset[args.year][s]['files'][d]]
        for s in fullmcfileset[args.year]:
            fileset[s] = []
            for d in fullmcfileset[args.year][s]['datasets']:
                fileset[s] += [redirector+x for x in fullmcfileset[args.year][s]['files'][d]]

    #if args.test:
    #    for s in fileset: fileset[s] = fileset[s][:1]

    print('Will process: ',list(fileset.keys()))
    if args.test:
        print(fileset)


    if args.processor in processor_map:
        corrections = load(f'corrections_{args.year}.coffea')
        processor_instance = processor_map[args.processor](
            year=args.year,
            corrections=corrections,
        )
    elif os.path.exists(args.processor):
        # compiled wont pickle?
        processor_instance = load(args.processor)
    else:
        print(f'Cannot understand {args.processor}.')

    if args.dask:
        executor = processor.dask_executor
        executor_args = {'client': client, 'compression': 1, 'savemetrics': True, 'flatten': True,}
    else:
        executor = processor.futures_executor
        #executor = processor.iterative_executor
        executor_args = {'workers': args.workers,'flatten': True,}

    accumulator = processor.run_uproot_job(
        fileset,
        treename='Events',
        processor_instance=processor_instance,
        executor=executor,
        executor_args=executor_args,
        chunksize=200000,
    )
    
    save(accumulator, args.output)
    if args.test: print(accumulator)
