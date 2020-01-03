#!/usr/bin/env python
import os
import json
import argparse
import glob
import logging

logging.basicConfig(filename='_run_processor.log', level=logging.INFO, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
#logging.basicConfig(level=logging.INFO, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
rootLogger = logging.getLogger()
logging.captureWarnings(True)

import uproot
import numpy as np
from coffea import hist, processor
from coffea.util import load, save

from utilities import python_mkdir

from parsl_config import parsl_condor_config, parsl_local_config

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run coffea file')
    parser.add_argument('baseprocessor', type=str, default=None, help='The analysis processor tag')
    parser.add_argument('year', type=str, default=None, help='The analysis year')
    parser.add_argument('fileset', default='', help='Fileset to process')
    parser.add_argument('--output', default='', help='Output histogram filename')
    parser.add_argument('-j', '--workers', type=int, default=1, help='Number of workers to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
    parser.add_argument('-n', '--maxfiles', type=int, default=-1, help='Maximum number of files to process')
    scheduler = parser.add_mutually_exclusive_group()
    scheduler.add_argument('--dask', action='store_true', help='Use dask to distribute')
    scheduler.add_argument('--parsl', action='store_true', help='Use parsl to distribute')
    parser.add_argument('--condor', action='store_true', help='Use distributed, but with condor')
    parser.add_argument('--debug', action='store_true', help='Add debug verbosity')
    args = parser.parse_args()

    dataset = os.path.splitext(os.path.basename(args.fileset))[0]
    processorpath = f'processors/{args.baseprocessor}_{args.year}.coffea'

    if not args.output:
        args.output = f'hists/{args.baseprocessor}/{args.year}/{dataset}.coffea'
    if os.path.dirname(args.output):
        os.makedirs(os.path.dirname(args.output), exist_ok=True)

    if os.path.exists(processorpath):
        processor_instance = load(processorpath)
    else:
        raise ValueError(f'Cannot understand {processorpath}.')

    #redirector = 'root://cms-xrd-global.cern.ch/'
    #redirector = 'root://xrootd-cms.infn.it/'
    redirector = 'root://cmsxrootd.fnal.gov/'
    #redirector = 'root://cmsxrootd.hep.wisc.edu/'

    with open(args.fileset) as f:
        fileset = json.load(f)
        for ds in fileset:
            fileset[ds] = [redirector+x if x.startswith('/store') else x for x in fileset[ds]]
            if args.maxfiles>0: fileset[ds] = fileset[ds][:args.maxfiles]

    rootLogger.info('Will process: '+' '.join(list(fileset.keys())))

    if args.dask:
        from dask.distributed import Client, LocalCluster
        if args.condor:
            pass
        else:
            cluster = LocalCluster(args.workers,local_dir='/scratch/dntaylor/dask-worker-space')
            client = Client(cluster)

    if args.parsl:
        import parsl
        if args.condor:
            htex = parsl_condor_config(workers=args.workers)
        else:
            htex = parsl_local_config(workers=args.workers)
        parsl.load(htex)
            

    executor_args = {
        'savemetrics': True, 'flatten':True, 
        'desc': f'Processing {args.baseprocessor} {args.year} {dataset}',
        'retries': 1, 'skipbadfiles': True, 'xrootdtimeout':120,
        'tailtimeout': 600,
    }
    pre_args = {
        'desc': f'Preprocessing {args.baseprocessor} {args.year} {dataset}',
        'tailtimeout': 600,
    }
    if args.dask:
        executor = processor.dask_executor
        executor_args['client'] = client
    elif args.parsl:
        executor = processor.parsl_executor
    else:
        if args.workers<=1:
            executor = processor.iterative_executor
        else:
            executor = processor.futures_executor
        executor_args['workers'] = args.workers

    accumulator = processor.run_uproot_job(
        fileset,
        treename='Events',
        processor_instance=processor_instance,
        executor=executor,
        executor_args=executor_args,
        pre_args=pre_args,
        chunksize=300000, # 200000 good for condor 1000 MB, request 2000 MB/core
    )
    
    save(accumulator, args.output)

    if args.parsl:
        parsl.dfk().cleanup()
        parsl.clear()
