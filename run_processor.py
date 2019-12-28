#!/usr/bin/env python
import os
import json
import argparse
import glob
import logging

logging.basicConfig(filename='run_processor.log', level=logging.INFO)
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
    parser.add_argument('processor', type=str, default=None, help='The analysis, a precompiled processor')
    parser.add_argument('--fileset', default='', help='Fileset to process')
    parser.add_argument('--output', default='', help='Output histogram filename')
    parser.add_argument('-j', '--workers', type=int, default=1, help='Number of workers to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
    scheduler = parser.add_mutually_exclusive_group()
    scheduler.add_argument('--dask', action='store_true', help='Use dask to distribute')
    scheduler.add_argument('--parsl', action='store_true', help='Use parsl to distribute')
    parser.add_argument('--condor', action='store_true', help='Use distributed, but with condor')
    parser.add_argument('--test', action='store_true', help='Only process a few files')
    parser.add_argument('--debug', action='store_true', help='Add debug verbosity')
    args = parser.parse_args()

    if not args.output:
        args.output = 'hists/' + os.path.splitext(os.path.basename(args.processor))[0]\
                      + '/' + os.path.splitext(os.path.basename(args.fileset))[0] + '.coffea'
    python_mkdir(os.path.dirname(args.output))

    if args.dask:
        from dask.distributed import Client, LocalCluster
        if args.condor:
            # try using dask cluster manager
            # doesn't work yet
            from dask_jobqueue.htcondor import HTCondorCluster
            cluster = HTCondorCluster(cores=1, memory="1GB", disk="2GB")
            cluster.scale(jobs=args.workers)
            client = Client(cluster)
            # manually creating cluster
            #client = Client(os.environ['DASK_SCHEDULER'])
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
            

    #redirector = 'root://cms-xrd-global.cern.ch/'
    #redirector = 'root://xrootd-cms.infn.it/'
    redirector = 'root://cmsxrootd.fnal.gov/'
    #redirector = 'root://cmsxrootd.hep.wisc.edu/'

    if args.fileset:
        with open(args.fileset) as f:
            fileset = json.load(f)
            for dataset in fileset:
                fileset[dataset] = [redirector+x if x.startswith('/store') else x for x in fileset[dataset]]
    elif args.test:
        fileset = {'DY': ['dy_2018.root'], 'HZZ': ['hzz_2018.root'], 'DoubleMuon': ['DoubleMuon_2018.root'],}

    rootLogger.info('Will process: '+' '.join(list(fileset.keys())))

    if os.path.exists(args.processor):
        processor_instance = load(args.processor)
    else:
        rootLogger.warning(f'Cannot understand {args.processor}.')

    executor_args = {'savemetrics': True, 'flatten':True, 'retries': 1, 'skipbadfiles': True, 'xrootdtimeout':120,}
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
        chunksize=300000, # 200000 good for condor 1000 MB, request 2000 MB/core
    )
    
    save(accumulator, args.output)
    if args.test: print(accumulator)
