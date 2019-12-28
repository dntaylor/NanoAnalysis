#!/usr/bin/env python
import os
import json
import argparse
import glob
import logging

logging.basicConfig(filename='run_processors.log', level=logging.INFO)
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
    parser.add_argument('-j', '--workers', type=int, default=1, help='Number of workers to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
    scheduler = parser.add_mutually_exclusive_group()
    scheduler.add_argument('--dask', action='store_true', help='Use dask to distribute')
    scheduler.add_argument('--parsl', action='store_true', help='Use parsl to distribute')
    parser.add_argument('--condor', action='store_true', help='Use distributed, but with condor')
    args = parser.parse_args()

    # load the distributed computing interface
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

    # setup jobs
    baseprocessors = ['hzzProcessor']
    years = ['2016','2017','2018']

    # TODO: will need to add feature that checks if a job needs to be updated.
    # idea: if processor/fileset timestep newer than the output hists timestep, rerun
    # will additionally want to load, stop the condor jobs if no new stuff needs to be done

    for baseprocessor in baseprocessors:
        for year in years:
            
            processorpath = f'processors/{baseprocessor}_{year}.coffea'

            filesetpaths = glob.glob(f'filesets/{year}/*.json')
            for filesetpath in filesetpaths:
                # TODO: only run datasets requested by the given processor
                dataset = os.path.splitext(os.path.basename(filesetpath))[0]
                if dataset in ['all','mc','data']: continue

                output = f'hists/{baseprocessor}/{year}/{dataset}.coffea'
                python_mkdir(os.path.dirname(output))


                with open(filesetpath) as f:
                    fileset = json.load(f)
                    for dataset in fileset:
                        fileset[dataset] = [redirector+x if x.startswith('/store') else x for x in fileset[dataset]]

                rootLogger.info(f'Will process: {baseprocessor} {year} {dataset}')

                if os.path.exists(processorpath):
                    processor_instance = load(processorpath)
                else:
                    rootLogger.warning(f'Cannot understand {processorpath}.')

                try:
                    executor_args = {
                        'savemetrics': True, 'flatten':True, 
                        'desc': f'Processing {baseprocessor} {year} {dataset}',
                        'retries': 1, 'skipbadfiles': True, 'xrootdtimeout':120,
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

                    pre_args = {
                        'desc': f'Preprocessing {baseprocessor} {year} {dataset}',
                    }

                    accumulator = processor.run_uproot_job(
                        fileset,
                        treename='Events',
                        processor_instance=processor_instance,
                        executor=executor,
                        executor_args=executor_args,
                        pre_args=pre_args,
                        chunksize=300000, # 200000 good for condor 1000 MB, request 2000 MB/core
                    )
                    
                    save(accumulator, output)

                # TODO: special handling, for example, preventing all future instances of a given processor...
                # for now, just raise the error.
                # eventually would like this all to be done continuously and never kill the parent
                except Exception as e:
                    raise e
