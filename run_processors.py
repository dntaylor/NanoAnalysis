#!/usr/bin/env python
import os
import json
import argparse
import glob
import time
import concurrent.futures
from functools import partial
import logging

logging.basicConfig(filename='_run_processors.log', level=logging.INFO, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
rootLogger = logging.getLogger()
logging.captureWarnings(True)

import uproot
import numpy as np
import cloudpickle
import lz4.frame as lz4f
from coffea import hist, processor
from coffea.util import load, save


from parsl_config import parsl_condor_config, parsl_local_config

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run coffea file')
    parser.add_argument('-j', '--workers', type=int, default=1, help='Number of workers to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
    scheduler = parser.add_mutually_exclusive_group()
    scheduler.add_argument('--parsl', action='store_true', help='Use parsl to distribute')
    parser.add_argument('--condor', action='store_true', help='Use distributed, but with condor')
    args = parser.parse_args()

    # load the distributed computing interface
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
    
    jobs = {}

    index = 0
    for baseprocessor in baseprocessors:
        for year in years:

            processorpath = f'processors/{baseprocessor}_{year}.coffea'
            if os.path.exists(processorpath):
                processor_instance = load(processorpath)
            else:
                rootLogger.warning(f'Cannot understand {processorpath}.')

            
            filesetpaths = glob.glob(f'filesets/{year}/*.json')
            for filesetpath in filesetpaths:
                # TODO: only run datasets requested by the given processor
                dataset = os.path.splitext(os.path.basename(filesetpath))[0]
                if dataset in ['all','mc','data']: continue

                output = f'hists/{baseprocessor}/{year}/{dataset}.coffea'
    
                rootLogger.info(f'Will process: {baseprocessor} {year} {dataset}')
    
                executor_args = {
                    'savemetrics': True, 'flatten':True, 
                    'desc': f'Processing {baseprocessor} {year} {dataset}', #'position': -1*(2*index+1),
                    'retries': 1, 'skipbadfiles': True, 'xrootdtimeout':120,
                }
                pre_args = {
                    'desc': f'Preprocessing {baseprocessor} {year} {dataset}',
                    #'position': -1*(2*index+2),
                }
    
                with open(filesetpath) as f:
                    fileset = json.load(f)
                    for dataset in fileset:
                        fileset[dataset] = [redirector+x if x.startswith('/store') else x for x in fileset[dataset]]

                pi_to_send = lz4f.compress(cloudpickle.dumps(processor_instance), compression_level=1)
    

                if index>4: break
                jobs[output] = {'args': [fileset], 'kwargs': {'pi_to_send': pi_to_send, 'executor_args': executor_args, 'pre_args': pre_args,}, }
                index += 1

    def parallel_run_uproot_job(fileset, **kwargs):
        pi_to_send = kwargs.pop('pi_to_send', None)
        executor_args = kwargs.pop('executor_args', {})
        pre_args = kwargs.pop('pre_args', {})

        from coffea import processor
        import cloudpickle
        import lz4.frame as lz4f

        processor_instance = cloudpickle.loads(lz4f.decompress(pi_to_send))

        executor = processor.parsl_executor
    
        return processor.run_uproot_job(fileset,
            treename='Events',
            processor_instance=processor_instance,
            executor=executor,
            executor_args=executor_args,
            pre_args=pre_args,
            chunksize=200000, # 200000 good for condor 1000 MB, request 2000 MB/core
        )


    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:
        futures = {}
        for output in jobs:
            function = partial(parallel_run_uproot_job, **jobs[output]['kwargs'])
            futures[output] = pool.submit(function, *jobs[output]['args'])

        while futures:
            finished = set(output for output,future in futures.items() if future.done())
            while finished:
                output = finished.pop()
                accumulator = futures.pop(output).result()
                os.makedirs(os.path.dirname(output),exist_ok=True)
                save(accumulator,output)
            time.sleep(0.5)


    parsl.dfk().cleanup()
    parsl.clear()
