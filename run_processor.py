#!/usr/bin/env python
import os
import getpass
import json
import argparse
import glob
import logging
from tqdm import tqdm

logging.basicConfig(filename='processor.log', level=logging.INFO)
rootLogger = logging.getLogger()
logging.captureWarnings(True)


import uproot
import numpy as np
from coffea import hist, processor
from coffea.util import load, save

from hzzProcessor import HZZProcessor

UID = os.getuid()
UNAME = getpass.getuser()

processor_map = {
    'HZZ': HZZProcessor,
}





def parsl_condor_config(args):
    from parsl.providers import CondorProvider
    from parsl.channels import LocalChannel
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.addresses import address_by_hostname

    x509_proxy = f'x509up_u{UID}'
    grid_proxy_dir = '/tmp'

    cores_per_job = 1
    mem_per_core = 2000
    mem_request = mem_per_core * cores_per_job
    init_blocks = args.workers
    min_blocks = args.workers
    max_blocks = 40*args.workers
    htex_label='coffea_parsl_condor_htex'
    log_dir = 'parsl_logs'

    worker_init = f'''
echo "Setting up environment"
tar -zxf columnar.tar.gz
source columnar/bin/activate
export PATH=columnar/bin:$PATH
export PYTHONPATH=columnar/lib/python3.6/site-packages:$PYTHONPATH
export X509_USER_PROXY={x509_proxy}
echo "Environment ready"
'''

    # requirements for T2_US_Wisconsin (HAS_CMS_HDFS forces to run a T2 node not CHTC)
    scheduler_options = f'''
transfer_output_files   = {log_dir}/{htex_label}
RequestMemory           = {mem_request}
RequestCpus             = {cores_per_job}
+RequiresCVMFS          = True
Requirements            = TARGET.HAS_CMS_HDFS && TARGET.Arch == "X86_64"
priority                = 10
'''

    transfer_input_files = ['columnar.tar.gz', os.path.join(grid_proxy_dir, x509_proxy)]

    htex = Config(
        executors=[
            HighThroughputExecutor(
                label=htex_label,
                address=address_by_hostname(),
                prefetch_capacity=0,
                cores_per_worker=1,
                max_workers=cores_per_job,
                worker_logdir_root=log_dir,
                provider=CondorProvider(
                    channel=LocalChannel(),
                    init_blocks=init_blocks,
                    min_blocks=min_blocks,
                    max_blocks=max_blocks,
                    nodes_per_block=1,
                    worker_init=worker_init,
                    transfer_input_files=transfer_input_files,
                    scheduler_options=scheduler_options
                ),
            )
        ],
        strategy='simple',
        run_dir=os.path.join(log_dir,'runinfo'),
        retries = 2, # in case of xrootd errors, wish we can restrict to only xrootd error
    )

    return htex

def parsl_local_config(args):
    from parsl.providers import LocalProvider
    from parsl.channels import LocalChannel
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor

    log_dir = 'parsl_logs'

    htex = Config(
        executors=[
            HighThroughputExecutor(
                label="coffea_parsl_default",
                cores_per_worker=1,
                max_workers=args.workers,
                worker_logdir_root=log_dir,
                provider=LocalProvider(
                    channel=LocalChannel(),
                    init_blocks=1,
                    max_blocks=1,
                ),
            )
        ],
        strategy=None,
        run_dir=os.path.join(log_dir,'runinfo'),
        retries = 2,
    ) 
    return htex

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run coffea file')
    parser.add_argument('processor', type=str, default=None, help='The analysis, either precompiled or a string in processor_map')
    parser.add_argument('year', choices=['2016','2017','2018'], default='2018', help='Data taking year (default: %(default)s)')
    parser.add_argument('--fileset', default='', help='Fileset to process')
    parser.add_argument('--output', default='hists.coffea', help='Output histogram filename (default: %(default)s)')
    parser.add_argument('-j', '--workers', type=int, default=1, help='Number of workers to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
    scheduler = parser.add_mutually_exclusive_group()
    scheduler.add_argument('--dask', action='store_true', help='Use dask to distribute')
    scheduler.add_argument('--parsl', action='store_true', help='Use parsl to distribute')
    parser.add_argument('--condor', action='store_true', help='Use distributed, but with condor')
    parser.add_argument('--test', action='store_true', help='Only process a few files')
    parser.add_argument('--debug', action='store_true', help='Add debug verbosity')
    args = parser.parse_args()

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
            htex = parsl_condor_config(args)
        else:
            htex = parsl_local_config(args)
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
    else:
        with open('data.json') as f:
            fulldatafileset = json.load(f)
        with open('mc.json') as f:
            fullmcfileset = json.load(f)

        fileset = {}
        for s in fulldatafileset[args.year]:
            fileset[s] = []
            for d in fulldatafileset[args.year][s]['datasets']:
                fileset[s] += [redirector+x for x in fulldatafileset[args.year][s]['files'][d]]
        for s in fullmcfileset[args.year]:
            fileset[s] = []
            for d in fullmcfileset[args.year][s]['datasets']:
                fileset[s] += [redirector+x for x in fullmcfileset[args.year][s]['files'][d]]

    rootLogger.info('Will process: '+' '.join(list(fileset.keys())))


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
