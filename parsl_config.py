import os
import getpass

from parsl.providers import LocalProvider, CondorProvider
from parsl.channels import LocalChannel
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_hostname

UID = os.getuid()
UNAME = getpass.getuser()


def parsl_local_config(workers=1):
    log_dir = 'parsl_logs'

    htex = Config(
        executors=[
            HighThroughputExecutor(
                label="coffea_parsl_default",
                cores_per_worker=1,
                max_workers=workers,
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
        #retries = 2,
    )
    return htex

def parsl_condor_config(workers=1):

    x509_proxy = f'x509up_u{UID}'
    grid_proxy_dir = '/tmp'

    cores_per_job = 1
    mem_per_core = 2000
    mem_request = mem_per_core * cores_per_job
    init_blocks = 1
    min_blocks = 1
    max_blocks = workers
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
            ),
            # TODO: works, but really isn't helpful since half of the tasks get shipped to the condor
            # executor and don't flock back when the local executor is empty
            # an alternative could be to preprocess locally and process on the grid
            # add a local executor so stuff starts fast
            #HighThroughputExecutor(
            #    label="coffea_parsl_default",
            #    cores_per_worker=1,
            #    max_workers=1, # TODO: multicore local?
            #    worker_logdir_root=log_dir,
            #    provider=LocalProvider(
            #        channel=LocalChannel(),
            #        init_blocks=1,
            #        max_blocks=1,
            #    ),
            #),
        ],
        strategy='simple',
        run_dir=os.path.join(log_dir,'runinfo'),
        #retries = 2,
    )

    return htex
