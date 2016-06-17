"""
Module defining the operation scheduler class
"""
import logging
import multiprocessing
import os
import sys
import imp
from factor.lib.context import Timer

log = logging.getLogger('factor:scheduler')


def call_generic_pipeline(op_name, direction_name, parset, config, logbasename,
    genericpipeline_executable):
    """
    Creates a GenericPipeline object and runs the pipeline

    Parameters
    ----------
    op_name : str
        Name of operation
    direction_name : str
        Name of direction
    parset : str
        Name of pipeline parset file
    config : str
        Name of pipeline config file
    logbasename : str
        Log file base name
    genericpipeline_executable : str
        Path to genericpipeline.py executable

    """
    from lofarpipe.support.pipelinelogging import getSearchingLogger
    from factor.lib.context import RedirectStdStreams
    genericpipeline_path = os.path.dirname(genericpipeline_executable)
    loader = imp.load_source('loader', os.path.join(genericpipeline_path,
        'loader.py'))
    gp = imp.load_source('gp', genericpipeline_executable)

    # Initalize pipeline object
    pipeline = gp.GenericPipeline()

    # Add needed attr/methods
    pipeline.name = '{0}_{1}'.format(op_name, direction_name)
    pipeline.logger = getSearchingLogger(pipeline.name)
    pipeline.inputs['args'] = [parset]
    pipeline.inputs['config'] = config
    pipeline.inputs['job_name'] = direction_name

    # Set pipeline logging to DEBUG level
    logging.root.setLevel(logging.DEBUG)
    pipeline.logger.setLevel(logging.DEBUG)
    for handler in pipeline.logger.handlers:
        handler.setLevel(logging.DEBUG)

    # Run the pipeline, redirecting screen output to log files
    log.info('<-- Operation {0} started (direction: {1})'.format(op_name,
        direction_name))
    with open("{0}.out.log".format(logbasename), "wb") as out, \
        open("{0}.err.log".format(logbasename), "wb") as err:
        with RedirectStdStreams(stdout=out, stderr=err):
            status = pipeline.run(pipeline.name)

    return (op_name, direction_name, status)


class Scheduler(object):
    """
    The scheduler runs all jobs sent to it in parallel

    Parameters
    ----------
    genericpipeline_executable : str
        Path to genericpipeline.py executable
    max_procs : int, optional
        Limit the number of parallel processes to this number
    name : str, optional
        Name of the scheduler
    dry_run : bool, optional
        If True, the pipelines are not run but all parsets and config files
        are made as normal

    """
    def __init__(self, genericpipeline_executable, max_procs=1, name='scheduler',
        dry_run=False):
        self.genericpipeline_executable = genericpipeline_executable
        self.max_procs = max_procs
        self.name = name
        self.dry_run = dry_run
        self.success = True


    def allocate_resources(self):
        """
        Divide up nodes and cpus among the operations to be run in parallel

        The following values are set:
        nimg_per_node : number of imagers per node
        max_cpus_per_node : maximum number of cores that the pipeline
            should use
        max_cpus_per_img : maximum number of threads per multi-process
            imager call
        max_io_proc_per_node : maximum number of IO-intensive processes
        max_cpus_per_chunk : number of threads in NDPPP calls that are run "ntimes" times
        max_cpus_per_band : number of threads in NDPPP calls that are run "nfiles" times
        max_percent_memory : percentage of memory to use in imagers that are only run once (e.g. imaging the full facet)
        max_percent_memory_per_img : percentage of memory to use in multi-process imager calls

        """
        node_list = self.operation_list[0].node_list
        all_directions = [op.direction for op in self.operation_list]

        ndir_simul = self.max_procs
        for i in range(int(np.ceil(len(all_directions)/float(ndir_simul)))):
            directions = all_directions[i*ndir_simul:(i+1)*ndir_simul]

            if len(directions) >= len(node_list):
                for i in range(len(directions)-len(node_list)):
                    node_list.append(node_list[i])
                hosts = [[n] for n in node_list]
            else:
                parts = len(directions)
                hosts = [node_list[i*len(node_list)//parts:
                    (i+1)*len(node_list)//parts] for i in range(parts)]

            # Find duplicates and divide up available nodes and cores
            h_flat = []
            for h in hosts:
                h_flat.extend(h)
            c = Counter(h_flat)
            for d, h in zip(directions, hosts):
                d.hosts = h
                if len(h) == 1:
                    ndir_per_node = min(ndir_per_node, c[h[0]])
                else:
                    ndir_per_node = 1
                d.nimg_per_node = nimg_per_node
                d.max_cpus_per_node =  max(1, int(round(ncpu_max / float(ndir_per_node))))
                d.max_cpus_per_img =  max(1, int(round(ncpu_max / float(nimg_per_node))))
                d.max_io_proc_per_node = int(np.ceil(np.sqrt(ncpu_max)))
                nchunks_per_node = max(1, int(round(float(len(self.bands[0].nfiles)) / len(d.hosts))))
                d.max_cpus_per_chunk = int(round(d.max_cpus_per_node / nchunks_per_node))
                d.max_cpus_per_band = max(1, int(round(d.max_cpus_per_node *
                    len(d.hosts) / float(nbands))))
                d.max_percent_memory = fmem_max / float(ndir_per_node) * 100.0
                d.max_percent_memory_per_img = fmem_max / float(ndir_per_node) / float(nimg_per_node) * 100.0
                d.save_state()


    def result_callback(self, result):
        """
        Callback function for apply_async result
        """
        op_name, direction_name, status = result

        # Identify the current operation from the direction name
        try:
            this_op_indx = [op.direction.name for op in self.operation_list].index(direction_name)
            this_op =  self.operation_list[this_op_indx]
        except ValueError:
            this_op = None

        if status == 0:
            log.info('--> Operation {0} completed (direction: '
                '{1})'.format(op_name, direction_name))

            # Finalize the operation
            if this_op is not None:
                this_op.finalize()
                this_op.set_completed()
            else:
                log.error('Operation {0} (direction: {1}) not in list of active '
                    'operations'.format(op_name, direction_name))
                self.success = False
        else:
            log.error('Operation {0} failed due to an error (direction: '
                '{1})'.format(op_name, direction_name))
            self.success = False

        # Give the completed op's nodes to the next one in line (if any)
        if this_op is not None and len(self.queued_ops) > 0:
            next_op = self.queued_ops.pop(0)
            next_op.direction.hosts = this_op.direction.hosts[:]
            next_op.setup()


    def run(self, operation_list):
        """
        Runs a list of operations in parallel

        Parameters
        ----------
        operation_list : Operation instance or list of Operation instances
            List of operations to process

        """
        if type(operation_list) != list:
            operation_list = [operation_list]
        self.operation_list = [op for op in operation_list if not op.check_completed()]
        if len(self.operation_list) == 0:
            return

        # Run the operation(s)
        self.num_complete = 0
        self.allocate_resources()
        with Timer(log, 'operation'):
            pool = multiprocessing.Pool(processes=self.max_procs)
            self.queued_ops = self.operation_list[self.max_procs:]
            for op in self.operation_list:
                if not self.dry_run:
                    # Only run incomplete operations (and only if this is
                    # not a dry run)
                    op.setup()
                    op.set_started()
                    pool.apply_async(call_generic_pipeline, (op.name,
                        op.direction.name, op.pipeline_parset_file,
                        op.pipeline_config_file, op.logbasename,
                        self.genericpipeline_executable),
                        callback=self.result_callback)
                else:
                    # For completed operations or dry runs, run finalize() to
                    # be sure that all attributes are set properly
                    op.finalize()
                    if not self.dry_run:
                        op.set_completed()
            pool.close()
            pool.join()

        if not self.success:
            log.error('One or more operations failed due to an error. Exiting...')
            sys.exit(1)
