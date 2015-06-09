"""
Module defining the operation scheduler class for multiproccessing
"""
import logging
import multiprocessing
import os
import sys
from factor.lib.context import Timer
import factor._logging


def call_generic_pipeline(op_name, direction_name, parset, config, logbasename):
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

    """
    from genericpipeline.bin import genericpipeline as gp
    from lofarpipe.support.pipelinelogging import getSearchingLogger
    from factor.lib.context import RedirectStdStreams
    import sys

    # Initalize pipeline object
    pipeline = gp.GenericPipeline()

    # Add needed attr/methods
    pipeline.name = '{0}_{1}'.format(op_name, direction_name)
    pipeline.logger = getSearchingLogger(pipeline.name)
    pipeline.inputs['args'] = [parset]
    pipeline.inputs['config'] = config
    pipeline.inputs['job_name'] = direction_name

    # Run the pipeline, redirecting screen output to log files
    with open("{0}.out.log".format(logbasename), "wb") as out, \
        open("{0}.err.log".format(logbasename), "wb") as err:
        with RedirectStdStreams(stdout=out, stderr=err):
            status = pipeline.run(pipeline.name)

    return (op_name, direction_name, status)


class Scheduler(object):
    """
    The scheduler runs all jobs sent to it in parallel
    """
    def __init__(self, max_procs=1, name='scheduler', op_parset=None, dry_run=False):
        """
        Create Scheduler object

        Parameters
        ----------
        max_procs : int, optional
            Limit the number of parallel processes to this number
        name : str, optional
            Name of the scheduler
        op_parset : dict, optional
            Dict of operation parameters
        dry_run : bool, optional
            If True, the pipelines are not run, but all parsets and config files
            are made as normal

        """
        self.max_procs = max_procs
        self.name = name
        self.op_parset = op_parset
        self.dry_run = dry_run
        self.log = logging.getLogger('factor.{0}'.format(name))
        self.success = True


    def result_callback(self, result):
        """
        Callback function for apply_async result
        """
        op_name, direction_name, status = result

        if status == 0:
            self.log.info('--> Operation {0} completed successfully (direction: '
                '{1})'.format(op_name, direction_name))
        else:
            self.log.error('Operation {0} failed (direction: {1})'.format(op_name,
                direction_name))
            self.success = False


    def run(self, operation_list):
        """
        Runs a list of operations in parallel

        Each operation is checked for whether it was already successfully
        completed, and, if so, it will not be scheduled. However, the finalize()
        method will be called for all operations, whether or not they were
        scheduled. This is to ensure proper handling of direction object
        attributes.

        Parameters
        ----------
        operation_list : Operation instance or list of Operation instances
            List of operations to process

        """
        if type(operation_list) != list:
            operation_list = [operation_list]

        # Check state of each operation
        operations_to_run = [op for op in operation_list if not op.check_completed()]

        # Set up the operation(s)
        for op in operations_to_run:
             op.setup()

        # Run the operation(s)
        if not self.dry_run and len(operations_to_run) > 0:
            with Timer(self.log, 'operation'):
                pool = multiprocessing.Pool(processes=self.max_procs)
                for op in operations_to_run:
                    self.log.info('<-- Operation {0} started (direction: {1})'.
                        format(op.name, op.direction.name))
                    pool.apply_async(call_generic_pipeline, (op.name,
                    	op.direction.name, op.pipeline_parset_file,
                    	op.pipeline_config_file, op.logbasename),
                    	callback=self.result_callback)
                pool.close()
                pool.join()

            if not self.success:
                self.log.error('One or more operations failed. Exiting...')
                sys.exit(1)

        # Finalize the operations
        for op in operation_list:
            op.finalize()
            if not self.dry_run:
                op.set_completed()
