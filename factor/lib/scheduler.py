"""
Module defining the operation scheduler class
"""
import logging
import multiprocessing
import os
import sys
import imp
from factor.lib.context import Timer
import factor._logging


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
    gp = imp.load_source('gp', os.path.join(genericpipeline_path,
        'genericpipeline.py'))

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
    def __init__(self, genericpipeline_executable, max_procs=1, name='scheduler',
        dry_run=False):
        """
        Create Scheduler object

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
        self.genericpipeline_executable = genericpipeline_executable
        self.max_procs = max_procs
        self.name = name
        self.dry_run = dry_run
        self.log = logging.getLogger('factor.{0}'.format(name))
        self.success = True


    def result_callback(self, result):
        """
        Callback function for apply_async result
        """
        op_name, direction_name, status = result

        if status == 0:
            self.log.info('--> Operation {0} completed (direction: '
                '{1})'.format(op_name, direction_name))
        else:
            self.log.error('Operation {0} failed due to an error (direction: '
                '{1})'.format(op_name, direction_name))
            self.success = False


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

        # Set up the operation(s)
        for op in operation_list:
             op.setup()

        # Run the operation(s)
        if not self.dry_run and len(operation_list) > 0:
            with Timer(self.log, 'operation'):
                pool = multiprocessing.Pool(processes=self.max_procs)
                for op in operation_list:
                    self.log.info('<-- Operation {0} started (direction: {1})'.
                        format(op.name, op.direction.name))
                    pool.apply_async(call_generic_pipeline, (op.name,
                    	op.direction.name, op.pipeline_parset_file,
                    	op.pipeline_config_file, op.logbasename,
                    	self.genericpipeline_executable),
                    	callback=self.result_callback)
                pool.close()
                pool.join()

            if not self.success:
                self.log.error('One or more operations failed due to an error. Exiting...')
                sys.exit(1)

        # Finalize the operations
        for op in operation_list:
            op.finalize()
            if not self.dry_run:
                op.set_completed()
