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


    def result_callback(self, result):
        """
        Callback function for apply_async result
        """
        op_name, direction_name, status = result

        if status == 0:
            log.info('--> Operation {0} completed (direction: '
                '{1})'.format(op_name, direction_name))

            # Identify the current operation from the direction name
            this_op = None
            for op in self.operation_list:
                if op.direction.name == direction_name:
                    this_op = op
                    break

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
        self.operation_list = operation_list

        # Run the operation(s)
        if len(operation_list) > 0:
            with Timer(log, 'operation'):
                pool = multiprocessing.Pool(processes=self.max_procs)
                for op in operation_list:
                    if not self.dry_run and not op.check_completed():
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
