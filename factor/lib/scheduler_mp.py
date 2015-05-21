"""
Module defining the operation scheduler class for multiproccessing
"""
import logging
import multiprocessing
import os
from factor.lib.context import Timer


def call_generic_pipeline(executable, parset, config):
    """
    Creates a GenericPipeline object and runs the pipeline
    """
    from genericpipeline.bin import genericpipeline as gp
    from lofarpipe.support.pipelinelogging import getSearchingLogger
    from factor.lib.context import RedirectStdStreams
    import sys

    # Initalize pipeline object
    pipeline = gp.GenericPipeline()

    # Add needed attr/methods
    pipeline.name = os.path.splitext(os.path.basename(executable))[0]
    pipeline.logger = getSearchingLogger(pipeline.name)
    pipeline.inputs['args'] = [parset]
    pipeline.inputs['config'] = config

    # Run the pipeline
    status = pipeline.run(pipeline.name)

    return status


class Scheduler(object):
    """
    The scheduler runs all jobs sent to it in parallel
    """
    def __init__(self, max_procs=1, name='', op_parset=None):
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

        """
        self.max_procs = max_procs
        self.name = name
        self.op_parset = op_parset
        self.log = logging.getLogger('factor.{0}'.format(name))


    def run(self, action_list):
        """
        Runs a list of actions in parallel

        If a distributed file system is used, all non-localhost files are
        synchronized with the distributed versions

        Parameters
        ----------
        action_list : Action or list of Actions
            Single action or list of actions to run

        Returns
        -------
        results : Datamap or list of Datamaps
            Single result or list of results returned by the actions

        """
        if type(action_list) != list:
            single = True
            action_list = [action_list]
        else:
            single = False

        # Sync the local node to the remote nodes
        if self.op_parset['cluster_specific']['distribute']:
            host_list = self.op_parset['cluster_specific']['node_list']
            working_dir = self.op_parset['dir_working']
            for host in host_list:
                if host != self.hostname:
                    os.system('rsync -az --delete {0}/ {1}:{0}'.format(
                        working_dir, host))

        # Run the action(s)
        self.log.debug('Running up to {0} actions in parallel'.format(self.max_procs))
        with Timer(self.log, 'action'):
            pool = multiprocessing.Pool(processes=self.max_procs)
            for act in action_list:
                pool.apply_async(call_generic_pipeline, (act.pipeline_executable,
                    act.pipeline_parset_file, act.pipeline_config_file))
            pool.close()
            pool.join()

        # Sync the remote nodes to the local node
        if self.op_parset['cluster_specific']['distribute']:
            host_list = self.op_parset['cluster_specific']['node_list']
            working_dir = self.op_parset['dir_working']
            for host in host_list:
                if host != self.hostname:
                    os.system('rsync -az --delete {1}:{0}/ {0}'.format(
                        working_dir, host))

        # Get results
        results = []
        for action in action_list:
            results.append(action.get_results())
        if single:
            results = results[0]

        return results
