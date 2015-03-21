"""
Module defining the operation scheduler class
"""
import logging
from Queue import Queue
from threading import Thread
from factor.lib.context import Timer


class Scheduler(object):
    """
    The scheduler runs all jobs sent to it in parallel
    """
    def __init__(self, max_threads=1, name='', op_parset=None):
        """
        Create Scheduler object

        Parameters
        ----------
        max_threads : int, optional
            Limit the number of parallel process to this number
        name : str, optional
            Name of the scheduler
        op_parset : dict, optional
            Dict of operation parameters

        """
        self.max_threads = max_threads
        self.name = name
        self.op_parset = op_parset
        self.log = logging.getLogger(name)

        def worker(queue):
            for cmd in iter(queue.get, None):
                cmd.run()

        self.q = Queue()
        self.threads = [Thread(target=worker, args=(self.q,)) for _ in
            range(self.max_threads)]
        for i, t in enumerate(self.threads): # start workers
            t.daemon = True
            t.start()


    def startup(self):
        """
        Starts up scheduler
        """
        def worker(queue):
            for cmd in iter(queue.get, None):
                cmd.run()
        self.q = Queue()
        self.threads = [Thread(target=worker, args=(self.q,)) for _ in
            range(self.max_threads)]
        for i, t in enumerate(self.threads): # start workers
            t.daemon = True
            t.start()


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
        self.startup()
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
        with Timer(self.log, 'action'):
            for act in action_list:
                self.q.put_nowait(act)
            for _ in self.threads:
                self.q.put(None) # signal no more commands
            for t in self.threads:
                t.join() # wait for completion

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
