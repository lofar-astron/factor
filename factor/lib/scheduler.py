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
    def __init__(self, max_threads=1, name=''):
        """
        Create Scheduler object

        Parameters
        ----------
        max_threads : int, optional
            Limit the number of parallel process to this number
        name : str, optional
            Name of the scheduler

        """
        self.max_threads = max_threads
        self.name = name
        self.log = logging.getLogger(name)

        def worker(queue):
            for cmd in iter(queue.get, None):
                cmd.run()

        self.q = Queue()
        self.threads = [Thread(target=worker, args=(self.q,)) for _ in range(self.max_threads)]
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
        self.threads = [Thread(target=worker, args=(self.q,)) for _ in range(self.max_threads)]
        for i, t in enumerate(self.threads): # start workers
            t.daemon = True
            t.start()


    def run(self, op_list):
        """
        Runs a list of operations in parallel

        Parameters
        ----------
        op_list : list of Operation objects
            List of operations to run

        """
        self.startup()
        if type(op_list) != list:
            op_list = [op_list]
        with Timer(self.log):
            for op in op_list:
                self.q.put_nowait(op)
            for _ in self.threads:
                self.q.put(None) # signal no more commands
            for t in self.threads:
                t.join() # wait for completion
