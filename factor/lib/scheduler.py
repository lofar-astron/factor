import logging
from Queue import Queue
from threading import Thread
from factor.lib.context import Timer


class Scheduler(object):
    """
    Schedule jobs
    The scheduler runs all jobs sent to it in parallel
    """
    def __init__(self, max_threads=1, name=''):
        """
        max_threads: limit the number of parallel process to this number
        name: name of the operation requesting a scheduler
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


    def run_parallel(self, op_list):
        """
        Runs a list of operations in parallel
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


    def get_result(self):
        """
        """
        pass
