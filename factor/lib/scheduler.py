import logging
from Queue import Queue
from threading import Thread
from factor.lib.context import timer

class scheduler( object ):
    """
    Schedule jobs
    The scheduler run all actions sent to it in parallel
    """
    def __init__(self, max_threads = 1, name=''):
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

    def run_action_parallel(self, action_list):
        """
        List of actions to run in parallel
        """
        self.log.debug("Run %i action(s) on %i cpu(s)." % (len(action_list), self.max_threads))
        
        with timer(action_list[0].log): # start timer
            for action in action_list:
                self.q.put_nowait(action)
            for _ in self.threads: self.q.put(None) # signal no more commands
            for t in self.threads: t.join() # wait for completion

    def get_result(self):
        """
        """
        pass
