import logging
from Queue import Queue
from threading import Thread

class scheduler( object ):
    """
    Schedule jobs
    The scheduler run all actions sent to it in parallel
    """
    def __init__(self, max_threads = 1):
        self.max_threads = max_threads

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
        logging.debug("Run %i action(s) on %i cpu(s)." % (len(action_list), self.max_threads))
        for action in action_list:
            self.q.put_nowait(action)
        for _ in self.threads: self.q.put(None) # signal no more commands
        for t in self.threads: t.join() # wait for completion

    def get_result(self):
        """
        """
        pass
