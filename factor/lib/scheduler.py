from Queue import Queue
from multiprocessing import Pool

class scheduler():
    """
    Schedule jobs
    The scheduler run all actions sent to it in parallel
    """
    def __init__(self, max_threads = 1):
        self.q = Queue()

    def put(self, command):
        self.q.put(command)

    def run_action_parallel(self, action_list):
        """
        List of actions to run in parallel
        """

    def run_action(self, command):
        """
        run a single action
        """
