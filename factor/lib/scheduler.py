from Queue import Queue
from multiprocessing import Pool

class scheduler():
    """
    Schedule jobs
    The scheduler run all actions sent to it in parallel
    """
    def __init__(self, max_threads = 1):
        self.Pool = Pool(max_threads)

    def run_action_parallel(self, action_list):
        """
        List of actions to run in parallel
        """
        for action in action_list:
            self.Pool(target=action.run())


    def run_action(self, command):
        """
        run a single action
        """
