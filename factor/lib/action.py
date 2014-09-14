"""
General action library, contains the master class for actions
Commands contained in an action are run in sequential mode.
"""

import logging
import scheduler

class action():
    """
    Generic action class
    """
    def __init__(self):
        self.scheduler = scheduler.scheduler()

    def run(self):
        raise NotImplementedError


