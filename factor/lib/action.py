"""
General action library, contains the master class for actions
Commands contained in an action are run in sequential mode.
"""

import logging

class action( object ):
    """
    Generic action class
    """
    def __init__(self, name = None):
        self.name = name

    def run(self):
        raise NotImplementedError


