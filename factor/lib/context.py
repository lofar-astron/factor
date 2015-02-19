"""
Definition of context managers (with statements) used for operations
"""
import time
import logging

class Timer():
    """
    context manager used to time the operations
    """

    def __init__(self, log=''):
        """
        log: is a logging istance to print the correct log format
        if nothing is passed, root is used
        """
        if log == '': self.log = logging
        else: self.log = log

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, type, value, tb):

        if type is not None:
            raise type, value, tb

        elapsed = (time.time() - self.start)
        self.log.debug('Time for operation: %i sec' % (elapsed))

