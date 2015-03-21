"""
Definition of context managers (with statements) used for actions and operations
"""
import time
import logging

class Timer():
    """
    Context manager used to time actions and operations
    """
    def __init__(self, log=None, type='action'):
        """
        Create object

        Parameters
        ----------
        log : logging instance
            The logging instance to use. If None, root is used
        """
        if log is None:
            self.log = logging
        else:
            self.log = log
        self.type = type


    def __enter__(self):
        self.start = time.time()


    def __exit__(self, type, value, tb):
        if type is not None:
            raise type, value, tb

        elapsed = (time.time() - self.start)
        self.log.debug('Time for {0}: {1} sec'.format(self.type, int(elapsed)))

