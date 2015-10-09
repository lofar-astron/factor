"""
Definition of context managers
"""
import time
import datetime
import logging

class Timer(object):
    """
    Context manager used to time operations
    """
    def __init__(self, log=None, type='operation'):
        """
        Create object

        Parameters
        ----------
        log : logging instance
            The logging instance to use. If None, root is used
        type : str, optional
            Type of operation
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

        elapsed = time.time() - self.start
        self.log.debug('Time for {0}: {1:0>8}'.format(self.type,
            datetime.timedelta(seconds=elapsed)))


class RedirectStdStreams(object):
    """
    Context manager used to redirect streams
    """
    def __init__(self, stdout=None, stderr=None):
        import sys
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr

    def __enter__(self):
        import sys
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        import sys
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr
