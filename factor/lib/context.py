"""
Definition of context managers (with statements) used for operations
"""

class op_timer:
    """
    context manager used to time the operations
    """
    import time

    def __enter__(self):
        self.start = time.clock()

    def __exit__(self, type, value, tb):
        if type is not None:
            pass
        elapsed = (time.clock() - self.start)
        logging.info('Time for operation: %i sec'.format(elapsed))

class op_logger:
    """
    context manager used to log the operations
    """

    def __enter__(self, operation):

        logging.info()
        return operation(param)

    def __exit__(self, type, value, tb):
        if type is not None:
            logging.error('')
        logging.info()
