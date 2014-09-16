"""
Definition of context managers (with statements) used for operations
"""
import time
import logging

class timer():
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

class op_init():
    """
    context manager used to initialize an operations
    """

    def __init__(self, name, parset, direction = None):
        """
        name: name of the operation
        parset: parset dict
        direction: optional direction object
        """

        self.name = name
        self.parset = parset
        self.direction = direction

    def __enter__(self):

        import factor.operations as op
        # find the module and create the object from the name
        op_class = getattr(getattr(op, self.name), self.name)
        if self.direction != None: self.op_obj = op_class(self.parset, self.name, self.direction)
        else: self.op_obj = op_class(self.parset, self.name)
        self.op_obj.setup()
        return  self.op_obj

    def __exit__(self, type, value, tb):

        # catch any exceptions from this operation
        if type is not None:
            logging.error('Problems running operation: %s' % (self.name))
            raise type, value, tb
        self.op_obj.finalize()
