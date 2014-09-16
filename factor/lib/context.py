"""
Definition of context managers (with statements) used for operations
"""
import time
import logging

class op_timer():
    """
    context manager used to time the operations
    """

    def __init__(self, operation):
        self.name = operation.name

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, type, value, tb):

        if type is not None:
            raise type, value, tb

        elapsed = (time.time() - self.start)
        log = logging.getLogger(self.name)
        log.debug('Time for operation: %i sec' % (elapsed))

class op_init():
    """
    context manager used to initialize an operations
    """

    def __init__(self, name, parset, direction = None):

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
