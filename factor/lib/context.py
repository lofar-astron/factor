"""
Definition of context managers (with statements) used for operations
"""
import logging, time

class op_timer():
    """
    context manager used to time the operations
    """

    def __enter__(self):
        self.start = time.clock()

    def __exit__(self, type, value, tb):

        if type is not None:
            return 1

        elapsed = (time.clock() - self.start)
        logging.debug('Time for operation: %i sec' % (elapsed))

class op_init():
    """
    context manager used to initialize an operations
    """

    def __init__(self,name, parset, direction = None):

        self.name = name
        self.parset = parset
        self.direction = direction

    def __enter__(self):

        import factor.operations as op
        # find the module and create the object from the name
        op_class = getattr(getattr(op, self.name), self.name)
        if self.direction != None: self.op_obj = op_class(self.parset, self.direction)
        else: self.op_obj = op_class(self.parset)
        self.op_obj.setup()
        return  self.op_obj

    def __exit__(self, type, value, tb):

        # catch any exceptions from this operation
        if type is not None:
            logging.error('Problems running operation: %s' % (self.name))
            raise(value)
        self.op_obj.finalize()
