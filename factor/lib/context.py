"""
Definition of context managers (with statements) used for operations
"""

class op_timer():
    """
    context manager used to time the operations
    """

    def __enter__(self):
        import time
        self.start = time.clock()

    def __exit__(self, type, value, tb):
        import logging, time
        if type is not None:
            pass
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
        #obj = getattr(operations[name], name)
        obj = getattr(getattr(op, self.name), self.name)
        if self.direction != None: return obj(self.parset, self.direction)
        else: return obj(self.parset)

    def __exit__(self, type, value, tb):
        import logging
        if type is not None:
            logging.error('Problems running operation: %s' % (self.name))
            raise(type)
        logging.info('Operation %s, terminated successfully.' % (self.name))



