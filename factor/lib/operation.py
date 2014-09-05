"""
General operation library, contains the master class for operations
"""
import logging

class operation():
    """
    General class for operations. All operations should be in a separate module.
    Every module must have a class called in the same way of the module which
    inherits from this class.
    """

    def __init__(self, parset):
        self.name = None
        self.parset = parset
    
    def setup(self):
        logging.info('<-- Operation %s started.' % self.name)
        
    def run(self):
        raise(NotImplementedError)

    def finalize(self):
        logging.info('--> Operation %s terminated.' % self.name)

class operation_facet(operation):
    """
    Extend the operation class for facet-specific operations
    """

    def __init__(self, parset, direction):
        operation(parset)
        self.direction = direction


