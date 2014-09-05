"""
General operation library, contains the master class for operations
"""

class operation():
    """
    General class for operations. All operations should be in a separate module.
    Every module must have a class called in the same way of the module which
    inherits from this class.
    """
    import logging

    def __init__(self, parset):
        self.parset = parset

class operation_facet(operation):
    """
    Extend the operation class for facet-specific operations
    """

    def __init__(self, parset, direction):
        self.operation(parset)
        self.direction = direction
