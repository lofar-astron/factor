"""
Operation: init_subtract
Implements the tassellation of the whole field given a set of calibrators.
"""

from factor.lib.operation import *
import logging

class init_subtract(operation):

    def __init__(self, parset):
        operation(parset)
        self.name = 'Initial subtraction'

    def run(self):
        pass
