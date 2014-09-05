"""
Operation: tessellation
Implements the tassellation of the whole field given a set of calibrators.
"""

from factor.lib.operation import *
import logging

class tessellation(operation):

    def __init__(self, parset):
        operation(parset)
        self.name = 'Tessellation'

    def run(self):
        pass
