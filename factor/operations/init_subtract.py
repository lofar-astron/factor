"""
Operation: init_subtract
Implements the initial creation  initial sky model and empty MS by imaging at high resolution  and low resolution and subtracting these from the datasets
"""

from factor.lib.operation import *
import logging

class init_subtract(operation):

    def __init__(self, parset):
        operation(parset)
        self.name = 'Initial subtraction'

    def run(self):
        pass
