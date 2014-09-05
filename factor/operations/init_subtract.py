"""
Operation: init_subtract
Implements the tassellation of the whole field given a set of calibrators.
"""

from factor.lib.operations import *
import logging

class init_subtract(operation):

    def __init__(self, parset):
        operation(parset)
        logging.info('Operation init_subtract started!')
