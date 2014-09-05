"""
Operation: facet_setup
Implements the initial setup for a facet
"""

from factor.lib.operation import *
import logging

class facet_setup(operation_facet):

    def __init__(self, parset, direction):
        operation_facet(parset, direction)
        self.name = 'Facet setup'

    def run(self):
        pass
