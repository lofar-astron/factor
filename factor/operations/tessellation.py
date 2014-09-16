"""
Operation: tessellation
Implements the tassellation of the whole field given a set of calibrators.
"""

from factor.lib.operation import operation
import factor.actions as a

class tessellation(operation):

    def run(self):

        import logging
        log = logging.getLogger(self.name)


