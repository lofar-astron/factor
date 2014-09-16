"""
Operation: example
An example operation, to be used as a template
"""

from factor.lib.operation import operation
import factor.actions as a

class example_operation(operation):

    def run(self):

        import logging
        log = logging.getLogger(self.name)

        mss = self.parset['mss']

        # this is a typical multi-thread block
        log.info('Starting example procedure...')
        actions = [ a.action_name(ms, direction, arg1, arg2) for ms in mss ] # create multiple action objects
        self.s.run_action_parallel( actions ) # call the scheduler (it's in self.s)
        result = self.s.get_result() # collect results

        # this is a single call
        log.info('Starting example procedure...')
        self.s.run_action_parallel( a.action_name(ms, arg1, arg2) )

