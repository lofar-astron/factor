"""
Operation: facet_setup
Implements the initial setup for a facet
"""

from factor.lib.operation import operation
import factor.actions as a

class facet_setup(operation):

    def run(self):
        
        import logging
        log = logging.getLogger(self.name)

        mss = self.parset['mss']

        # add the central source to the visibilities
        log.info('Adding calibrator to vis...')
        #actions = [ a.add_cal(ms, direction) for ms in mss ]
        #self.s.run_action_parallel( actions )
        #result = self.s.get_result()

        # phase shift in the DD-calibrator direction
        log.info('Phase shifting...')
        actions = [ a.phase_shifter.phase_shifter(ms, self.direction) for ms in mss ]
        self.s.run_action_parallel( actions ) 
        result = self.s.get_result()
        
        # averaging to...
        log.info('Averaging...')
        #actions = [ actions.avg(ms, time='1', freq='1') for ms in mss ]
        #self.s.run_action_parallel( actions ) 
        #result = self.s.get_result()
        
        # apply direction independent?
        log.info('Apply DDI...')
        #actions = [ actions.apply_die(ms, die_parmdb='') for ms in mss ]
        #self.s.run_action_parallel( actions ) 
        #result = self.s.get_result()

