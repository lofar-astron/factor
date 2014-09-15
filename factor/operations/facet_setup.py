"""
Operation: facet_setup
Implements the initial setup for a facet
"""

import logging
from factor.lib.operation import operation
import factor.actions as a

class facet_setup(operation):

    def run(self):

        mss = self.parset['mss']

        # add the central source to the visibilities
        actions = [ a.add_cal(ms, direction) for ms in mss ]
        self.s.run_action_parallel( actions )
        result = self.s.get_result()

        # phase shift in the DD-calibrator direction
        actions = [ a.phase_shifter.phase_shifter(ms, self.direction) for ms in mss ]
        self.s.run_action_parallel( actions ) 
        result = self.s.get_result()
        
        # averaging to...
        #actions = [ actions.avg(ms, time='1', freq='1') for ms in mss ]
        #self.s.run_action_parallel( actions ) 
        #result = self.s.get_result()
        
        # apply direction independent?
        #actions = [ actions.apply_die(ms, die_parmdb='') for ms in mss ]
        #self.s.run_action_parallel( actions ) 
        #result = self.s.get_result()

