"""
Operation: facet_setup
Implements the initial setup for a facet
"""

import logging
from factor.lib.operation import *
import factor.actions as actions

class facet_setup(operation_facet):

    def __init__(self, parset, direction):
        operation_facet(parset, direction)
        self.name = 'Facet setup'

    def run(self):

        mss = parset['mss']

        # add the central source to the visibilities
        actions = [ action.add_cal(ms, direction) for ms in mss ]
        self.s.run_action_parallel( actions )
        result = self.s.get_result()

        # phase shift in the DD-calibrator direction
        actions = [ action.phase_shift(ms, direction) for ms in mss ]
        self.s.run_action_parallel( actions ) 
        result = self.s.get_result()
        
        # averaging to...
        actions = [ action.avg(ms, time='1', freq='1') for ms in mss ]
        self.s.run_action_parallel( actions ) 
        result = self.s.get_result()
        
        # apply direction independent?
        actions = [ action.apply_die(ms, die_parmdb='') for ms in mss ]
        self.s.run_action_parallel( actions ) 
        result = self.s.get_result()

