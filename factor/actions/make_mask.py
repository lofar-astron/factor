"""
Action: make_mask
Make a mask using PyBDSM from a given image
Used parameters:
* threshpix
* threshisl
Return:
* mask name
"""

from factor.lib.action import action
from factor.lib.action_lib import makeimagebasename

class make_mask(action):
    """
    Implment the make_mask action
    """

    def __init__(self, op_name, p, prefix = None, direction = None, clean=True):
        super(make_mask, self).__init__(op_name, name = 'make_mask')
        self.imagebasename = makeimagebasename(self.ms, prefix, direction)

    def run(self):
        # run pybdsm....


    def get_results(self):
        """
        Return mask name
        """
        return self.mask
