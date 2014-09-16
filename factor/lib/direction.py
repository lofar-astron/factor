"""
Definition of direction-related classes
"""
import logging

class direction( object ):
    """
    general direciton class
    """
    def __init__(self, name, ra, dec, reg, multiscale, solint_a, solint_p, hdr):
        logging.debug("Setting up direction %s." % name)
        self.name = name
        self.ra = ra
        self.dec = dec
        self.reg = reg
        self.multiscale = multiscale
        self.solint_a = solint_a
        self.solint_p = solint_p
        self.hdr = hdr

