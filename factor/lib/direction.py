"""
Definition of the direction class
"""
import logging

class Direction(object):
    """
    Generic direction class
    """
    def __init__(self, name, ra, dec, clean_reg, multiscale, solint_a, solint_p,
        make_final_image, cal_radius, apparent_flux):
        """
        Create Direction object

        Parameters
        ----------
        name : str
            Name of direction
        ra : float
            RA in degrees of direction center
        dec : float
            Dec in degrees of direction center
        clean_reg : str
            CASA region file for clean
        multiscale : bool
            Use multiscale clean?
        solint_a : int
            Solution interval for amplitude calibration (# of time slots)
        solint_p : int
            Solution interval for phase calibration (# of time slots)
        make_final_image : bool
            Make final image of this direction, after all directions have been
            selfcaled?
        cal_radius : float
            Radius in degrees of calibrator source
        apparent_flux : float
            Apparent flux in mJy of calibrator source
        """
        self.name = name
        self.ra = ra
        self.dec = dec
        self.clean_reg = clean_reg
        self.multiscale = multiscale
        self.solint_a = solint_a
        self.solint_p = solint_p
        self.make_final_image = make_final_image
        self.cal_radius_deg = cal_radius / 60.0
        self.apparent_flux_mjy = apparent_flux
        self.nchannels = 1
