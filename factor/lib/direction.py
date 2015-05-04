"""
Definition of the direction class
"""
import os
import logging

class Direction(object):
    """
    Generic direction class
    """
    def __init__(self, name, ra, dec, atrous_do, mscale_field_do, cal_imsize,
        solint_p, solint_a, field_imsize, dynamic_range, region_selfcal,
        region_field, peel_skymodel, outlier_do, factor_working_dir,
        make_final_image=False, cal_radius_deg=None, cal_flux_jy=None):
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
        atrous_do : bool
            Fit to wavelet images in PyBDSM?
        mscale_field_do : bool
            Use multiscale clean for facet field?
        cal_imsize : int
            Size of calibrator image in 1.5 arcsec pixels
        solint_p : int
            Solution interval for phase calibration (# of time slots)
        solint_a : int
            Solution interval for amplitude calibration (# of time slots)
        factor_working_dir : str
            Full path of working directory
        make_final_image : bool
            Make final image of this direction, after all directions have been
            selfcaled?
        cal_radius_deg : float
            Radius in degrees of calibrator source
        cal_flux_jy : float
            Apparent flux in Jy of calibrator source
        """
        self.name = name
        self.ra = ra
        self.dec = dec
        self.atrous_do = atrous_do
        self.mscale_field_do = mscale_field_do
        self.cal_imsize = cal_imsize
        self.cal_radius_deg = cal_imsize * 1.5 / 3660.0
        self.solint_p = solint_p
        self.solint_a = solint_a
        self.field_imsize = field_imsize
        self.region_selfcal = region_selfcal
        if self.region_selfcal.lower() == 'empty':
            self.region_selfcal = None
        self.region_field = region_field
        if self.region_field.lower() == 'empty':
            self.region_field = None
        self.peel_skymodel = peel_skymodel
        if self.peel_skymodel.lower() == 'empty':
            self.peel_skymodel = None
        self.make_final_image = make_final_image
        self.dynamic_range = dynamic_range
        self.apparent_flux_mjy = cal_flux_jy * 1000.0
        self.nchannels = 1

        self.completed_operations = []
        self.save_file = os.path.join(factor_working_dir, 'state',
            self.name+'_save.pkl')


    def save_state(self):
        """
        Saves the direction state to a file
        """
        import pickle

        with open(self.save_file, 'wb') as f:
            pickle.dump(self.__dict__, f)


    def load_state(self):
        """
        Loads the direction state from a file

        Returns
        -------
        success : bool
            True if state was successfully loaded, False if not
        """
        import pickle

        try:
            with open(self.save_file, 'r') as f:
                self.__dict__ = pickle.load(f)
            return True
        except:
            return False
