"""
Definition of the direction class
"""
import os
import logging

class Direction(object):
    """
    Generic direction class
    """
    def __init__(self, name, ra, dec, atrous_do=False, mscale_field_do=False, cal_imsize=0,
        solint_p=0, solint_a=0, field_imsize=0, dynamic_range='LD', region_selfcal='',
        region_field='', peel_skymodel='', outlier_do=False, factor_working_dir='',
        make_final_image=False, cal_size_deg=None, cal_flux_jy=None):
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
        field_imsize : int
            Size of facet image in 1.5 arcsec pixels
        dynamic_range : str
            LD (low dynamic range) or HD (high dynamic range)
        region_selfcal : str
            Region for clean mask for calibrator selfcal
        region_field : str
            Region for clean mask for facet image
        peel_skymodel : str
            Sky model for peeling
        outlier_do : bool
            If True, peel source without selfcal
        factor_working_dir : str
            Full path of working directory
        make_final_image : bool, optional
            Make final image of this direction, after all directions have been
            selfcaled?
        cal_size_deg : float, optional
            Size in degrees of calibrator source(s)
        cal_flux_jy : float, optional
            Apparent flux in Jy of calibrator source
        """
        # Handle input args
        self.name = name
        self.ra = ra
        self.dec = dec
        self.atrous_do = atrous_do
        self.mscale_field_do = mscale_field_do
        self.cal_imsize = cal_imsize
        self.solint_p = solint_p
        self.solint_a = solint_a
        self.facet_imsize = field_imsize * 1.15
        self.dynamic_range = dynamic_range
        self.region_selfcal = region_selfcal
        if self.region_selfcal.lower() == 'empty':
            self.region_selfcal = None
        self.region_field = region_field
        if self.region_field.lower() == 'empty':
            self.region_field = None
        self.peel_skymodel = peel_skymodel
        if self.peel_skymodel.lower() == 'empty':
            self.peel_skymodel = None
        self.outlier_do = outlier_do
        self.make_final_image = make_final_image
        if cal_flux_jy is not None:
            self.apparent_flux_mjy = cal_flux_jy * 1000.0
        else:
            self.apparent_flux_mjy = None

        # Initialize some parameters
        self.loop_amp_selfcal = False
        self.selfcal_ok = True # whether selfcal is still improving after first amp cal
        self.max_residual_val = 0.5 # maximum residual in Jy for facet subtract test
        self.nchannels = 1 # set number of wide-band channels
        self.use_new_sub_data = False # set flag that tells which subtracted-data column to use
        self.cellsize_selfcal_deg = 0.000417 # selfcal cell size
        self.cellsize_verify_deg = 0.00833 # verify subtract cell size

        # Set the size of the calibrator (used to filter source lists)
        if cal_size_deg is None:
            # Get from cal imsize assuming 50% padding
            self.cal_size_deg = cal_imsize * self.cellsize_selfcal_deg / 1.5
        else:
            self.cal_size_deg = cal_size_deg
        self.cal_rms_box = self.cal_size_deg / self.cellsize_selfcal_deg

        # Set number of wplanes for casapy imaging
        self.wplanes = 1
        if self.cal_imsize > 512:
            self.wplanes = 64
        if self.cal_imsize > 799:
            self.wplanes = 96
        if self.cal_imsize > 1023:
            self.wplanes = 128
        if self.cal_imsize > 1599:
            self.wplanes = 256
        if self.cal_imsize > 2047:
            self.wplanes = 384
        if self.cal_imsize > 3000:
            self.wplanes = 448
        if self.cal_imsize > 4095:
            self.wplanes = 512

        # Scale solution intervals by apparent flux. The scaling is done so that
        # sources with flux densities of 100 mJy have a fast interval of 4 time
        # slots and a slow interval of 240 time slots. The scaling is currently
        # linear with flux (and thus we accept lower-SNR solutions for the
        # fainter sources). Ideally, these value should also scale with the
        # bandwidth
        if self.apparent_flux_mjy is not None:
            ref_flux = 100.0
            self.solint_p = max(1, int(round(4 * ref_flux / self.apparent_flux_mjy)))
            self.solint_a = max(30, int(round(240 * ref_flux / self.apparent_flux_mjy)))

        # Define some directories, etc.
        self.working_dir = factor_working_dir
        self.completed_operations = []
        self.cleanup_mapfiles = []
        self.save_file = os.path.join(self.working_dir, 'state',
            self.name+'_save.pkl')
        self.pipeline_dir = os.path.join(self.working_dir, 'pipeline')


    def set_image_sizes(self, test_run=False):
        """
        Sets sizes for various images

        Parameters
        ----------
        test_run : bool, optional
            If True, use test sizes

        """
        if not test_run:
            # Set selfcal and facet image sizes
            if hasattr(self, 'width'):
                self.facet_imsize = max(512, self.get_optimum_size(self.width
                    / self.cellsize_selfcal_deg * 1.15)) # full facet has 15% padding
            else:
                self.facet_imsize = None
            self.cal_imsize = max(512, self.get_optimum_size(self.cal_size_deg
                / self.cellsize_selfcal_deg * 1.2)) # cal size has 20% padding
        else:
            self.facet_imsize = self.get_optimum_size(128)
            self.cal_imsize = self.get_optimum_size(128)


    def get_optimum_size(self, size):
        """
        Gets the nearest optimum image size

        Taken from the casa source code (cleanhelper.py)

        Parameters
        ----------
        size : int
            Target image size in pixels

        Returns
        -------
        optimum_size : int
            Optimum image size nearest to target size

        """
        import numpy

        def prime_factors(n, douniq=True):
            """ Return the prime factors of the given number. """
            factors = []
            lastresult = n
            sqlast=int(numpy.sqrt(n))+1
            if n == 1:
                return [1]
            c=2
            while 1:
                 if (lastresult == 1) or (c > sqlast):
                     break
                 sqlast=int(numpy.sqrt(lastresult))+1
                 while 1:
                     if(c > sqlast):
                         c=lastresult
                         break
                     if lastresult % c == 0:
                         break
                     c += 1

                 factors.append(c)
                 lastresult /= c

            if (factors==[]): factors=[n]
            return  numpy.unique(factors).tolist() if douniq else factors

        n = int(size)
        if (n%2 != 0):
            n+=1
        fac=prime_factors(n, False)
        for k in range(len(fac)):
            if (fac[k] > 7):
                val=fac[k]
                while (numpy.max(prime_factors(val)) > 7):
                    val +=1
                fac[k]=val
        newlarge=numpy.product(fac)
        for k in range(n, newlarge, 2):
            if ((numpy.max(prime_factors(k)) < 8)):
                return k
        return newlarge


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


    def reset_state(self):
        """
        Resets the direction to initial state to allow reprocessing

        Currently, this means just deleting the facetselfcal results directory,
        but it could be changed to delete only a subset of selfcal steps (by
        modifying the selfcal pipeline statefile).
        """
        import glob

        operations = ['facetselfcal']
        for op in operations:
            # Remove entry in completed_operations
            self.completed_operations.remove(op)

            # Delete pipeline state
            action_dirs = glob.glob(os.path.join(self.pipeline_dir, op, '*'))
            for action_dir in action_dirs:
                facet_dir = os.path.join(action_dir, self.name)
                if os.path.exists(facet_dir):
                    os.system('rm -rf {0}'.format(facet_dir))

        self.save_state()


    def cleanup(self):
        """
        Cleans up unneeded data
        """
        from lofarpipe.support.data_map import DataMap

        for mapfile in self.cleanup_mapfiles:
            try:
                datamap = DataMap.load(mapfile)
                for item in datamap:
                    if os.path.exists(item.file):
                        os.system('rm -rf {0}'.format(item.file))
            except IOError:
                pass
