"""
Definition of the direction class
"""
import os
import logging
from astropy.coordinates import Angle
import numpy as np
from lsmtool.operations_lib import radec2xy
import matplotlib.path as mplPath


class Direction(object):
    """
    Generic direction class

    A direction object holds all the parameters needed for an operation in a
    given direction (facet).

    Note:
    All attributes needed by the pipeline templates should be set on the class
    instance so that they can be passed with self.__dict__

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
    def __init__(self, name, ra, dec, atrous_do=False, mscale_field_do=False, cal_imsize=512,
        solint_p=1, solint_a=30, field_imsize=2048, dynamic_range='LD', region_selfcal='',
        region_field='', peel_skymodel='', outlier_do=False, factor_working_dir='',
        make_final_image=False, cal_size_deg=None, cal_flux_jy=None):
        # Handle input args
        self.name = name
        if type(ra) is str:
            ra = Angle(ra).to('deg').value
        if type(dec) is str:
            dec = Angle(dec).to('deg').value
        self.ra = ra
        self.dec = dec
        self.atrous_do = atrous_do
        self.mscale_field_do = mscale_field_do
        self.cal_imsize = cal_imsize
        self.solint_p = solint_p
        self.solint_a = solint_a
        self.facet_imsize = field_imsize * 1.15
        self.dynamic_range = dynamic_range
        if region_selfcal.lower() == 'empty':
            # Set to empty list (casapy format)
            self.region_selfcal = '[]'
        else:
            self.region_selfcal = '["{0}"]'.format(region_selfcal)
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

        # Initialize some parameters to default/initial values
        self.loop_amp_selfcal = False
        self.selfcal_ok = False # whether selfcal succeeded
        self.skip_add_subtract = None # whether to skip add/subtract in facetsub op
        self.max_residual_val = 0.5 # maximum residual in Jy for facet subtract test
        self.nchannels = 1 # set number of wide-band channels
        self.use_new_sub_data = False # set flag that tells which subtracted-data column to use
        self.cellsize_selfcal_deg = 0.000417 # selfcal cell size
        self.cellsize_verify_deg = 0.00833 # verify subtract cell size
        self.target_rms_rad = 0.2 # preaverage target rms
        self.subtracted_data_colname = 'SUBTRACTED_DATA_ALL'
        self.pre_average = False
        self.blavg_weight_column = 'WEIGHT_SPECTRUM'
        self.started_operations = []
        self.completed_operations = []
        self.cleanup_mapfiles = []

        # Set the size of the calibrator (used to filter source lists)
        if cal_size_deg is None:
            # Get from cal imsize assuming 50% padding
            self.cal_size_deg = cal_imsize * self.cellsize_selfcal_deg / 1.5
        else:
            self.cal_size_deg = cal_size_deg
        self.cal_radius_deg = self.cal_size_deg / 2.0
        self.cal_rms_box = self.cal_size_deg / self.cellsize_selfcal_deg

        # Define some directories and files
        self.working_dir = factor_working_dir
        self.save_file = os.path.join(self.working_dir, 'state',
            self.name+'_save.pkl')
        self.vertices_file = self.save_file


    def set_imaging_parameters(self, nbands, nbands_per_channel, initial_skymodel,
        test_run=False):
        """
        Sets various parameters for images in facetselfcal and facetimage pipelines

        Parameters
        ----------
        nbands : int
            Number of bands
        nbands_per_channel : int
            Number of bands per output channel (WSClean only)
        initial_skymodel : LSMTool SkyModel object
            Sky model used to check source sizes
        test_run : bool, optional
            If True, use test sizes

        """
        if not test_run:
            # Set selfcal and facet image sizes
            if hasattr(self, 'width'):
                self.facet_imsize = max(512, self.get_optimum_size(self.width
                    / self.cellsize_selfcal_deg * 1.3)) # full facet has 30% padding
            else:
                self.facet_imsize = None
            self.cal_imsize = max(512, self.get_optimum_size(self.cal_size_deg
                / self.cellsize_selfcal_deg * 1.2)) # cal size has 20% padding
        else:
            self.facet_imsize = self.get_optimum_size(128)
            self.cal_imsize = self.get_optimum_size(128)

        self.cal_wplanes = self.set_wplanes(self.cal_imsize)
        self.facet_wplanes = self.set_wplanes(self.facet_imsize)

        # Set number of channels for wide-band imaging
        if nbands > 5:
            self.nchannels = int(round(float(nbands)/
                float(nbands_per_channel)))
        else:
            self.nchannels = 1

        if self.nchannels > 1:
            self.nterms = 2
            self.casa_suffix = '.tt0'
            self.wsclean_suffix = '-MFS-image.fits'
        else:
            self.nterms = 1
            self.casa_suffix = None
            self.wsclean_suffix = '-image.fits'

        # Set number of iterations for full facet image, scaled to the number
        # of bands
        self.wsclean_full_image_niter = int(5000 * (np.sqrt(np.float(nbands))))
        self.wsclean_full_image_threshold_jy =  1.5e-3 * 0.7 / (np.sqrt(np.float(nbands)))
        self.casa_full_image_niter = int(2000 * (np.sqrt(np.float(nbands))))
        self.casa_full_image_threshold_mjy = "{}mJy".format(1.5 * 0.7 / (np.sqrt(np.float(nbands))))

        # Set multiscale imaging mode
        x, y, midRA, midDec = initial_skymodel._getXY()
        xv, yv = radec2xy(self.vertices[0], self.vertices[1], midRA, midDec)
        xyvertices = np.array([[xp, yp] for xp, yp in zip(xv, yv)])
        bbPath = mplPath.Path(xyvertices)
        inside = np.zeros(len(initial_skymodel), dtype=bool)
        for i in range(len(initial_skymodel)):
            inside[i] = bbPath.contains_point((x[i], y[i]))
        initial_skymodel.select(inside, force=True)

        # Get source sizes and check for large sources (anything above 1 arcmin)
        large_size_arcmin = 1.0
        sizes = s.getPatchSizes(units='arcmin')
        if any([s > large_size_arcmin for s in sizes]):
            self.mscale_field_do = True

        if self.mscale_field_do:
            self.casa_multiscale = '[0, 3, 7, 25, 60, 150]'
            self.wsclean_multiscale = '-multiscale,'
        else:
            self.casa_multiscale = '[]'
            self.wsclean_multiscale = ''


    def set_wplanes(self, imsize):
        """
        Sets number of wplanes for casa clean

        Parameters
        ----------
        imsize : int
            Image size in pixels

        """
        wplanes = 1
        if imsize > 512:
            wplanes = 64
        if imsize > 799:
            wplanes = 96
        if imsize > 1023:
            wplanes = 128
        if imsize > 1599:
            wplanes = 256
        if imsize > 2047:
            wplanes = 384
        if imsize > 3000:
            wplanes = 448
        if imsize > 4095:
            wplanes = 512

        return wplanes


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


    def set_averaging_steps_and_solution_intervals(self, chan_width_hz, nchan,
        timestep_sec, ntimes, nbands, pre_average=False):
        """
        Sets the averaging step sizes and solution intervals for selfcal

        The solution-interval scaling is done so that sources with total flux
        densities below 1.4 Jy at the highest frequency have a fast interval of
        8 time slots and a slow interval of 240 time slots for a bandwidth of 4
        bands. The fast intervals are scaled with the bandwidth and flux as
        nbands^-0.5 and flux^2. The slow intervals are scaled as flux^2.

        Note: the frequency step for averaging must be an even divisor of the
        number of channels

        Parameters
        ----------
        chan_width_hz : float
            Channel width
        nchan : int
            Number of channels per band
        timestep_sec : float
            Time step
        ntimes : int
            Number of timeslots per band
        nbands : int
            Number of bands
        pre_average : bool, optional
            If True, use baseline-dependent averaging and solint_p = 1 for
            phase-only calibration

        """
        # Set name of column to use for data and averaged weights
        self.pre_average = pre_average
        if self.pre_average:
            self.data_column = 'BLAVG_DATA'
            self.blavg_weight_column = 'BLAVG_WEIGHT_SPECTRUM'
        else:
            self.data_column = 'DATA'
            self.blavg_weight_column = 'WEIGHT_SPECTRUM'

        # For initsubtract, average to 0.5 MHz per channel and 20 sec per time
        # slot. Since each band is imaged separately and the smearing and image
        # sizes both scale linearly with frequency, a single frequency and time
        # step is valid for all bands
        self.initsubtract_freqstep = max(1, min(int(round(0.5 * 1e6 / chan_width_hz)), nchan))
        while nchan % self.initsubtract_freqstep:
            self.initsubtract_freqstep += 1
        self.initsubtract_timestep = max(1, int(round(20.0 / timestep_sec)))

        # For selfcal, average to 2 MHz per channel and 120 s per time slot for
        # an image of 512 pixels
        target_bandwidth_mhz = 2.0 * 512.0 / self.cal_imsize
        target_timewidth_s = 120 * 512.0 / self.cal_imsize # used for imaging only
        self.facetselfcal_freqstep = max(1, min(int(round(target_bandwidth_mhz * 1e6 / chan_width_hz)), nchan))
        while nchan % self.facetselfcal_freqstep:
            self.facetselfcal_freqstep += 1
        self.facetselfcal_timestep = max(1, int(round(target_timewidth_s / timestep_sec)))

        # For facet imaging, average to 0.5 MHz per channel and 30 sec per time
        # slot for an image of 2048 pixels
        target_bandwidth_mhz = 0.5 * 2048.0 / self.facet_imsize
        target_timewidth_s = 30 * 2048.0 / self.facet_imsize
        self.facetimage_freqstep = max(1, min(int(round(target_bandwidth_mhz * 1e6 / chan_width_hz)), nchan))
        while nchan % self.facetimage_freqstep:
            self.facetimage_freqstep += 1
        self.facetimage_timestep = max(1, int(round(target_timewidth_s / timestep_sec)))

        # For selfcal verify, average to 2 MHz per channel and 60 sec per time
        # slot
        self.verify_freqstep = max(1, min(int(round(2.0 * 1e6 / chan_width_hz)), nchan))
        while nchan % self.verify_freqstep:
            self.verify_freqstep += 1
        self.verify_timestep = max(1, int(round(60.0 / timestep_sec)))

        # Set time intervals for selfcal solve steps
        if self.apparent_flux_mjy is not None:
            ref_flux = 1400.0 * (4.0 / nbands)**0.5
            if self.pre_average:
                # Set solution interval to 1 timeslot and vary the target rms per
                # solution interval instead (which affects the width of the
                # preaveraging Gaussian)
                self.solint_p = 1
                self.target_rms_rad = int(round(0.5 * (ref_flux / self.apparent_flux_mjy)**2))
                if self.target_rms_rad < 0.2:
                    self.target_rms_rad = 0.2
                if self.target_rms_rad > 0.5:
                    self.target_rms_rad = 0.5
            else:
                self.solint_p = int(round(8 * (ref_flux / self.apparent_flux_mjy)**2))
                if self.solint_p < 1:
                    self.solint_p = 1
                if self.solint_p > 8:
                    self.solint_p = 8

            # Amplitude solve is per band, so don't scale with number of bands
            ref_flux = 1400.0
            self.solint_a = int(round(240 * (ref_flux / self.apparent_flux_mjy)**2))
            if self.solint_a < 30:
                self.solint_a = 30
            if self.solint_a > 240:
                self.solint_a = 240

        # Set chunk width for time chunking to the amplitude solution time
        # interval (minus one time slot to ensure that we don't get a very short
        # solution interval at the end of the chunk) so that it's close to
        # ~ 200 time slots (to avoid memory/performance issues)
        self.chunk_width = self.solint_a - 1
        while self.chunk_width < 200:
            self.chunk_width += self.solint_a - 1

        # Set frequency interval for selfcal solve steps to the number of
        # channels in a band after averaging
        self.solint_freq = nchan / self.facetselfcal_freqstep


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

        Note: only state attributes are loaded to avoid overwritting
        non-state attributes

        Returns
        -------
        success : bool
            True if state was successfully loaded, False if not
        """
        import pickle

        try:
            with open(self.save_file, 'r') as f:
                d = pickle.load(f)

                # Load list of started operations
                if 'started_operations' in d:
                    self.started_operations = d['started_operations']

                # Load list of completed operations
                if 'completed_operations' in d:
                    self.completed_operations = d['completed_operations']

                # Load mapfiles needed for facetsubreset
                if ('diff_models_field_datamap' in d and
                    'input_bands_datamap' in d and
                    'subtracted_data_colname' in d):
                    self.diff_models_field_datamap = d['diff_models_field_datamap']
                    self.input_bands_datamap = d['input_bands_datamap']
                    self.subtracted_data_colname = d['subtracted_data_colname']
            return True
        except:
            return False


    def reset_state(self):
        """
        Resets the direction to allow reprocessing

        Currently, this means just deleting the results directories,
        but it could be changed to delete only a subset of selfcal steps (by
        modifying the selfcal pipeline statefile).
        """
        # Reset selfcal flag
        self.selfcal_ok = False

        # Remove operation name from lists of started and completed operations
        # and delete the results directories
        for op in set(self.completed_operations[:] + self.started_operations[:]):
            if op in self.completed_operations:
                self.completed_operations.remove(op)
            if op in self.started_operations:
                self.started_operations.remove(op)

            # Delete results directory for this operation
            op_dir = os.path.join(self.working_dir, 'results', op, self.name)
            if os.path.exists(op_dir):
                os.system('rm -rf {0}'.format(op_dir))

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
                    # Handle case in which item.file is a Python list
                    if item.file[0] == '[' and item.file[-1] == ']':
                        files = item.file.strip('[]').split(',')
                    else:
                        files = [item.file]
                    for f in files:
                        if os.path.exists(f):
                            os.system('rm -rf {0}'.format(f))
            except IOError:
                pass
