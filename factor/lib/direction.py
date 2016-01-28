"""
Definition of the direction class
"""
import os
import logging
from astropy.coordinates import Angle
import numpy as np
from lsmtool.operations_lib import radec2xy
import matplotlib.path as mplPath
from scipy.ndimage import gaussian_filter


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
        RA in degrees of calibrator center
    dec : float
        Dec in degrees of calibrator center
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
        self.log = logging.getLogger('factor:{0}'.format(self.name))

        if type(ra) is str:
            ra = Angle(ra).to('deg').value
        if type(dec) is str:
            dec = Angle(dec).to('deg').value
        self.ra = ra
        self.dec = dec
        self.atrous_do = atrous_do
        self.mscale_field_do = mscale_field_do
        self.cal_imsize = cal_imsize
        self.solint_time_p = solint_p
        self.solint_time_a = solint_a
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
        self.is_outlier = outlier_do
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
        self.subtracted_data_colname = 'SUBTRACTED_DATA_ALL' # name of empty data column
        self.pre_average = False # whether to use baseline averaging
        self.blavg_weight_column = 'WEIGHT_SPECTRUM' # name of weights column
        self.started_operations = []
        self.completed_operations = []
        self.cleanup_mapfiles = []
        self.do_reset = False # whether to reset this direction
        self.is_patch = False # whether direction is just a patch (not full facet)
        self.nchunks = 1
        self.num_selfcal_groups = 1

        # Set the size of the calibrator (used to filter source lists)
        if cal_size_deg is None:
            # Try to get from cal_imsize assuming 50% padding
            if self.cal_imsize == 0:
                self.log.error('The cal_imsize must be specified in the directions '
                    'file')
                sys.exit(1)
            else:
                self.cal_size_deg = self.cal_imsize * self.cellsize_selfcal_deg / 1.5
        else:
            self.cal_size_deg = cal_size_deg
            if self.cal_imsize == 0:
                self.cal_imsize = max(512, self.get_optimum_size(self.cal_size_deg
                    / self.cellsize_selfcal_deg * 1.2)) # cal size has 20% padding

        self.cal_radius_deg = self.cal_size_deg / 2.0
        self.cal_rms_box = self.cal_size_deg / self.cellsize_selfcal_deg

        # Define some directories and files
        self.working_dir = factor_working_dir
        self.save_file = os.path.join(self.working_dir, 'state',
            self.name+'_save.pkl')
        self.vertices_file = self.save_file


    def set_imaging_parameters(self, nbands, nbands_per_channel,
        initial_skymodel=None, test_run=False):
        """
        Sets various parameters for images in facetselfcal and facetimage pipelines

        Parameters
        ----------
        nbands : int
            Number of bands
        nbands_per_channel : int
            Number of bands per output channel (WSClean only)
        initial_skymodel : LSMTool SkyModel object, optional
            Sky model used to check source sizes
        test_run : bool, optional
            If True, use test sizes

        """
        if not test_run:
            # Set facet image size
            if self.facet_imsize == 0:
                if hasattr(self, 'width'):
                    self.facet_imsize = max(512, self.get_optimum_size(self.width
                        / self.cellsize_selfcal_deg * 1.3)) # full facet has 30% padding
                else:
                    self.facet_imsize = None
        else:
            self.facet_imsize = self.get_optimum_size(128)
            self.cal_imsize = self.get_optimum_size(128)

        self.cal_wplanes = self.set_wplanes(self.cal_imsize)
        self.facet_wplanes = self.set_wplanes(self.facet_imsize)

        # Determine whether the total bandwidth is large enough that wide-band
        # imaging is needed
        if nbands > 5:
            self.use_wideband = True
        else:
            self.use_wideband = False

        # Set number of channels for wide-band imaging with WSClean and nterms
        # for the CASA imager. Also define the image suffixes (which depend on
        # whether or not wide-band clean is done)
        if self.use_wideband:
            self.nchannels = int(round(float(nbands)/
                float(nbands_per_channel)))
            self.nterms = 2
            self.casa_suffix = '.tt0'
            self.wsclean_suffix = '-MFS-image.fits'
        else:
            self.nchannels = 1
            self.nterms = 1
            self.casa_suffix = None
            self.wsclean_suffix = '-image.fits'

        # Set number of iterations and threshold for full facet image, scaled to
        # the number of bands
        scaling_factor = np.sqrt(np.float(nbands))
        self.wsclean_full_image_niter = int(5000 * scaling_factor)
        self.wsclean_full_image_threshold_jy =  1.5e-3 * 0.7 / scaling_factor
        self.casa_full_image_niter = int(2000 * scaling_factor)
        self.casa_full_image_threshold_mjy = "{}mJy".format(1.5 * 0.7 / scaling_factor)

        # Set multiscale imaging mode: Get source sizes and check for large
        # sources (anything above 2 arcmin -- the CC sky model was convolved
        # with a Gaussian of 1 arcmin, so unresolved sources have sizes of ~
        # 1 arcmin)
        if initial_skymodel is not None:
            sizes = self.get_source_sizes(initial_skymodel)
            large_size_arcmin = 2.0
            if any([s > large_size_arcmin for s in sizes]):
                self.mscale_field_do = True
            else:
                self.mscale_field_do = False
            if self.mscale_field_do:
                self.casa_multiscale = '[0, 3, 7, 25, 60, 150]'
                self.wsclean_multiscale = '-multiscale,'
            else:
                self.casa_multiscale = '[0]'
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


    def get_source_sizes(self, skymodel):
        """
        Returns list of source sizes in arcmin

        Parameters
        ----------
        skymodel : LSMTool SkyModel object
            CC sky model used to determine source sizes. The sky model is
            filtered to include only those sources within the direction
            facet

        """
        x, y, midRA, midDec = skymodel._getXY()
        xv, yv = radec2xy(self.vertices[0], self.vertices[1], midRA, midDec)
        xyvertices = np.array([[xp, yp] for xp, yp in zip(xv, yv)])
        bbPath = mplPath.Path(xyvertices)
        inside = np.zeros(len(skymodel), dtype=bool)
        for i in range(len(skymodel)):
            inside[i] = bbPath.contains_point((x[i], y[i]))
        skymodel.select(inside, force=True)
        sizes = skymodel.getPatchSizes(units='arcmin', weight=True)

        return sizes


    def get_cal_fluxes(self, skymodel, fwhmArcsec=25.0, threshold=0.1):
        """
        Returns total flux density in Jy and max peak flux density in
        Jy per beam for calibrator

        Parameters
        ----------
        skymodel : LSMTool SkyModel object
            CC sky model used to determine source fluxes. The sky model is
            filtered to include only those sources within the calibrator region

        Returns
        -------
        tot_flux_jy, peak_flux_jy_bm : float, float
            Total flux density in Jy and max peak flux density in
            Jy per beam for calibrator

        """
        dist = skymodel.getDistance(self.ra, self.dec)
        skymodel.select(dist < self.cal_radius_deg)

        # Generate image grid with 1 pix = FWHM / 4
        x, y, midRA, midDec  = skymodel._getXY(crdelt=fwhmArcsec/4.0/3600.0)
        fluxes_jy = skymodel.getColValues('I', units='Jy')
        sizeX = int(np.ceil(1.2 * (max(x) - min(x)))) + 1
        sizeY = int(np.ceil(1.2 * (max(y) - min(y)))) + 1
        image = np.zeros((sizeX, sizeY))
        xint = np.array(x, dtype=int)
        xint += -1 * min(xint)
        yint = np.array(y, dtype=int)
        yint += -1 * min(yint)
        for xi, yi, f in zip(xint, yint, fluxes_jy):
            image[xi, yi] = f

        # Convolve with Gaussian of FWHM = 4 pixels
        image_blur = gaussian_filter(image, [4.0/2.35482, 4.0/2.35482])
        beam_area_pix = 1.1331*(4.0)**2

        return np.sum(fluxes_jy), np.max(image_blur)*beam_area_pix


    def set_averaging_steps_and_solution_intervals(self, chan_width_hz, nchan,
        timestep_sec, ntimes, nbands, initial_skymodel=None,
        preaverage_flux_jy=0.0):
        """
        Sets the averaging step sizes and solution intervals

        The solution-interval scaling is done so that sources with total flux
        densities below 1.4 Jy at the highest frequency have a fast interval of
        8 time slots and a slow interval of 240 time slots for a bandwidth of 4
        bands. The fast intervals are scaled with the bandwidth and flux as
        nbands^-0.5 and flux^2. The slow intervals are scaled as flux^2.

        When multiple sources are combined into a single calibrator, the flux
        density of each source is obviously lower than the total, and hence
        the model will be less well constrained. To compensate for this effect,
        we scale the solution intervals by the number of sources in the
        calibrator.

        Note: the frequency step for averaging must be an even divisor of the
        number of channels

        Parameters
        ----------
        chan_width_hz : float
            Channel width in Hz
        nchan : int
            Number of channels per band
        timestep_sec : float
            Time step
        ntimes : int
            Number of timeslots per band
        nbands : int
            Number of bands
        initial_skymodel : LSMTool SkyModel object, optional
            Sky model used to check source sizes
        preaverage_flux_jy : bool, optional
            Use baseline-dependent averaging and solint_time_p = 1 for phase-only
            calibration for sources below this flux value

        """
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
        #
        # Calculate the effective flux density. This is the one used to set the
        # intervals. It is the peak flux density adjusted to account for cases
        # in which the total flux density is larger than the peak flux density
        # would indicate (either due to source being extended or to multiple
        # calibrator sources). In these cases, we can use a higher effective
        # flux density to set the intervals. A scaling with a power of 1/1.5
        # seems to work well
        if initial_skymodel is not None:
            # The initial skymodel is not used for field directions, so steps
            # below are skipped
            total_flux_jy, peak_flux_jy_bm = self.get_cal_fluxes(initial_skymodel)
            effective_flux_jy = peak_flux_jy_bm * (total_flux_jy / peak_flux_jy_bm)**0.667
            ref_flux_jy = 1.4 * (4.0 / nbands)**0.5
            self.log.debug('Total flux density of calibrator: {} Jy'.format(total_flux_jy))
            self.log.debug('Peak flux density of calibrator: {} Jy/beam'.format(peak_flux_jy_bm))
            self.log.debug('Effective flux density of calibrator: {} Jy'.format(effective_flux_jy))

            # Set baseline-dependent pre-averaging flag
            if effective_flux_jy < preaverage_flux_jy:
                self.pre_average = True
            else:
                self.pre_average = False

            # Set fast (phase-only) solution interval
            if self.solint_time_p == 0:
                if self.pre_average:
                    # Set solution interval to 1 timeslot and vary the target rms per
                    # solution interval instead (which affects the width of the
                    # preaveraging Gaussian)
                    self.solint_time_p = 1
                    self.target_rms_rad = int(round(0.5 * (ref_flux_jy / effective_flux_jy)**2))
                    if self.target_rms_rad < 0.2:
                        self.target_rms_rad = 0.2
                    if self.target_rms_rad > 0.5:
                        self.target_rms_rad = 0.5
                else:
                    self.solint_time_p = int(round(8 * (ref_flux_jy / effective_flux_jy)**2))
                    if self.solint_time_p < 1:
                        self.solint_time_p = 1
                    if self.solint_time_p > 2:
                        self.solint_time_p = 2

            # Set slow (gain) solution interval
            if self.solint_time_a == 0:
                # Amplitude solve is per band, so don't scale with number of bands
                ref_flux = 1400.0
                self.solint_time_a = int(round(240 * (ref_flux_jy / effective_flux_jy)**2))
                if self.solint_time_a < 30:
                    self.solint_time_a = 30
                if self.solint_time_a > 120:
                    self.solint_time_a = 120

            self.log.debug('Using solution intervals of {0} (fast) and {1} '
                '(slow) time slots'.format(self.solint_time_p, self.solint_time_a))

            # Set chunk width for time chunking to the amplitude solution time
            # interval (minus one time slot to ensure that we don't get a very short
            # solution interval at the end of the chunk) so that it's close to
            # ~ 200 time slots (to avoid memory/performance issues)
            self.chunk_width = self.solint_time_a - 1
            while self.chunk_width < 200:
                self.chunk_width += self.solint_time_a - 1
            self.nchunks = int(np.ceil((np.float(ntimes) / np.float(self.chunk_width))))

            # Set frequency interval for selfcal solve steps. The interval for
            # slow (amp) selfcal should be the number of channels in a band after
            # averaging. The interval for fast (phase) selfcal should be the
            # number of channels in 20 MHz or less
            num_chan_per_band_after_avg = nchan / self.facetselfcal_freqstep
            self.solint_freq_a = num_chan_per_band_after_avg
            num_cal_blocks = np.ceil(nchan * nbands * chan_width_hz/1e6
                / 20.0)
            self.solint_freq_p = int(np.ceil(num_chan_per_band_after_avg * nbands
                / num_cal_blocks))

        # Set name of column to use for data and averaged weights
        if self.pre_average:
            self.data_column = 'BLAVG_DATA'
            self.blavg_weight_column = 'BLAVG_WEIGHT_SPECTRUM'
        else:
            self.data_column = 'DATA'
            self.blavg_weight_column = 'WEIGHT_SPECTRUM'


    def save_state(self):
        """
        Saves the direction state to a file

        """
        import pickle

        with open(self.save_file, 'wb') as f:
            save_dict = self.__dict__.copy()
            save_dict.pop('log')
            pickle.dump(save_dict, f)


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
                if ('diff_models_field_mapfile' in d and
                    'input_bands_mapfile' in d and
                    'subtracted_data_colname' in d):
                    self.diff_models_field_mapfile = d['diff_models_field_mapfile']
                    self.input_bands_mapfile = d['input_bands_mapfile']
                    self.subtracted_data_colname = d['subtracted_data_colname']
            return True
        except:
            return False


    def reset_state(self, op_names=None):
        """
        Resets the direction to allow reprocessing

        Currently, this means just deleting the results directories,
        but it could be changed to delete only a subset of selfcal steps (by
        modifying the selfcal pipeline statefile).

        Parameters
        ----------
        op_names : list of str, optional
            Name of operation to reset. If None, all started and completed
            operations are reset

        """
        if op_names is None:
            op_names = self.completed_operations[:] + self.started_operations[:]
        elif type(op_names) is str:
            op_names = [op_names]
        self.log.info('Resetting state for operation(s): {}'.format(', '.join(op_names)))

        # Reset selfcal flag
        if 'facetselfcal' in op_names:
            self.selfcal_ok = False

        # Remove operation name from lists of started and completed operations
        # and delete the results directories
        for op_name in op_names:
            while op_name in self.completed_operations:
                self.completed_operations.remove(op_name)
            while op_name in self.started_operations:
                self.started_operations.remove(op_name)

            # Delete results directory for this operation
            op_dir = os.path.join(self.working_dir, 'results', op_name, self.name)
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
