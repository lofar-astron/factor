"""
Definition of the direction class
"""
import os
import logging
from astropy.coordinates import Angle
import numpy as np
import lsmtool
from lsmtool.operations_lib import radec2xy
import matplotlib.path as mplPath
from scipy.ndimage import gaussian_filter
from scipy.special import erf
import sys
import glob


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
    cal_size_deg : float, optional
        Size in degrees of calibrator source(s)

    """
    def __init__(self, name, ra, dec, atrous_do=False, mscale_field_do=False,
    	cal_imsize=512, solint_p=1, solint_a=30, dynamic_range='LD',
    	region_selfcal='empty', region_field='empty', peel_skymodel='empty',
    	outlier_do=False, factor_working_dir='', cal_size_deg=None):
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
        self.dynamic_range = dynamic_range
        if self.dynamic_range.lower() not in ['ld', 'hd']:
            self.log.error('Dynamic range is "{}" but must be either "LD" or "HD".'.format(self.dynamic_range))
            sys.exit(1)
        if region_selfcal.lower() == 'empty':
            # Set to empty list (casa format)
            self.region_selfcal = '[]'
        else:
            if not os.path.exists(region_selfcal):
                self.log.error('Calibrator clean-mask region file {} not found.'.format(region_selfcal))
                sys.exit(1)
            self.region_selfcal = '["{0}"]'.format(region_selfcal)
            self.log.info('Using calibrator clean-mask region file {}'.format(self.region_selfcal))
        self.region_field = region_field
        if self.region_field.lower() == 'empty':
            self.region_field = '[]'
        elif not os.path.exists(self.region_field):
            self.log.error('Facet region file {} not found.'.format(self.region_field))
            sys.exit(1)
        else:
            self.log.info('Using facet clean-mask region file {}'.format(self.region_field))
        self.peel_skymodel = peel_skymodel
        if self.peel_skymodel.lower() == 'empty':
            self.peel_skymodel = None
        elif not os.path.exists(self.peel_skymodel):
            self.log.error('Peel sky model file {} not found.'.format(self.peel_skymodel))
            sys.exit(1)
        else:
            self.log.info('Using sky model file {} for selfcal/peeling'.format(self.peel_skymodel))
        self.is_outlier = outlier_do
        self.cal_size_deg = cal_size_deg

        # Initialize some parameters to default/initial values
        self.selfcal_ok = False # whether selfcal succeeded
        self.wsclean_nchannels = 1 # set number of wide-band channels
        self.use_new_sub_data = False # set flag that tells which subtracted-data column to use
        self.cellsize_verify_deg = 0.00833 # verify subtract cell size
        self.target_rms_rad = 0.2 # preaverage target rms
        self.subtracted_data_colname = 'SUBTRACTED_DATA_ALL' # name of empty data column
        self.pre_average = False # whether to use baseline averaging
        self.blavg_weight_column = 'WEIGHT_SPECTRUM' # name of weights column
        self.peel_calibrator = False # whether to peel calibrator before imaging
        self.solve_all_correlations = False # whether to solve for all corrs for slow gain
        self.do_reset = False # whether to reset this direction
        self.is_patch = False # whether direction is just a patch (not full facet)
        self.skymodel = None # direction's sky model
        self.use_existing_data = False # whether to use existing data for reimaging
        self.full_res_facetimage_freqstep = None # frequency step of existing data
        self.full_res_facetimage_timestep = None # time step of existing data
        self.average_image_data = False # whether to average the existing data before imaging them
        self.facet_imsize = None # size of facet image (None for patch and field directions)
        self.started_operations = []
        self.completed_operations = []
        self.reset_operations = []
        self.cleanup_mapfiles = []

        # Define some directories and files
        self.working_dir = factor_working_dir
        self.save_file = os.path.join(self.working_dir, 'state',
            self.name+'_save.pkl')
        self.vertices_file = self.save_file


    def set_cal_size(self, selfcal_cellsize_arcsec):
        """
        Sets the calibrator image size from the calibrator size (or vice versa)

        Parameters
        ----------
        selfcal_cellsize_arcsec : float
            Cellsize for selfcal imaging
        padding : float, optional
            Padding factor for image. Padded regions to

        """
        self.cellsize_selfcal_deg = selfcal_cellsize_arcsec / 3600.0

        if self.cal_size_deg is None:
            # Set calibrator size from cal_imsize assuming 50% padding
            if self.cal_imsize == 0:
                self.log.error('The cal_imsize must be specified in the directions '
                    'file if cal_size_deg is not specified')
                sys.exit(1)
            else:
                self.cal_size_deg = self.cal_imsize * self.cellsize_selfcal_deg / 1.5
        else:
            if self.cal_imsize == 0:
                # Set image size to size of calibrator, padded to 40% extra
                # (the padded region is not cleaned)
                self.cal_imsize = max(512, self.get_optimum_size(self.cal_size_deg
                    / self.cellsize_selfcal_deg * 1.0 / 0.6))

        self.cal_radius_deg = self.cal_size_deg / 2.0
        self.cal_rms_box = self.cal_size_deg / self.cellsize_selfcal_deg


    def set_imcal_parameters(self, parset, bands, facet_cellsize_arcsec=None,
        facet_robust=None, facet_taper_arcsec=None, facet_min_uv_lambda=None,
        imaging_only=False, use_existing_data=False, existing_data_freqstep=None,
        existing_data_timestep=None):
        """
        Sets various parameters for imaging and calibration

        Parameters
        ----------
        parset : dict
            Parset of operation
        bands : list of Band objects
            Bands for this operation
        facet_cellsize_arcsec : float, optional
            Pixel size in arcsec for facet imaging
        facet_robust : float, optional
            Briggs robust parameter for facet imaging
        facet_taper_arcsec : float, optional
            Taper in arcsec for facet imaging
        imaging_only : bool, optional
            If True, set only imaging-related parameters
        use_existing_data : bool, optional
            If True, existing, potentially averaged, data are to be used instead
            of the unaveraged data
        existing_data_freqstep : int, optional
            The freqstep used to average the existing data
        existing_data_timestep : int, optional
            The timestep used to average the existing data

        """
        mean_freq_mhz = np.mean([b.freq for b in bands]) / 1e6
        min_peak_smearing_factor = 1.0 - parset['imaging_specific']['max_peak_smearing']
        padding = parset['imaging_specific']['wsclean_image_padding']
        self.wsclean_patch_model_padding = parset['imaging_specific']['wsclean_patch_model_padding']
        if parset['imaging_specific']['skip_facet_imaging']:
            self.wsclean_facet_model_padding = self.wsclean_patch_model_padding
        else:
            self.wsclean_facet_model_padding = parset['imaging_specific']['wsclean_facet_model_padding']
        wsclean_nchannels_factor = parset['imaging_specific']['wsclean_nchannels_factor']
        chan_width_hz = bands[0].chan_width_hz
        nchan = bands[0].nchan
        timestep_sec = bands[0].timepersample
        ntimes = bands[0].minSamplesPerFile
        if (use_existing_data and existing_data_freqstep is not None and
            existing_data_timestep is not None):
            # Adjust the above values to match the existing data if necessary
            chan_width_hz *= existing_data_freqstep
            nchan /= existing_data_freqstep
            timestep_sec *= existing_data_timestep
            ntimes /= existing_data_timestep

        nbands = self.get_nbands(bands)
        preaverage_flux_jy = parset['calibration_specific']['preaverage_flux_jy']
        tec_block_mhz = parset['calibration_specific']['tec_block_mhz']
        peel_flux_jy = parset['calibration_specific']['peel_flux_jy']

        self.robust_selfcal = parset['imaging_specific']['selfcal_robust']
        self.solve_min_uv_lambda = parset['calibration_specific']['solve_min_uv_lambda']
        self.selfcal_min_uv_lambda = parset['imaging_specific']['selfcal_min_uv_lambda']
        self.use_selfcal_clean_threshold = parset['imaging_specific']['selfcal_clean_threshold']
        self.use_selfcal_adaptive_threshold = parset['imaging_specific']['selfcal_adaptive_threshold']

        if facet_cellsize_arcsec is None:
            facet_cellsize_arcsec = parset['imaging_specific']['selfcal_cellsize_arcsec']
        self.cellsize_facet_deg = facet_cellsize_arcsec / 3600.0

        if facet_robust is None:
            facet_robust = parset['imaging_specific']['selfcal_robust']
        self.robust_facet = facet_robust

        if facet_taper_arcsec is None:
            facet_taper_arcsec = 0.0
        self.taper_facet_arcsec = facet_taper_arcsec

        if facet_min_uv_lambda is None:
            facet_min_uv_lambda = parset['imaging_specific']['selfcal_min_uv_lambda']
        self.facet_min_uv_lambda = facet_min_uv_lambda

        self.set_imaging_parameters(nbands, padding)
        self.set_averaging_steps_and_solution_intervals(chan_width_hz, nchan,
            timestep_sec, ntimes, nbands, mean_freq_mhz, self.skymodel,
            preaverage_flux_jy, min_peak_smearing_factor, tec_block_mhz,
            peel_flux_jy, imaging_only=imaging_only)

        # Set channelsout for wide-band imaging with WSClean. Note that the
        # number of WSClean channels must be an even divisor of the total number
        # of channels in the full bandwidth after averaging to prevent
        # mismatches during the predict step on the unaveraged data. For selfcal,
        # we use nchannels = nbands and fit a third-order polynomial
        # using 3 averaged channels
        #
        # Also define the image suffixes (which depend on whether or not
        # wide-band clean is done)
        if self.use_wideband:
            self.wsclean_nchannels = max(1, int(np.ceil(nbands / float(wsclean_nchannels_factor))))
            nchan_after_avg = nchan * nbands / self.facetimage_freqstep
            self.nband_pad = 0 # padding to allow self.wsclean_nchannels to be a divisor
            if parset['imaging_specific']['wsclean_add_bands']:
                while nchan_after_avg % self.wsclean_nchannels:
                    self.nband_pad += 1
                    nchan_after_avg = nchan * (nbands + self.nband_pad) / self.facetimage_freqstep
            else:
                while nchan_after_avg % self.wsclean_nchannels:
                    self.wsclean_nchannels += 1
            if self.wsclean_nchannels > nbands:
                self.wsclean_nchannels = nbands

            self.wsclean_nchannels_selfcal = nbands
            self.wsclean_suffix = '-MFS-image.fits'
        else:
            self.wsclean_nchannels = 1
            self.wsclean_nchannels_selfcal = 1
            self.nband_pad = 0
            self.wsclean_suffix = '-image.fits'

        # Set the baseline-averaging limit for WSClean, which depends on the
        # integration time given the specified maximum allowed smearing. We scale
        # it from the imaging cell size assuming normal sampling as:
        #
        # max baseline in nwavelengths = 1 / theta_rad ~= 1 / (cellsize_deg * 3 * pi / 180)
        #
        # nwavelengths = max baseline in nwavelengths * 2 * pi * integration time in seconds / (24 * 60 * 60)
        #
        if (hasattr(self, 'facetselfcal_timestep_sec') and
            parset['imaging_specific']['wsclean_bl_averaging']):
            max_baseline = 1 / (3 * self.cellsize_selfcal_deg * np.pi / 180)
            self.facetselfcal_wsclean_nwavelengths = int(max_baseline * 2*np.pi * self.facetselfcal_timestep_sec / (24*60*60))
        else:
            self.facetselfcal_wsclean_nwavelengths = 0
        if (hasattr(self, 'facetimage_timestep_sec') and
            parset['imaging_specific']['wsclean_bl_averaging']):
            max_baseline = 1 / (3 * self.cellsize_facet_deg * np.pi / 180)
            self.facetimage_wsclean_nwavelengths = int(max_baseline * 2*np.pi * self.facetimage_timestep_sec / (24*60*60))
        else:
            self.facetimage_wsclean_nwavelengths = 0


    def set_imaging_parameters(self, nbands, padding=1.05):
        """
        Sets various parameters for images in facetselfcal and facetimage pipelines

        Parameters
        ----------
        nbands : int
            Number of bands
        padding : float, optional
            Padding factor by which size of facet is multiplied to determine
            the facet image size

        """
        # Set facet image size
        if hasattr(self, 'width'):
            self.facet_imsize_nopadding = int(self.width / self.cellsize_facet_deg)
            self.facet_imsize = max(512, self.get_optimum_size(self.width
                / self.cellsize_facet_deg * padding))

        # Determine whether the total bandwidth is large enough that wide-band
        # imaging is needed
        if nbands > 5:
            self.use_wideband = True
        else:
            self.use_wideband = False

        # Set number of iterations and threshold for full facet image, scaled to
        # the number of bands. We use 6 times more iterations for the full2
        # image to ensure the imager has a reasonable chance to reach the
        # threshold first (which is set by the masking step)
        scaling_factor = np.sqrt(np.float(nbands))
        self.wsclean_full1_image_niter = int(2000 * scaling_factor)
        self.wsclean_full1_image_threshold_jy =  1.5e-3 * 0.7 / scaling_factor
        self.wsclean_full2_image_niter = int(12000 * scaling_factor)

        # Set multiscale imaging mode for facet imaging: Get source sizes and
        # check for large sources (anything above 4 arcmin -- the CC sky model
        # was convolved with a Gaussian of 1 arcmin, so unresolved sources have
        # sizes of ~ 1 arcmin)
        large_size_arcmin = 6.0
        if self.mscale_field_do is None:
            sizes_arcmin = self.get_source_sizes()
            if any([s > large_size_arcmin for s in sizes_arcmin]):
                self.mscale_field_do = True
            else:
                self.mscale_field_do = False
        if self.mscale_field_do:
            self.wsclean_multiscale = '-multiscale,'
            self.wsclean_full1_image_niter /= 2 # fewer iterations are needed
            self.wsclean_full2_image_niter /= 2 # fewer iterations are needed
        else:
            self.wsclean_multiscale = ''

        # Set whether to use wavelet module in calibrator masking
        if self.atrous_do is None:
            sizes_arcmin = self.get_source_sizes(cal_only=True)
            if any([s > large_size_arcmin for s in sizes_arcmin]):
                self.atrous_do = True
            else:
                self.atrous_do = False

        # If wavelet module is activated, also activate multiscale clean
        if self.atrous_do:
            self.wsclean_selfcal_multiscale = True
        else:
            self.wsclean_selfcal_multiscale = False


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


    def set_skymodel(self, skymodel):
        """
        Sets the direction sky model

        The sky model is filtered to include only those sources within the
        direction facet

        Parameters
        ----------
        skymodel : LSMTool SkyModel object
            CC sky model

        """
        x, y, midRA, midDec = skymodel._getXY()
        xv, yv = radec2xy(self.vertices[0], self.vertices[1], midRA, midDec)
        xyvertices = np.array([[xp, yp] for xp, yp in zip(xv, yv)])
        bbPath = mplPath.Path(xyvertices)
        inside = np.zeros(len(skymodel), dtype=bool)
        for i in range(len(skymodel)):
            inside[i] = bbPath.contains_point((x[i], y[i]))
        skymodel.select(inside, force=True)
        self.skymodel = skymodel


    def get_source_sizes(self, cal_only=False):
        """
        Returns list of source sizes in arcmin
        """
        skymodel = self.skymodel.copy()
        if cal_only:
            dist = skymodel.getDistance(self.ra, self.dec, byPatch=True)
            skymodel.select(dist < self.cal_radius_deg, aggregate=True)
        sizes = skymodel.getPatchSizes(units='arcmin', weight=False)

        return sizes


    def get_cal_fluxes(self, fwhmArcsec=25.0, threshold=0.1):
        """
        Returns total flux density in Jy and max peak flux density in
        Jy per beam for calibrator

        Parameters
        ----------
        fwhmArcsec : float, optional
            Smoothing scale
        threshold : float, optional
            Threshold

        Returns
        -------
        tot_flux_jy, peak_flux_jy_bm : float, float
            Total flux density in Jy and max peak flux density in
            Jy per beam for calibrator

        """
        dist = self.skymodel.getDistance(self.ra, self.dec)
        skymodel = self.skymodel.copy()
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
        timestep_sec, ntimes_min, nbands, mean_freq_mhz, initial_skymodel=None,
        preaverage_flux_jy=0.0, min_peak_smearing_factor=0.95, tec_block_mhz=10.0,
        peel_flux_jy=25.0, imaging_only=False):
        """
        Sets the averaging step sizes and solution intervals

        The averaging is set by the need to keep product of bandwidth and time
        smearing below 1 - min_peak_smearing_factor. Averaging is limited to
        a maximum of ~ 120 sec and ~ 2 MHz.

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
            Channel width in Hz
        nchan : int
            Number of channels per band
        timestep_sec : float
            Time step
        ntimes_min : int
            Minimum number of timeslots in a chunk
        nbands : int
            Number of bands
        mean_freq_mhz : float
            Mean frequency in MHz of full bandwidth
        initial_skymodel : LSMTool SkyModel object, optional
            Sky model used to check source sizes
        preaverage_flux_jy : bool, optional
            Use baseline-dependent averaging and solint_time_p = 1 for phase-only
            calibration for sources below this flux value
        min_peak_smearing_factor : float, optional
            Min allowed peak flux density reduction due to smearing at the mean
            frequency (facet imaging only)
        tec_block_mhz : float, optional
            Size of frequency block in MHz over which a single TEC solution is
            fit
        peel_flux_jy : float, optional
            Peel cailbrators with fluxes above this value
        imaging_only : bool, optional
            If True, set only imaging-related parameters

        """
        # generate a (numpy-)array with the divisors of nchan
        tmp_divisors = []
        for step in range(nchan, 0, -1):
            if (nchan % step) == 0:
                tmp_divisors.append(step)
        freq_divisors = np.array(tmp_divisors)

        # For selfcal, use the size of the calibrator to set the averaging
        # steps
        if not imaging_only and self.cal_size_deg is not None:
            # Set min allowable smearing reduction factor for bandwidth and time
            # smearing so that they are equal and their product is 0.85
            min_peak_smearing_factor_selfcal = 0.85
            peak_smearing_factor = np.sqrt(min_peak_smearing_factor_selfcal)

            # Get target time and frequency averaging steps
            delta_theta_deg = self.cal_size_deg / 2.0
            self.log.debug('Calibrator is {} deg across'.format(delta_theta_deg*2.0))
            resolution_deg = 3.0 * self.cellsize_selfcal_deg # assume normal sampling of restoring beam
            target_timewidth_s = min(120.0, self.get_target_timewidth(delta_theta_deg,
                resolution_deg, peak_smearing_factor))
            if self.dynamic_range.lower() == 'hd':
                # For high-dynamic range calibration, we use 0.2 MHz per channel
                target_bandwidth_mhz = 0.2
            else:
                target_bandwidth_mhz = min(2.0, self.get_target_bandwidth(mean_freq_mhz,
                    delta_theta_deg, resolution_deg, peak_smearing_factor))
            self.log.debug('Target timewidth for selfcal is {} s'.format(target_timewidth_s))
            self.log.debug('Target bandwidth for selfcal is {} MHz'.format(target_bandwidth_mhz))

            # Find averaging steps for given target values
            self.facetselfcal_freqstep = max(1, min(int(round(target_bandwidth_mhz * 1e6 / chan_width_hz)), nchan))
            self.facetselfcal_freqstep = freq_divisors[np.argmin(np.abs(freq_divisors - self.facetselfcal_freqstep))]
            self.facetselfcal_timestep = max(1, int(round(target_timewidth_s / timestep_sec)))
            self.facetselfcal_timestep_sec = self.facetselfcal_timestep * timestep_sec
            self.log.debug('Using averaging steps of {0} channels and {1} time slots '
                'for selfcal'.format(self.facetselfcal_freqstep, self.facetselfcal_timestep))

            # For selfcal verify, average to 2 MHz per channel and 120 sec per time
            # slot
            self.verify_freqstep = max(1, min(int(round(2.0 * 1e6 / chan_width_hz)), nchan))
            self.verify_freqstep = freq_divisors[np.argmin(np.abs(freq_divisors - self.verify_freqstep))]
            self.verify_timestep = max(1, int(round(120.0 / timestep_sec)))

        # For facet imaging, use the facet image size (before padding) to set the averaging steps
        if self.facet_imsize is not None:
            # Set min allowable smearing reduction factor for bandwidth and time
            # smearing so that they are equal and their product is
            # min_peak_smearing_factor
            peak_smearing_factor = np.sqrt(min_peak_smearing_factor)

            # Get target time and frequency averaging steps
            delta_theta_deg = self.facet_imsize_nopadding * self.cellsize_facet_deg / 2.0
            self.log.debug('Facet image before padding is {0} x {0} pixels ({1} x {1} deg)'.format(
                self.facet_imsize_nopadding, delta_theta_deg*2.0))
            resolution_deg = 3.0 * self.cellsize_facet_deg # assume normal sampling of restoring beam
            target_timewidth_s = min(120.0, self.get_target_timewidth(delta_theta_deg,
                resolution_deg, peak_smearing_factor))
            target_bandwidth_mhz = min(2.0, self.get_target_bandwidth(mean_freq_mhz,
                delta_theta_deg, resolution_deg, peak_smearing_factor))
            self.log.debug('Target timewidth for facet imaging is {} s'.format(target_timewidth_s))
            self.log.debug('Target bandwidth for facet imaging is {} MHz'.format(target_bandwidth_mhz))

            # Find averaging steps for given target values
            self.facetimage_freqstep = max(1, min(int(round(target_bandwidth_mhz * 1e6 / chan_width_hz)), nchan))
            self.facetimage_freqstep = freq_divisors[np.argmin(np.abs(freq_divisors - self.facetimage_freqstep))]
            self.facetimage_timestep = max(1, int(round(target_timewidth_s / timestep_sec)))
            self.facetimage_timestep_sec = self.facetimage_timestep * timestep_sec
            self.log.debug('Using averaging steps of {0} channels and {1} time slots '
                'for facet imaging'.format(self.facetimage_freqstep, self.facetimage_timestep))

        # Set time intervals for selfcal solve steps
        if not imaging_only:
            # Calculate the effective flux density. This is the one used to set the
            # intervals. It is the peak flux density adjusted to account for cases
            # in which the total flux density is larger than the peak flux density
            # would indicate (either due to source being extended or to multiple
            # calibrator sources). In these cases, we can use a higher effective
            # flux density to set the intervals. A scaling with a power of 1/1.5
            # seems to work well
            total_flux_jy, peak_flux_jy_bm = self.get_cal_fluxes()
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

            # Set fast (phase-only) solution time interval
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
                    target_timewidth_s = 16.0 * (ref_flux_jy / effective_flux_jy)**2
                    if target_timewidth_s > 16.0:
                        target_timewidth_s = 16.0
                    self.solint_time_p = max(1, int(round(target_timewidth_s / timestep_sec)))

            # Set slow (gain) solution time interval
            if self.solint_time_a == 0:
                # Slow gain solve is per band, so don't scale the interval with
                # the number of bands but only with the effective flux. Also,
                # avoid cases in which the last solution interval is much smaller
                # than the target interval (for the smallest time chunk, this
                # assumes that the chunks are all about the same length)
                ref_flux_jy = 1.4
                target_timewidth_s = 1200.0 * (ref_flux_jy / effective_flux_jy)**2
                if target_timewidth_s < 240.0:
                    target_timewidth_s = 240.0
                if target_timewidth_s > 1200.0:
                    target_timewidth_s = 1200.0
                solint_time_a = int(round(target_timewidth_s / timestep_sec))
                solint_time_a_lower = solint_time_a
                solint_time_a_upper = solint_time_a
                while (ntimes_min % solint_time_a_lower > 0 and
                    ntimes_min % solint_time_a_lower < solint_time_a_lower / 2.0):
                    solint_time_a_lower -= 1
                while (ntimes_min % solint_time_a_upper > 0 and
                    ntimes_min % solint_time_a_upper < solint_time_a_upper / 2.0):
                    solint_time_a_upper += 1
                if solint_time_a - solint_time_a_lower <= solint_time_a_upper - solint_time_a:
                    self.solint_time_a = solint_time_a_lower
                else:
                    self.solint_time_a = solint_time_a_upper

            self.log.debug('Using solution intervals of {0} (fast) and {1} '
                '(slow) time slots'.format(self.solint_time_p, self.solint_time_a))

            # Set frequency intervals for selfcal solve steps. The interval for
            # slow (amp) selfcal should be the number of channels in a band after
            # averaging.
            num_chan_per_band_after_avg = nchan / self.facetselfcal_freqstep
            if self.dynamic_range.lower() == 'hd':
                # For high-dynamic range solve, the interval for slow (amp) selfcal
                # should be 1 (every channel)
                self.solint_freq_a = 1
            else:
                self.solint_freq_a = num_chan_per_band_after_avg

            # The interval for fast (phase) selfcal should be the number of
            # channels in tec_block_mhz, but no less than 2 MHz
            min_block_mhz = 2.0
            if tec_block_mhz < min_block_mhz:
                self.log.warn('Setting TEC block size to minimum allowed value of {} MHz'.format(min_block_mhz))
                tec_block_mhz = min_block_mhz
            mhz_per_chan_after_avg = self.facetselfcal_freqstep * chan_width_hz / 1e6
            total_bandwidth_mhz = nchan * nbands * chan_width_hz / 1e6
            num_cal_blocks = np.ceil(total_bandwidth_mhz / tec_block_mhz)
            nchan_per_block = np.ceil(num_chan_per_band_after_avg * nbands /
                num_cal_blocks)

            # Check for a partial block, and adjust the number to ensure that
            # it is at least half of the desired block size
            num_cal_blocks_lower = num_cal_blocks
            num_cal_blocks_upper = num_cal_blocks
            partial_block_mhz = self.calc_partial_block(num_chan_per_band_after_avg,
                nbands, num_cal_blocks, mhz_per_chan_after_avg)
            while (partial_block_mhz > 0.0 and partial_block_mhz < tec_block_mhz/2.0):
                num_cal_blocks_lower -= 1
                partial_block_mhz = self.calc_partial_block(num_chan_per_band_after_avg,
                    nbands, num_cal_blocks_lower, mhz_per_chan_after_avg)
            partial_block_mhz = self.calc_partial_block(num_chan_per_band_after_avg,
                nbands, num_cal_blocks, mhz_per_chan_after_avg)
            while (partial_block_mhz > 0.0 and partial_block_mhz < tec_block_mhz/2.0):
                num_cal_blocks_upper += 1
                partial_block_mhz = self.calc_partial_block(num_chan_per_band_after_avg,
                    nbands, num_cal_blocks_upper, mhz_per_chan_after_avg)
            if num_cal_blocks - num_cal_blocks_lower < num_cal_blocks_upper - num_cal_blocks:
                num_cal_blocks = num_cal_blocks_lower
            else:
                num_cal_blocks = num_cal_blocks_upper
            if num_cal_blocks < 1:
                num_cal_blocks = 1
            self.num_cal_blocks = num_cal_blocks
            self.num_bands_per_cal_block = int(np.ceil(nbands / float(num_cal_blocks)))
            self.solint_freq_p = int(np.ceil(num_chan_per_band_after_avg * nbands /
                float(num_cal_blocks)))

        # Set name of column to use for data and averaged weights
        if self.pre_average:
            self.data_column = 'BLAVG_DATA'
            self.blavg_weight_column = 'BLAVG_WEIGHT_SPECTRUM'
        else:
            self.data_column = 'DATA'
            self.blavg_weight_column = 'WEIGHT_SPECTRUM'


    def calc_partial_block(self, num_chan_per_band_after_avg, nbands,
        num_cal_blocks, mhz_per_chan_after_avg):
        """
        Returns size of partial block in MHz
        """
        nchan_per_block = np.ceil(num_chan_per_band_after_avg * nbands /
            num_cal_blocks)
        partial_block_mhz = (num_chan_per_band_after_avg * nbands %
            nchan_per_block) * mhz_per_chan_after_avg
        return partial_block_mhz


    def get_target_timewidth(self, delta_theta, resolution, reduction_factor):
        """
        Returns the time width for given peak flux density reduction factor

        Parameters
        ----------
        delta_theta : float
            Distance from phase center
        resolution : float
            Resolution of restoring beam
        reduction_factor : float
            Ratio of pre-to-post averaging peak flux density

        Returns
        -------
        delta_time : float
            Time width in seconds for target reduction_factor

        """
        delta_time = np.sqrt( (1.0 - reduction_factor) /
            (1.22E-9 * (delta_theta / resolution)**2.0) )

        return delta_time


    def get_bandwidth_smearing_factor(self, freq, delta_freq, delta_theta, resolution):
        """
        Returns peak flux density reduction factor due to bandwidth smearing

        Parameters
        ----------
        freq : float
            Frequency at which averaging will be done
        delta_freq : float
            Bandwidth over which averaging will be done
        delta_theta : float
            Distance from phase center
        resolution : float
            Resolution of restoring beam

        Returns
        -------
        reduction_factor : float
            Ratio of pre-to-post averaging peak flux density

        """
        beta = (delta_freq/freq) * (delta_theta/resolution)
        gamma = 2*(np.log(2)**0.5)
        reduction_factor = ((np.pi**0.5)/(gamma * beta)) * (erf(beta*gamma/2.0))

        return reduction_factor


    def get_target_bandwidth(self, freq, delta_theta, resolution, reduction_factor):
        """
        Returns the bandwidth for given peak flux density reduction factor

        Parameters
        ----------
        freq : float
            Frequency at which averaging will be done
        delta_theta : float
            Distance from phase center
        resolution : float
            Resolution of restoring beam
        reduction_factor : float
            Ratio of pre-to-post averaging peak flux density

        Returns
        -------
        delta_freq : float
            Bandwidth over which averaging will be done
        """
        # Increase delta_freq until we drop below target reduction_factor
        delta_freq = 1e-3 * freq
        while self.get_bandwidth_smearing_factor(freq, delta_freq, delta_theta,
            resolution) > reduction_factor:
            delta_freq *= 1.1

        return delta_freq


    def find_peel_skymodel(self):
        """
        Searches for an appropriate sky model for peeling
        """
        if self.peel_skymodel is not None:
            return

        max_separation_arcmin = 1.0
        factor_lib_dir = os.path.dirname(os.path.abspath(__file__))
        skymodel_dir = os.path.join(os.path.split(factor_lib_dir)[0], 'skymodels')
        skymodels = glob.glob(os.path.join(skymodel_dir, '*.skymodel'))
        for skymodel in skymodels:
            try:
                s = lsmtool.load(skymodel)
                dist_deg = s.getDistance(self.ra, self.dec)
                if any(dist_deg*60.0 < max_separation_arcmin):
                    self.peel_skymodel = skymodel
                    break
            except IOError:
                pass


    def get_nbands(self, bands):
        """
        Returns total number of bands including missing ones

        Parameters
        ----------
        bands : list of Band objects
            Bands for this operation

        """
        freqs_hz = [b.freq for b in bands]
        chan_width_hz = bands[0].chan_width_hz
        nchan = bands[0].nchan
        freq_width_hz = chan_width_hz * nchan

        # Find gaps, if any
        missing_bands = []
        for i, (freq1, freq2) in enumerate(zip(freqs_hz[:-1], freqs_hz[1:])):
            ngap = int(round((freq2 - freq1)/freq_width_hz))
            missing_bands.extend([i + j + 1 for j in range(ngap-1)])

        return len(bands) + len(missing_bands)


    def save_state(self):
        """
        Saves the direction state to a file

        """
        import pickle

        with open(self.save_file, 'wb') as f:
            # Remove log and skymodel objects, as they cannot be pickled
            save_dict = self.__dict__.copy()
            save_dict.pop('log')
            save_dict.pop('skymodel')
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
                if 'dir_dep_parmdb_mapfile' in d:
                    self.dir_dep_parmdb_mapfile = d['dir_dep_parmdb_mapfile']
                if 'facet_model_mapfile' in d:
                    self.facet_model_mapfile = d['facet_model_mapfile']
                if 'wsclean_modelimg_size_mapfile' in d:
                    self.wsclean_modelimg_size_mapfile = d['wsclean_modelimg_size_mapfile']
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
            List of names of operations to reset. Reset is done only if the
            operation name appears in self.reset_operations. If None, all
            started and completed operations that are in self.reset_operations
            are reset

        """
        if op_names is None:
            op_names = self.completed_operations[:] + self.started_operations[:]
        elif type(op_names) is str:
            op_names = [op_names]
        op_names_reset = [op for op in op_names if op in self.reset_operations]
        if len(op_names_reset) > 0:
            self.log.info('Resetting state for operation(s): {}'.format(', '.join(op_names_reset)))
        else:
            return

        # Reset selfcal flag
        if 'facetselfcal' in op_names_reset:
            self.selfcal_ok = False

        # Remove operation name from lists of started and completed operations
        # and delete the results directories
        for op_name in op_names_reset:
            while op_name.lower() in self.completed_operations:
                self.completed_operations.remove(op_name.lower())
            while op_name.lower() in self.started_operations:
                self.started_operations.remove(op_name.lower())

            # Delete results directory for this operation
            op_dir = os.path.join(self.working_dir, 'results', op_name.lower(), self.name)
            if os.path.exists(op_dir):
                os.system('rm -rf {0}'.format(op_dir))

        self.save_state()


    def cleanup(self):
        """
        Cleans up unneeded data
        """
        from lofarpipe.support.data_map import DataMap
        import glob

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

                        # Deal with special case of f being a WSClean image
                        if f.endswith('MFS-image.fits'):
                            # Search for related images and delete if found
                            image_root = f.split('MFS-image.fits')[0]
                            extra_files = glob.glob(image_root+'*.fits')
                            for e in extra_files:
                                if os.path.exists(e):
                                    os.system('rm -rf {0}'.format(e))
                        elif f.endswith('-image.fits'):
                            # Search for related images and delete if found
                            image_root = f.split('-image.fits')[0]
                            extra_files = glob.glob(image_root+'*.fits')
                            for e in extra_files:
                                if os.path.exists(e):
                                    os.system('rm -rf {0}'.format(e))
            except IOError:
                pass
