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
from scipy.special import erf
import sys


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
    make_final_image : bool, optional
        Make final image of this direction, after all directions have been
        selfcaled?
    cal_size_deg : float, optional
        Size in degrees of calibrator source(s)
    cal_flux_jy : float, optional
        Apparent flux density in Jy of calibrator source

    """
    def __init__(self, name, ra, dec, atrous_do=False, mscale_field_do=False,
    	cal_imsize=512, solint_p=1, solint_a=30, dynamic_range='LD',
    	region_selfcal='empty', region_field='empty', peel_skymodel='empty',
    	outlier_do=False, factor_working_dir='', make_final_image=False,
    	cal_size_deg=None, cal_flux_jy=None):

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
        self.dynamic_range = dynamic_range
        if self.dynamic_range.lower() not in ['ld', 'hd']:
            self.log.error('Dynamic range is "{}" but must be either "LD" or "HD".'.format(self.dynamic_range))
            sys.exit(1)
        if region_selfcal.lower() == 'empty':
            # Set to empty list (casapy format)
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
            self.log.info('Using sky model file {} for selfcal'.format(self.peel_skymodel))
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
        self.wsclean_nchannels = 1 # set number of wide-band channels
        self.use_new_sub_data = False # set flag that tells which subtracted-data column to use
        self.cellsize_selfcal_deg = 0.000417 # selfcal cell size
        self.cellsize_facet_deg = 0.000417 # facet image cell size
        self.cellsize_verify_deg = 0.00833 # verify subtract cell size
        self.target_rms_rad = 0.2 # preaverage target rms
        self.subtracted_data_colname = 'SUBTRACTED_DATA_ALL' # name of empty data column
        self.pre_average = False # whether to use baseline averaging
        self.blavg_weight_column = 'WEIGHT_SPECTRUM' # name of weights column
        self.peel_calibrator = False # whether to peel calibrator before imaging
        self.started_operations = []
        self.completed_operations = []
        self.cleanup_mapfiles = []
        self.do_reset = False # whether to reset this direction
        self.is_patch = False # whether direction is just a patch (not full facet)
        self.nchunks = 1
        self.num_selfcal_groups = 1
        self.timeSlotsPerParmUpdate = 100

        # Set the size of the calibrator
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
                    / self.cellsize_selfcal_deg * 1.2)) # cal imsize has 20% padding

        self.cal_radius_deg = self.cal_size_deg / 2.0
        self.cal_rms_box = self.cal_size_deg / self.cellsize_selfcal_deg

        # Define some directories and files
        self.working_dir = factor_working_dir
        self.save_file = os.path.join(self.working_dir, 'state',
            self.name+'_save.pkl')
        self.vertices_file = self.save_file


    def set_imcal_parameters(self, nbands_per_channel, chan_width_hz,
    	nchan, timestep_sec, ntimes, nbands, mean_freq_mhz, initial_skymodel=None,
    	preaverage_flux_jy=0.0, min_peak_smearing_factor=0.95, peel_flux_jy=25.0):
        """
        Sets various parameters for imaging and calibration

        Parameters
        ----------
        nbands_per_channel : int
            Number of bands per output channel (WSClean only)
        nchan_per_band : int
            Number of channels per band
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
        peel_flux_jy : float, optional
            Peel cailbrators with fluxes above this value

        """
        self.set_imaging_parameters(nbands, nbands_per_channel, nchan,
            initial_skymodel)
        self.set_averaging_steps_and_solution_intervals(chan_width_hz, nchan,
            timestep_sec, ntimes, nbands, mean_freq_mhz, initial_skymodel,
            preaverage_flux_jy, min_peak_smearing_factor, peel_flux_jy)


    def set_imaging_parameters(self, nbands, nbands_per_channel, nchan_per_band,
        initial_skymodel=None, padding=1.05):
        """
        Sets various parameters for images in facetselfcal and facetimage pipelines

        Parameters
        ----------
        nbands : int
            Number of bands
        nbands_per_channel : int
            Number of bands per output channel (WSClean only)
        nchan_per_band : int
            Number of channels per band
        initial_skymodel : LSMTool SkyModel object, optional
            Sky model used to check source sizes
        padding : float, optional
            Padding factor by which size of facet is multiplied to determine
            the facet image size

        """
        # Set facet image size
        if hasattr(self, 'width'):
            self.facet_imsize = max(512, self.get_optimum_size(self.width
                / self.cellsize_selfcal_deg * padding))
        else:
            self.facet_imsize = None

        self.cal_wplanes = self.set_wplanes(self.cal_imsize)
        self.facet_wplanes = self.set_wplanes(self.facet_imsize)

        # Determine whether the total bandwidth is large enough that wide-band
        # imaging is needed
        if nbands > 5:
            self.use_wideband = True
        else:
            self.use_wideband = False

        # Set number of channels for wide-band imaging with WSClean and nterms
        # for the CASA imager. Note that the number of WSClean channels must be
        # an even divisor of the total number of channels in the full bandwidth.
        # For now, we just set the number of channels to the number of bands to
        # avoid any issues with this setting (revisit once we switch to fitting
        # a spectral function during deconvolution).
        #
        # Also define the image suffixes (which depend on whether or not
        # wide-band clean is done)
        if self.use_wideband:
            self.wsclean_nchannels = nbands
            self.nterms = 2
            self.casa_suffix = '.tt0'
            self.wsclean_suffix = '-MFS-image.fits'
        else:
            self.wsclean_nchannels = 1
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
        # sources (anything above 4 arcmin -- the CC sky model was convolved
        # with a Gaussian of 1 arcmin, so unresolved sources have sizes of ~
        # 1 arcmin)
        if initial_skymodel is not None and self.mscale_field_do is None:
            sizes = self.get_source_sizes(initial_skymodel.copy())
            large_size_arcmin = 4.0
            if any([s > large_size_arcmin for s in sizes]):
                self.mscale_field_do = True
            else:
                self.mscale_field_do = False
        if self.mscale_field_do:
            self.casa_multiscale = '[0, 3, 7, 25, 60, 150]'
            self.wsclean_multiscale = '-multiscale,'
            self.wsclean_full_image_niter /= 2.0 # fewer iterations are needed
        else:
            self.casa_multiscale = '[0]'
            self.wsclean_multiscale = ''

        # Set wavelet source-finding mode
        if self.atrous_do is None:
            if self.mscale_field_do:
                self.atrous_do = True
            else:
                self.atrous_do = False


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
        timestep_sec, ntimes_min, nbands, mean_freq_mhz, initial_skymodel=None,
        preaverage_flux_jy=0.0, min_peak_smearing_factor=0.95, peel_flux_jy=25.0):
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
            Minimum number of timeslots per band, currently not used
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
        peel_flux_jy : float, optional
            Peel cailbrators with fluxes above this value

        """
        # generate a (numpy-)array with the divisors of nchan
        tmp_divisors = []
        for step in range(nchan, 0, -1):
            if (nchan % step) == 0:
                tmp_divisors.append(step)
        freq_divisors = np.array(tmp_divisors)

        # For selfcal, use the size of the calibrator to set the averaging
        # steps
        if self.cal_size_deg is not None:
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
            self.log.debug('Using averaging steps of {0} channels and {1} time slots '
                'for selfcal'.format(self.facetselfcal_freqstep, self.facetselfcal_timestep))

        # For facet imaging, use the facet image size to set the averaging steps
        if self.facet_imsize is not None:
            # Set min allowable smearing reduction factor for bandwidth and time
            # smearing so that they are equal and their product is
            # min_peak_smearing_factor
            peak_smearing_factor = np.sqrt(min_peak_smearing_factor)

            # Get target time and frequency averaging steps
            delta_theta_deg = self.facet_imsize * self.cellsize_facet_deg / 2.0
            self.log.debug('Facet image is {0} x {0} pixels ({1} x {1} deg)'.format(
                self.facet_imsize, delta_theta_deg*2.0))
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
            self.log.debug('Using averaging steps of {0} channels and {1} time slots '
                'for facet imaging'.format(self.facetimage_freqstep, self.facetimage_timestep))

        # For selfcal verify, average to 2 MHz per channel and 120 sec per time
        # slot
        self.verify_freqstep = max(1, min(int(round(2.0 * 1e6 / chan_width_hz)), nchan))
        self.verify_freqstep = freq_divisors[np.argmin(np.abs(freq_divisors - self.verify_freqstep))]
        self.verify_timestep = max(1, int(round(120.0 / timestep_sec)))

        # Set timeSlotsPerParmUpdate to an even divisor of the number of time slots
        # to work around a bug in DPPP ApplyCal
        self.timeSlotsPerParmUpdate = 100
        if ntimes_min<self.timeSlotsPerParmUpdate:
            self.timeSlotsPerParmUpdate = ntimes_min
        else:
            while ntimes_min % self.timeSlotsPerParmUpdate:
                self.timeSlotsPerParmUpdate += 1

        # Set time intervals for selfcal solve steps. The initial skymodel is
        # None for field-type directions, so the steps below are skipped
        if initial_skymodel is not None:
            # Calculate the effective flux density. This is the one used to set the
            # intervals. It is the peak flux density adjusted to account for cases
            # in which the total flux density is larger than the peak flux density
            # would indicate (either due to source being extended or to multiple
            # calibrator sources). In these cases, we can use a higher effective
            # flux density to set the intervals. A scaling with a power of 1/1.5
            # seems to work well
            total_flux_jy, peak_flux_jy_bm = self.get_cal_fluxes(initial_skymodel.copy())
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

            # Set peeling flag
            if effective_flux_jy < peel_flux_jy:
                self.peel_calibrator = True
            else:
                self.peel_calibrator = False

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
                # the number of bands but only with the effective flux
                ref_flux_jy = 1.4
                target_timewidth_s = 1200.0 * (ref_flux_jy / effective_flux_jy)**2
                if target_timewidth_s < 240.0:
                    target_timewidth_s = 240.0
                if target_timewidth_s > 1200.0:
                    target_timewidth_s = 1200.0
                self.solint_time_a = int(round(target_timewidth_s / timestep_sec))

            self.log.debug('Using solution intervals of {0} (fast) and {1} '
                '(slow) time slots'.format(self.solint_time_p, self.solint_time_a))

            # Set frequency intervals for selfcal solve steps. The interval for
            # slow (amp) selfcal should be the number of channels in a band after
            # averaging. The interval for fast (phase) selfcal should be the
            # number of channels in 20 MHz or less
            num_chan_per_band_after_avg = nchan / self.facetselfcal_freqstep
            if self.dynamic_range.lower() == 'hd':
                # For high-dynamic range solve, the interval for slow (amp) selfcal
                # should be 1 (every channel)
                self.solint_freq_a = 1
            else:
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
                    'input_files_single_mapfile' in d and
                    'subtracted_data_colname' in d):
                    self.diff_models_field_mapfile = d['diff_models_field_mapfile']
                    self.input_files_single_mapfile = d['input_files_single_mapfile']
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
