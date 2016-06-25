"""
Module that holds all parset-related functions
"""
import sys
import os
import glob
import logging
import ConfigParser
import numpy as np
from factor._logging import set_log_file

log = logging.getLogger('factor:parset')


def parset_read(parset_file, use_log_file=True):
    """
    Read a Factor-formatted parset file and return dict of parameters

    Parameters
    ----------
    parset_file : str
        Filename of Factor-formated parset file
    use_log_file : bool, optional
        Use a log file as well as outputing to the screen

    Returns
    -------
    parset_dict : dict
        Dict of parset parameters

    """
    if not os.path.isfile(parset_file):
        log.critical("Missing parset file ({}), I don't know what to do :'(".format(parset_file))
        sys.exit(1)

    log.info("Reading parset file: {}".format(parset_file))
    parset = ConfigParser.RawConfigParser()
    parset.read(parset_file)

    # Handle global parameters
    parset_dict = get_global_options(parset)

    # Handle calibration parameters
    parset_dict['calibration_specific'].update(get_calibration_options(parset))

    # Handle imaging parameters
    parset_dict['imaging_specific'].update(get_imaging_options(parset))

    # Handle directions-related parameters
    parset_dict['direction_specific'].update(get_directions_options(parset))

    # Handle cluster-specific parameters
    parset_dict['cluster_specific'].update(get_cluster_options(parset))

    # Handle checkfactor parameters
    parset_dict['checkfactor'].update(get_checkfactor_options(parset))

    # Set up working directory. All output will be placed in this directory
    if '+' in parset_dict['dir_working']:
        # Check if "+" is in path, as casapy buildmytasks will not work
        # correctly (due to its use of sed)
        log.critical("A '+' appears in the working dir path {}. FACTOR's custom "
            "CASA ft task will not work correctly".format(parset_dict['dir_working']))
        sys.exit(1)
    if not os.path.isdir(parset_dict['dir_working']):
        os.mkdir(parset_dict['dir_working'])
    try:
        os.chdir(parset_dict['dir_working'])
        for subdir in ['logs', 'state', 'results', 'regions', 'chunks']:
            subdir_path = os.path.join(parset_dict['dir_working'], subdir)
            if not os.path.isdir(subdir_path):
                os.mkdir(subdir_path)
    except Exception as e:
        log.critical("Cannot use the working dir {0}: {1}".format(parset_dict['dir_working'], e))
        sys.exit(1)
    if use_log_file:
        set_log_file(os.path.join(parset_dict['dir_working'], 'factor.log'))
    log.info("=========================================================\n")
    log.info("Working directory is {0}".format(parset_dict['dir_working']))

    # Get all the MS files in the input directory. These are identified by the
    # extensions 'ms', 'MS', 'dpppaverage' or 'ndppp_prep_target' (both used by
    # various versions of the pre-facet pipeline)
    ms_files = []
    for exten in ['MS', 'ms', 'dpppaverage', 'ndppp_prep_target']:
        ms_files += glob.glob(os.path.join(parset_dict['dir_ms'], '*.{}'.format(exten)))
    parset_dict['mss'] = sorted(ms_files)
    if len(parset_dict['mss']) == 0:
        log.error('No MS files found in {}!'.format(parset_dict['dir_ms']))
        sys.exit(1)
    log.info("Input MS directory is {0}".format(parset_dict['dir_ms']))
    log.info("Working on {} input files.".format(len(parset_dict['mss'])))

    # Handle MS-specific parameters
    parset_dict.update(get_ms_options(parset, parset_dict['mss'], parset_dict['skymodel_extension']) )

    # Check for unused sections
    given_sections = parset._sections.keys()
    msfiles = [os.path.basename(m) for m in parset_dict['mss']]
    allowed_sections = ['global', 'calibration', 'imaging', 'directions',
        'cluster', 'checkfactor'] + msfiles
    for section in given_sections:
        if section not in allowed_sections:
            log.warning('Section "{}" was given in the parset but is not a valid '
                'section name'.format(section))

    return parset_dict


def get_global_options(parset):
    """
    Handle the global options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all global options

    """
    parset_dict = parset._sections['global'].copy()
    parset_dict.update({'direction_specific': {}, 'calibration_specific': {},
                        'imaging_specific': {}, 'cluster_specific': {},
                        'checkfactor': {}})

    # Parmdb name for dir-indep. selfcal solutions (stored inside the input band
    # measurement sets, so path should be relative to those; default =
    # instrument_directionindependent)
    if 'parmdb_name' not in parset_dict:
        parset_dict['parmdb_name'] = 'instrument_directionindependent'

    # Extension that when concatenated with the 'extension-stripped' MS path gives
    # a path that is checked if it contains a skymodel
    if 'skymodel_extension' not in parset_dict:
        parset_dict['skymodel_extension'] = '.wsclean_low2-model.merge'

    # Size of time chunks in seconds (default = 2400; minimum allowed value is
    # 1200). Generally, the number of chunks should be at least the number of
    # available CPUs.
    if 'chunk_size_sec' in parset_dict:
        parset_dict['chunk_size_sec'] = parset.getfloat('global', 'chunk_size_sec')
    else:
        parset_dict['chunk_size_sec'] = 2400.0

    # Use interactive mode (default = False). Factor will ask for confirmation of
    # internally derived DDE calibrators and facets
    if 'interactive' in parset_dict:
        parset_dict['interactive'] = parset.getboolean('global', 'interactive')
    else:
        parset_dict['interactive'] = False

    # Make final mosaic (default = True)
    if 'make_mosaic' in parset_dict:
        parset._sections['imaging']['make_mosaic'] = parset_dict['make_mosaic']

    # Exit if selfcal fails for any direction (default = True). If False, processing
    # will continue and the failed direction will receive the selfcal solutions of
    # the nearest successful direction unless skip_selfcal_check is True, in which
    # case processing continues as if the selfcal succeeded
    if 'exit_on_selfcal_failure' in parset_dict:
        parset._sections['calibration']['exit_on_selfcal_failure'] = parset_dict['exit_on_selfcal_failure']
    if 'skip_selfcal_check' in parset_dict:
        parset._sections['calibration']['skip_selfcal_check'] = parset_dict['skip_selfcal_check']

    # Padding factor for WSClean images (default = 1.6)
    if 'wsclean_image_padding' in parset_dict:
        parset._sections['imaging']['wsclean_image_padding'] = parset_dict['wsclean_image_padding']

    # Padding factor for WSClean images (default = 1.4)
    if 'wsclean_model_padding' in parset_dict:
        parset._sections['imaging']['wsclean_model_padding'] = parset_dict['wsclean_model_padding']

    # Use WSClean or CASA for imaging of entire facet (default = wsclean). For large
    # bandwidths, the CASA imager is typically faster
    if 'facet_imager' in parset_dict:
        parset._sections['imaging']['facet_imager'] = parset_dict['facet_imager']

    # Keep calibrated data for each facet (default = True for averaged data and
    # False for unaveraged data). If a target is specified (see below), the averaged
    # data for the target is always kept, regardless of this setting. If the
    # averaged data are kept, reimaging will be dramatically faster if multiple
    # images per facet are made
    if 'keep_avg_facet_data' in parset_dict:
        parset_dict['keep_avg_facet_data'] = parset.getboolean('global', 'keep_avg_facet_data')
    else:
        parset_dict['keep_avg_facet_data'] = True
    if 'keep_unavg_facet_data' in parset_dict:
        parset_dict['keep_unavg_facet_data'] = parset.getboolean('global', 'keep_unavg_facet_data')
    else:
        parset_dict['keep_unavg_facet_data'] = False

    # Maximum number of cycles of the last step of selfcal to perform (default =
    # 10). The last step is looped until the number of cycles reaches this value or
    # until the improvement in dynamic range over the previous image is less than
    # 1.25%
    if 'max_selfcal_loops' in parset_dict:
        parset._sections['calibration']['max_selfcal_loops'] = parset_dict['max_selfcal_loops']

    # Use baseline-dependent preaveraging to increase the signal-to-noise of the
    # phase-only solve for sources below this flux (default = 0.0; i.e., disabled).
    # When activated, averaging in time is done to exploit the time coherence in the
    # TEC solutions
    if 'preaverage_flux_jy' in parset_dict:
        parset._sections['calibration']['preaverage_flux_jy'] = parset_dict['preaverage_flux_jy']

    # Peel the calibrator for sources above this flux density (default = 25.0).
    # When activated, the calibrator is peeled using a supplied sky model and
    # the facet is then imaged as normal. Note: a sky model must be specified in the
    # directions file in the peel_skymodel column for each source that should be
    # peeled
    if 'peel_flux_jy' in parset_dict:
        parset._sections['calibration']['peel_flux_jy'] = parset_dict['peel_flux_jy']

    # Use multi-scale selfcal that starts at 20 arcsec resolution and increases the
    # resolution in stages to the full resolution (default = False). This method may
    # improve convergence, especially when the starting model is poor
    if 'multiscale_selfcal' in parset_dict:
        log.warning('Option "multiscale_selfcal" is deprecated and should be changed to "multires_selfcal"')
        parset._sections['calibration']['multires_selfcal'] = parset_dict['multiscale_selfcal']

    # Max desired peak flux density reduction at center of the facet edges due to
    # bandwidth smearing (at the mean frequency) and time smearing (default = 0.15 =
    # 15% reduction in peak flux). Higher values result in shorter run times but
    # more smearing away from the facet centers. This value only applies to the
    # facet imaging (selfcal always uses a value of 0.15)
    if 'max_peak_smearing' in parset_dict:
        parset._sections['imaging']['max_peak_smearing'] = parset_dict['max_peak_smearing']

    # Size of frequency block in MHz over which a single TEC solution is fit
    # (default = 10)
    if 'tec_block_mhz' in parset_dict:
        parset._sections['calibration']['tec_block_mhz'] = parset_dict['tec_block_mhz']

    # Selfcal imaging parameters: pixel size in arcsec (default = 1.5) and Briggs
    # robust parameter (default = -0.25)
    if 'selfcal_cellsize_arcsec' in parset_dict:
        parset._sections['imaging']['selfcal_cellsize_arcsec'] = parset_dict['selfcal_cellsize_arcsec']
    if 'selfcal_robust' in parset_dict:
        parset._sections['imaging']['selfcal_robust'] = parset_dict['selfcal_robust']


    # Check for unused options
    given_options = parset.options('global')
    allowed_options = ['dir_working', 'dir_ms', 'parmdb_name', 'interactive',
        'make_mosaic', 'exit_on_selfcal_failure', 'skip_selfcal_check',
        'wsclean_nbands', 'facet_imager', 'keep_avg_facet_data', 'chunk_size_sec',
        'wsclean_image_padding', 'wsclean_model_padding', 'peel_flux_jy',
        'keep_unavg_facet_data', 'max_selfcal_loops', 'preaverage_flux_jy',
        'multiscale_selfcal', 'skymodel_extension', 'max_peak_smearing',
        'tec_block_mhz', 'selfcal_cellsize_arcsec', 'selfcal_robust']
    allowed_options.extend(['direction_specific', 'calibration_specific',
        'imaging_specific', 'cluster_specific']) # add dicts needed for deprecated options
    deprecated_options_imaging = ['make_mosaic', 'facet_imager',
        'max_peak_smearing', 'selfcal_cellsize_arcsec', 'selfcal_robust']
    deprecated_options_cal = ['exit_on_selfcal_failure',
        'skip_selfcal_check', 'max_selfcal_loops', 'preaverage_flux_jy',
        'multiscale_selfcal', 'tec_block_mhz', 'peel_flux_jy']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [global] section of the '
                'parset but is not a valid global option'.format(option))
        if option in deprecated_options_imaging:
            log.warning('Option "{}" was given in the [global] section of the '
                'parset but should be in the [imaging] section'.format(option))
        if option in deprecated_options_cal:
            log.warning('Option "{}" was given in the [global] section of the '
                'parset but should be in the [calibration] section'.format(option))

    return parset_dict


def get_calibration_options(parset):
    """
    Handle the calibration options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all calibration options

    """
    if 'calibration' in parset._sections.keys():
        parset_dict = parset._sections['calibration']
        given_options = parset.options('calibration')
    else:
        parset_dict = {}
        given_options = []

    # Exit if selfcal fails for any direction (default = True). If False, processing
    # will continue and the failed direction will receive the selfcal solutions of
    # the nearest successful direction unless skip_selfcal_check is True, in which
    # case processing continues as if the selfcal succeeded
    if 'exit_on_selfcal_failure' in parset_dict:
        parset_dict['exit_on_selfcal_failure'] = parset.getboolean('calibration',
            'exit_on_selfcal_failure')
    else:
        parset_dict['exit_on_selfcal_failure'] = True
    if 'skip_selfcal_check' in parset_dict:
        parset_dict['skip_selfcal_check'] = parset.getboolean('calibration',
            'skip_selfcal_check')
    else:
        parset_dict['skip_selfcal_check'] = False

    # Maximum number of cycles of the last step of selfcal to perform (default =
    # 10). The last step is looped until the number of cycles reaches this value or
    # until the improvement in dynamic range over the previous image is less than
    # 1.25%
    if 'max_selfcal_loops' in parset_dict:
        parset_dict['max_selfcal_loops'] = parset.getint('calibration', 'max_selfcal_loops')
    else:
        parset_dict['max_selfcal_loops'] = 10

    # Use baseline-dependent preaveraging to increase the signal-to-noise of the
    # phase-only solve for sources below this flux (default = 0.0; i.e., disabled).
    # When activated, averaging in time is done to exploit the time coherence in the
    # TEC solutions
    if 'preaverage_flux_jy' in parset_dict:
        parset_dict['preaverage_flux_jy'] = parset.getfloat('calibration', 'preaverage_flux_jy')
    else:
        parset_dict['preaverage_flux_jy'] = 0.0

    # Use multi-resolution selfcal that starts at 20 arcsec resolution and increases the
    # resolution in stages to the full resolution (default = False). This method may
    # improve convergence, especially when the starting model is poor
    if 'multires_selfcal' in parset_dict:
        parset_dict['multires_selfcal'] = parset.getboolean('calibration', 'multires_selfcal')
    elif 'multiscale_selfcal' in parset_dict:
        log.warning('Option "multiscale_selfcal" is deprecated and should be changed to "multires_selfcal"')
        parset_dict['multires_selfcal'] = parset.getboolean('calibration', 'multiscale_selfcal')
    else:
        parset_dict['multires_selfcal'] = False

    # Size of frequency block in MHz over which a single TEC solution is fit
    # (default = 10)
    if 'tec_block_mhz' in parset_dict:
        parset_dict['tec_block_mhz'] = parset.getfloat('calibration', 'tec_block_mhz')
    else:
        parset_dict['tec_block_mhz'] = 10.0

    # Peel the calibrator for sources above this flux density (default = 25.0).
    # When activated, the calibrator is peeled using a supplied sky model and
    # the facet is then imaged as normal. Note: a sky model must be specified in the
    # directions file in the peel_skymodel column for each source that should be
    # peeled
    if 'peel_flux_jy' in parset_dict:
        parset_dict['peel_flux_jy'] = parset.getfloat('calibration', 'peel_flux_jy')
    else:
        parset_dict['peel_flux_jy'] = 25.0

    # Minimum uv distance in lambda for calibration (default = 80)
    if 'solve_min_uv_lambda' in parset_dict:
        parset_dict['solve_min_uv_lambda'] = parset.getfloat('calibration', 'solve_min_uv_lambda')
    else:
        parset_dict['solve_min_uv_lambda'] = 80.0

    # Smooth amplitudes with spline fit + 2-D median (default = True, i.e., smooth
    # with a 1-D median only)
    if 'spline_smooth2d' in parset_dict:
        parset_dict['spline_smooth2d'] = parset.getboolean('calibration', 'spline_smooth2d')
    else:
        parset_dict['spline_smooth2d'] = True

    # Include XY and YX correlations during the slow gain solve for sources above
    # this flux density (default = 1000.0). Below this value, only the XX and YY
    # correlations are included. Note that spline_smooth2D must be True to use this
    # option
    if 'solve_all_correlations_flux_jy' in parset_dict:
        parset_dict['solve_all_correlations_flux_jy'] = parset.getfloat('calibration', 'solve_all_correlations_flux_jy')
    else:
        parset_dict['solve_all_correlations_flux_jy'] = 1000.0

    # Check for unused options
    allowed_options = ['exit_on_selfcal_failure', 'skip_selfcal_check',
        'max_selfcal_loops', 'preaverage_flux_jy', 'multiscale_selfcal', 'multires_selfcal',
        'tec_block_mhz', 'peel_flux_jy', 'solve_min_uv_lambda', 'spline_smooth2d',
        'solve_all_correlations_flux_jy']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [calibration] section of the '
                'parset but is not a valid calibration option'.format(option))

    return parset_dict


def get_imaging_options(parset):
    """
    Handle the imaging options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all imaging options

    """
    if 'imaging' in parset._sections.keys():
        parset_dict = parset._sections['imaging']
        given_options = parset.options('imaging')
    else:
        parset_dict = {}
        given_options = []


    # Make final mosaic (default = True)
    if 'make_mosaic' in parset_dict:
        parset_dict['make_mosaic'] = parset.getboolean('imaging', 'make_mosaic')
    else:
        parset_dict['make_mosaic'] = True

    # Re-image directions for which selfcal was successful (default = True)
    if 'reimage_selfcaled' in parset_dict:
        parset_dict['reimage_selfcaled'] = parset.getboolean('imaging',
            'reimage_selfcaled')
    elif 'reimage' in parset._sections['directions']:
        log.warning('Option "reimage" was given in the [directions] section of the '
            'parset but should be in the [imaging] section and should be changed to '
            '"reimage_selfcaled"')
        parset_dict['reimage_selfcaled'] = parset.getboolean('directions',
            'reimage')
    else:
        parset_dict['reimage_selfcaled'] = True

    # Max factor used to set the number of WSClean channel images when wide-band
    # clean is used (default = 4). The number of channel images is determined by
    # dividing the number of bands by the nearest divisor to this factor. Smaller
    # values produce better results but require longer run times. Wide-band clean is
    # activated when there are more than 5 bands
    if 'wsclean_nchannels_factor' in parset_dict:
        parset_dict['wsclean_nchannels_factor'] = parset.getint('imaging', 'wsclean_nchannels_factor')
        if parset_dict['wsclean_nchannels_factor'] < 1:
            parset_dict['wsclean_nchannels_factor'] = 1
    else:
        parset_dict['wsclean_nchannels_factor'] = 4

    # Allow flagged data to be added during WSClean imaging to allow
    # wsclean_nchannels_factor to be a divisor of the number bands (default = True).
    # Enabling this option can dramatically speed up imaging with WSClean when the
    # number of bands before padding does not allow wsclean_nchannels_factor to be
    # greater than 1 (e.g., wsclean_nchannels_factor must be 1 to be an even divisor
    # of 29 bands, so activating this option would add 1 band of flagged data to
    # produce 30 bands, which will work with wsclean_nchannels_factor = 3, 5, or 6)
    if 'wsclean_add_bands' in parset_dict:
        parset_dict['wsclean_add_bands'] = parset.getboolean('imaging', 'wsclean_add_bands')
    else:
        parset_dict['wsclean_add_bands'] = True

    # Use WSClean or CASA for imaging of entire facet (default = wsclean). For large
    # bandwidths, the CASA imager is typically faster
    if 'facet_imager' not in parset_dict:
        parset_dict['facet_imager'] = 'wsclean'

    # Max desired peak flux density reduction at center of the facet edges due to
    # bandwidth smearing (at the mean frequency) and time smearing (default = 0.15 =
    # 15% reduction in peak flux). Higher values result in shorter run times but
    # more smearing away from the facet centers. This value only applies to the
    # facet imaging (selfcal always uses a value of 0.15)
    if 'max_peak_smearing' in parset_dict:
        parset_dict['max_peak_smearing'] = parset.getfloat('imaging', 'max_peak_smearing')
    else:
        parset_dict['max_peak_smearing'] = 0.15

    # Selfcal imaging parameters: pixel size in arcsec (default = 1.5), Briggs
    # robust parameter (default = -0.25 for casa and -0.5 for wsclean), minimum uv
    # distance in lambda (default = 80), and multiscale clean scales (default = [0,
    # 3, 7, 25, 60, 150]). These settings apply both to selfcal images and to the
    # full facet image used to make the improved facet model that is subtracted from
    # the data
    if 'selfcal_cellsize_arcsec' in parset_dict:
        parset_dict['selfcal_cellsize_arcsec'] = parset.getfloat('imaging', 'selfcal_cellsize_arcsec')
    else:
        parset_dict['selfcal_cellsize_arcsec'] = 1.5
    if 'selfcal_robust' in parset_dict:
        parset_dict['selfcal_robust'] = parset.getfloat('imaging', 'selfcal_robust')
    else:
        parset_dict['selfcal_robust'] = -0.25
    if 'selfcal_robust_wsclean' in parset_dict:
        parset_dict['selfcal_robust_wsclean'] = parset.getfloat('imaging', 'selfcal_robust_wsclean')
    else:
        parset_dict['selfcal_robust_wsclean'] = -0.5
    if 'selfcal_min_uv_lambda' in parset_dict:
        parset_dict['selfcal_min_uv_lambda'] = parset.getfloat('imaging', 'selfcal_min_uv_lambda')
    else:
        parset_dict['selfcal_min_uv_lambda'] = 80.0
    if not 'selfcal_scales' in parset_dict:
        parset_dict['selfcal_scales'] = '[0, 3, 7, 25, 60, 150]'

    # Use a clean threshold during selfcal imaging (default = False). If False,
    # clean will always stop at 1000 iterations. If True, clean will go to 1 sigma
    # noise level
    if 'selfcal_clean_threshold' in parset_dict:
        parset_dict['selfcal_clean_threshold'] = parset.getboolean('imaging', 'selfcal_clean_threshold')
    else:
        parset_dict['selfcal_clean_threshold'] = False

    # Use an adaptive masking threshold during selfcal imaging (default = False). If
    # True, the masking threshold will be estimated using the negative peaks in the
    # image, which can help selfcal convergence in the presence of strong artifacts
    if 'selfcal_adaptive_threshold' in parset_dict:
        parset_dict['selfcal_adaptive_threshold'] = parset.getboolean('imaging', 'selfcal_adaptive_threshold')
    else:
        parset_dict['selfcal_adaptive_threshold'] = False

    # Facet imaging parameters: pixel size in arcsec, Briggs robust parameter, uv
    # taper in arcsec, and minimum uv distance in lambda. These parameters are used
    # only for making full facet images (and not for making improved models). One
    # set of images and one mosaic image will be made for each set of parameters. By
    # default, facets will be imaged using the selfcal imaging parameters above
    len_list = []
    if 'facet_cellsize_arcsec' in parset_dict:
        val_list = parset_dict['facet_cellsize_arcsec'].strip('[]').split(',')
        val_list = [float(v) for v in val_list]
        parset_dict['facet_cellsize_arcsec'] = val_list
        len_list.append(len(val_list))
    if 'facet_taper_arcsec' in parset_dict:
        val_list = parset_dict['facet_taper_arcsec'].strip('[]').split(',')
        val_list = [float(v) for v in val_list]
        parset_dict['facet_taper_arcsec'] = val_list
        len_list.append(len(val_list))
    if 'facet_robust' in parset_dict:
        val_list = parset_dict['facet_robust'].strip('[]').split(',')
        val_list = [float(v) for v in val_list]
        parset_dict['facet_robust'] = val_list
        len_list.append(len(val_list))
    if 'facet_min_uv_lambda' in parset_dict:
        val_list = parset_dict['facet_min_uv_lambda'].strip('[]').split(',')
        val_list = [float(v) for v in val_list]
        parset_dict['facet_min_uv_lambda'] = val_list
        len_list.append(len(val_list))

    # Check that all the above options have the same number of entries
    if len(set(len_list)) == 0:
        nvals = 1
    elif len(set(len_list)) == 1:
        nvals = len_list[0]
    else:
        log.error('The options facet_cellsize_arcsec, facet_taper_arcsec, facet_robust, and '
            'facet_min_uv_lambda must all have the same number of entires')
        sys.exit(1)

    # Set defaults for any that did not have entries
    if 'facet_cellsize_arcsec' not in parset_dict:
        parset_dict['facet_cellsize_arcsec'] = [parset_dict['selfcal_cellsize_arcsec']] * nvals
    if 'facet_taper_arcsec' not in parset_dict:
        parset_dict['facet_taper_arcsec'] = [0.0] * nvals
    if 'facet_robust' not in parset_dict:
        if parset_dict['facet_imager'] == 'wsclean':
            selfcal_robust = parset_dict['selfcal_robust_wsclean']
        else:
            selfcal_robust = parset_dict['selfcal_robust']

        parset_dict['facet_robust'] = [parset_dict['selfcal_robust']] * nvals
    if 'facet_min_uv_lambda' not in parset_dict:
        parset_dict['facet_min_uv_lambda'] = [parset_dict['selfcal_min_uv_lambda']] * nvals

   # Padding factor for WSClean images (default = 1.6)
    if 'wsclean_image_padding' in parset_dict:
        parset_dict['wsclean_image_padding'] = parset.getfloat('imaging', 'wsclean_image_padding')
    else:
        parset_dict['wsclean_image_padding'] = 1.6

    # Padding factor for WSClean images (default = 1.4)
    if 'wsclean_model_padding' in parset_dict:
        parset_dict['wsclean_model_padding'] = parset.getfloat('imaging', 'wsclean_model_padding')
    else:
        parset_dict['wsclean_model_padding'] = 1.4

    # Check for unused options
    allowed_options = ['make_mosaic', 'wsclean_nchannels_factor', 'facet_imager',
        'max_peak_smearing', 'selfcal_cellsize_arcsec', 'selfcal_robust',
        'selfcal_robust_wsclean', 'selfcal_clean_threshold', 'selfcal_adaptive_threshold',
        'facet_cellsize_arcsec', 'facet_taper_arcsec', 'facet_robust',
        'reimage_selfcaled', 'wsclean_image_padding', 'wsclean_model_padding',
        'selfcal_min_uv_lambda', 'facet_min_uv_lambda', 'wsclean_add_bands',
        'selfcal_robust_wsclean']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [imaging] section of the '
                'parset but is not a valid imaging option'.format(option))

    return parset_dict


def get_directions_options(parset):
    """
    Handle the directions options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all directions options

    """
    if 'directions' in parset._sections.keys():
        parset_dict = parset._sections['directions']
        given_options = parset.options('directions')
    else:
        parset_dict = {}
        given_options = []

    # Full path to sky model (in makesourcedb format) to be used for calibrator
    # selection and facet-boundary source avoidance (default is to use
    # direction-independent sky model of the highest-frequency band). The sky
    # model must be grouped into patches by source (in PyBDSM, this grouping can be
    # done by setting bbs_patches = 'source' in the write_catalog task)
    if 'faceting_skymodel' not in parset_dict:
        parset_dict['faceting_skymodel'] = None

    # Check whether any sources from the initial subtract sky model fall on facet
    # edges. If any are found, the facet regions are adjusted to avoid them (default
    # is True)
    if 'check_edges' in parset_dict:
        parset_dict['check_edges'] = parset.getboolean('directions', 'check_edges')
    else:
        parset_dict['check_edges'] = True

    # Parameters for selecting directions internally (radius from phase center
    # within which to consider sources as potential calibrators, min flux, max size
    # of a source, and max separation between sources below which they are grouped
    # into one direction; required if no directions_file is given). The number of
    # internally derived directions can be limited to a maximum number of directions
    # if desired with max_num (default = all). Lastly, faceting_radius_deg sets the
    # radius within which facets will be used; outside of this radius, small patches
    # are used that do not appear in the final mosaic.  These parameters will
    # determine the faceting of the field

    # Radius from phase center within which to consider sources as potential
    # calibrators (default = 2 * FWHM of primary beam of highest-frequency band)
    if 'max_radius_deg' in parset_dict:
        parset_dict['max_radius_deg'] = parset.getfloat('directions',
            'max_radius_deg')
    else:
        parset_dict['max_radius_deg'] = None

    # If no directions_file is given, the selection criteria for calibrator sources
    # that follow must be given. For merging of multiple sources into one calibrator
    # group, flux_min_for_merging_Jy (default = 0.1 Jy) and size_max_arcmin set the min
    # flux density and max size of individual sources to be considered for grouping,
    # and separation_max_arcmin sets the max separation between sources below which
    # they are grouped into one calibrator. After grouping, flux_min_Jy sets the
    # min total flux density of a source (or group) to be considered as a DDE
    # calibrator
    if 'flux_min_for_merging_jy' in parset_dict:
        parset_dict['flux_min_for_merging_jy'] = parset.getfloat('directions',
            'flux_min_for_merging_jy')
    else:
        parset_dict['flux_min_for_merging_jy'] = 0.1
    if 'size_max_arcmin' in parset_dict:
        parset_dict['size_max_arcmin'] = parset.getfloat('directions',
            'size_max_arcmin')
    else:
        parset_dict['size_max_arcmin'] = None
    if 'separation_max_arcmin' in parset_dict:
        parset_dict['separation_max_arcmin'] = parset.getfloat('directions',
            'separation_max_arcmin')
    else:
        parset_dict['separation_max_arcmin'] = None
    if 'flux_min_jy' in parset_dict:
        parset_dict['flux_min_jy'] = parset.getfloat('directions',
            'flux_min_jy')
    else:
        parset_dict['flux_min_jy'] = None

    # Number of internally derived directions can be limited to a maximum number
    # of directions if desired with max_num (default = all).
    if 'ndir_max' in parset_dict:
        parset_dict['ndir_max'] = parset.getint('directions', 'ndir_max')
    elif 'max_num' in parset_dict:
        log.warning('Option "max_num" is deprecated and should be changed to "ndir_max"')
        parset_dict['ndir_max'] = parset.getint('directions', 'max_num')
    else:
        parset_dict['ndir_max'] = None

    # Radius within which facets will be used (default = 1.25 * FWHM of primary beam
    # of highest-frequency band); outside of this radius, small patches are used
    # that do not appear in the final mosaic.
    if 'faceting_radius_deg' in parset_dict:
        parset_dict['faceting_radius_deg'] = parset.getfloat('directions',
            'faceting_radius_deg')
    else:
        parset_dict['faceting_radius_deg'] = None

    # Grouping of directions into groups that are selfcal-ed in parallel, defined as
    # grouping:n_total_per_grouping. For example, groupings = 1:5, 4:0 means two
    # groupings are used, with the first 5 directions put into groups of one (i.e.,
    # each direction processed in series) and the rest of the directions divided
    # into groups of 4 (i.e., 4 directions processed in parallel). Default is one at
    # a time (i.e., groupings = 1:0)
    if 'groupings' in parset_dict:
        groupings=[]
        keys = []
        vals = []
        kvs = parset_dict['groupings'].split(',')
        for kv in kvs:
            key, val = kv.split(':')
            keys.append(key.strip())
            vals.append(val.strip())
        for key, val in zip(keys, vals):
            groupings.append({key: int(val)})
        parset_dict['groupings'] = groupings
    else:
        parset_dict['groupings'] = [{'1':0}]
    log.info("Using the following groupings for directions: {}"
        .format(', '.join(['{0}:{1}'.format(n.keys()[0], n.values()[0])
        for n in parset_dict['groupings']])))

    # If groups are used to process more than one direction in parallel, reordering
    # of the directions in the groups can be done to maximize the flux-weighted
    # separation between directions in each group (default = True)
    if 'allow_reordering' in parset_dict:
        parset_dict['allow_reordering'] = parset.getboolean('directions',
            'allow_reordering')
    else:
        parset_dict['allow_reordering'] = True

    # Total number of directions to selfcal (default = all)
    if 'ndir_selfcal' in parset_dict:
        parset_dict['ndir_selfcal'] = parset.getint('directions', 'ndir_selfcal')
        if parset_dict['ndir_selfcal'] < 1:
            log.error('Total number of directions to selfcal must be 1 or more')
            sys.exit(1)
        log.info("Self calibrating up to %s direction(s)" % (parset_dict['ndir_selfcal']))
    else:
        parset_dict['ndir_selfcal'] = None

    # Total number of directions to process (default = all). If this number is
    # greater than ndir_selfcal, then the remaining directions will not be selfcal-
    # ed but will instead be imaged with the selfcal solutions from the nearest
    # direction for which selfcal succeeded (if a target is specified and
    # target_has_own_facet = True, it will be imaged in this way after ndir_total
    # number of directions are processed)
    if 'ndir_process' in parset_dict:
        parset_dict['ndir_process'] = parset.getint('directions', 'ndir_process')
        if parset_dict['ndir_process'] < 1:
            log.error('Total number of directions to process must be 1 or more')
            sys.exit(1)
        log.info("Processing up to %s direction(s) in total" % (parset_dict['ndir_process']))
    elif 'ndir_total' in parset_dict:
        log.warning('Option "ndir_total" is deprecated and should be changed to "ndir_process"')
        parset_dict['ndir_process'] = parset.getint('directions', 'ndir_total')
        if parset_dict['ndir_process'] < 1:
            log.error('Total number of directions to process must be 1 or more')
            sys.exit(1)
        log.info("Processing up to %s direction(s) in total" % (parset_dict['ndir_process']))
    else:
        parset_dict['ndir_process'] = None

    # A target can be specified to ensure that it falls entirely within a single
    # facet. The values should be those of a circular region that encloses the
    # source and not those of the target itself. Lastly, the target can be placed in
    # a facet of its own. In this case, it will not go through selfcal but will
    # instead use the selfcal solutions of the nearest facet for which selfcal was
    # done
    if 'target_ra' not in parset_dict:
        parset_dict['target_ra'] = None
    if 'target_dec' not in parset_dict:
        parset_dict['target_dec'] = None
    if 'target_radius_arcmin' in parset_dict:
        parset_dict['target_radius_arcmin'] = parset.getfloat('directions',
            'target_radius_arcmin')
    else:
        parset_dict['target_radius_arcmin'] = None
    if 'target_has_own_facet' in parset_dict:
        parset_dict['target_has_own_facet'] = parset.getboolean('directions',
            'target_has_own_facet')
    else:
        parset_dict['target_has_own_facet'] = False

    # Radius in degrees within which the direction-dependent solutions will be
    # transferred before starting selfcal (default = 0; i.e., disabled). If a
    # direction is within this distance of a calibrator for which selfcal was
    # successful, the dir-dep selfcal solutions from this calibrator will be used
    # instead of the dir-indep ones
    if 'transfer_radius_deg' in parset_dict:
        parset_dict['transfer_radius_deg'] = parset.getfloat('directions',
            'transfer_radius_deg')
    elif 'transfer_radius' in parset_dict:
        log.warning('Option "transfer_radius" is deprecated and should be changed to "transfer_radius_deg"')
        parset_dict['transfer_radius_deg'] = parset.getfloat('directions',
            'transfer_radius')
    else:
        parset_dict['transfer_radius_deg'] = 0.0

    # Check for unused options
    allowed_options = ['faceting_skymodel', 'directions_file', 'max_radius_deg',
        'flux_min_for_merging_jy', 'flux_min_jy', 'size_max_arcmin',
        'separation_max_arcmin', 'max_num', 'ndir_max',
        'faceting_radius_deg', 'check_edges', 'ndir_total', 'ndir_process',
        'ndir_selfcal', 'transfer_radius', 'transfer_radius_deg',
        'groupings', 'allow_reordering', 'target_ra', 'target_dec',
        'target_radius_arcmin', 'target_has_own_facet']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [directions] section of the '
                'parset but is not a valid directions option'.format(option))

    return parset_dict


def get_cluster_options(parset):
    """
    Handle the compute cluster options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all cluster options

    """
    if 'cluster' in parset._sections.keys():
        parset_dict = parset._sections['cluster']
        given_options = parset.options('cluster')
    else:
        parset_dict = {}
        given_options = []

    # Paths to the LOFAR software
    if 'lofarroot' not in parset_dict:
        if 'LOFARROOT' in os.environ:
            parset_dict['lofarroot'] = os.environ['LOFARROOT']
        else:
            log.critical("The LOFAR root directory cannot be determined. Please "
                "specify it in the [cluster] section of the parset as lofarroot")
            sys.exit(1)
    if 'lofarpythonpath' not in parset_dict:
        if parset_dict['lofarroot'] in os.environ['PYTHONPATH']:
            pypaths = os.environ['PYTHONPATH'].split(':')
            for pypath in pypaths:
                if parset_dict['lofarroot'] in pypath:
                    parset_dict['lofarpythonpath'] = pypath
                    break
        else:
            log.critical("The LOFAR Python root directory cannot be determined. "
                "Please specify it in the [cluster] section of the parset as "
                "lofarpythonpath")
            sys.exit(1)

    # Number of CPUs per node to be used.
    if 'ncpu' in parset_dict:
        parset_dict['ncpu'] = parset.getint('cluster', 'ncpu')
    else:
        import multiprocessing
        parset_dict['ncpu'] = multiprocessing.cpu_count()
    log.info("Using up to %s CPU(s) per node" % (parset_dict['ncpu']))

    # Maximum fraction of the total memory per node that WSClean may use (default =
    # 0.9)
    if 'wsclean_fmem' in parset_dict:
        parset_dict['wsclean_fmem'] = parset.getfloat('cluster', 'wsclean_fmem')
        if parset_dict['wsclean_fmem'] > 1.0:
            parset_dict['wsclean_fmem'] = 1.0
    elif 'fmem' in parset_dict:
        log.warning('Option "fmem" is deprecated and should be changed to "wsclean_fmem"')
        parset_dict['wsclean_fmem'] = parset.getfloat('cluster', 'fmem')
        if parset_dict['wsclean_fmem'] > 1.0:
            parset_dict['wsclean_fmem'] = 1.0
    else:
        parset_dict['wsclean_fmem'] = 0.9
    log.info("Using up to {0}% of the memory per node for WSClean jobs".format(parset_dict['wsclean_fmem']*100.0))

    # Maximum number of directions to process in parallel on each node (default =
    # 1). Note that the number of CPUs (set with the ncpu parameter) and the amount
    # of memory available to WSClean (set with the wsclean_fmem parameter) will be
    # divided among the directions on each node
    if 'ndir_per_node' in parset_dict:
        parset_dict['ndir_per_node'] = parset.getint('cluster',
            'ndir_per_node')
    else:
        parset_dict['ndir_per_node'] = 1
    log.info("Processing up to %i direction(s) in parallel per node" %
        (parset_dict['ndir_per_node']))

    # Maximum number of io-intensive threads to run per node (per
    # direction). If unset, defaults to sqrt of the number of compute
    # threads that will be used.

    if 'nthread_io' in parset_dict:
        parset_dict['nthread_io'] = parset.getint('cluster',
            'nthread_io')
    else:
        # code originally in lib/scheduler.py
        max_cpus_per_node =  max(1, int(round(parset_dict['ncpu'] / float(parset_dict['ndir_per_node']))))
        parset_dict['nthread_io'] = int(np.ceil(np.sqrt(max_cpus_per_node)))

    log.info("Will use up to %i IO-intensive thread(s) in parallel per node" %
        (parset_dict['nthread_io']))
    
    # Full path to cluster description file. Use clusterdesc_file = PBS to use the
    # PBS / torque reserved nodes. If not given, the clusterdesc file for a single
    # (i.e., local) node is used
    if 'clusterdesc_file' not in parset_dict:
        parset_dict['clusterdesc_file'] = parset_dict['lofarroot'] + '/share/local.clusterdesc'
        parset_dict['node_list'] = ['localhost']

    # Full path to a local disk on the nodes for I/O-intensive processing. The path
    # must be the same for all nodes. If not given, the default directory in the
    # working directory is used
    if 'dir_local' not in parset_dict:
        parset_dict['dir_local'] = None
    #elif parset_dict['clusterdesc_file'] != 'PBS':
        # The local directory only works when the dppp_scratch.py node recipe is
        # used, which is only done when clusterdesc_file = PBS, so exit if not
    #    log.critical('A local scratch directory can only be used when '
    #        'clusterdesc_file = PBS (i.e., on a cluster with many nodes)')
    #    sys.exit(1)

    # Check for unused options
    allowed_options = ['ncpu', 'fmem', 'wsclean_fmem', 'ndir_per_node',
        'clusterdesc_file', 'cluster_type', 'dir_local',
        'node_list', 'lofarroot', 'lofarpythonpath', 'nthread_io']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [cluster] section of the '
                'parset but is not a valid cluster option'.format(option))

    return parset_dict


def get_ms_options(parset, ms_files, skymodel_extension = '.wsclean_low2-model.merge'):
    """
    Handle the ms-specific options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset
    ms_files : list of str
        MS filenames
    skymodel_extension : str, optional
        search for skymodels that may have this extension

    Returns
    -------
    parset_dict : dict
        Dictionary with all ms-specific options

    """
    parset_dict = {'ms_specific': {}}
    for ms in ms_files:
        msbase = os.path.basename(ms)
        if msbase in parset._sections.keys():
            parset_dict['ms_specific'][msbase] = parset._sections[msbase]
        default_skymodel = os.path.splitext(ms)[0]+skymodel_extension
        if os.path.exists(default_skymodel):
            if msbase in parset_dict['ms_specific'].keys():
                if not 'init_skymodel' in parset_dict['ms_specific'][msbase].keys():
                    parset_dict['ms_specific'][msbase]['init_skymodel'] = default_skymodel
            elif msbase not in parset_dict['ms_specific'].keys():
                parset_dict['ms_specific'][msbase] = {'init_skymodel': default_skymodel}
    return parset_dict


def get_checkfactor_options(parset):
    """
    Handle the options for checkfactor

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset

    Returns
    -------
    parset_dict : dict
        Dictionary with all cluster options

    """
    if 'checkfactor' in parset._sections.keys():
        parset_dict = parset._sections['checkfactor']
        given_options = parset.options('checkfactor')
    else:
        parset_dict = {}
        given_options = []

    if 'facet_viewer' not in parset_dict:
        parset_dict['facet_viewer'] = 'casa'

    if 'ds9_limits' not in parset_dict:
        parset_dict['ds9_limits'] = None

    if 'image_display' not in parset_dict:
        parset_dict['image_display'] = 'display -geometry 800x600'
    elif parset_dict['image_display'] == 'display':
        parset_dict['image_display'] = 'display -geometry 800x600'

    if 'ds9_load_regions' in parset_dict:
        parset_dict['ds9_load_regions'] = parset.getboolean('checkfactor', 'ds9_load_regions')
    else:
        parset_dict['ds9_load_regions'] = False

    # Check for unused options
    allowed_options = ['facet_viewer', 'ds9_limits', 'image_display',
        'ds9_load_regions']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [checkfactor] section of the '
                'parset but is not a valid checkfactor option'.format(option))

    return parset_dict
