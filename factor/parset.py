"""
Module that holds all parset-related functions
"""
import sys
import os
import glob
import logging
import ConfigParser
from factor._logging import set_log_file

log = logging.getLogger('factor:parset')


def parset_read(parset_file):
    """
    Read a Factor-formatted parset file and return dict of parameters

    Parameters
    ----------
    parset_file : str
        Filename of Factor-formated parset file

    Returns
    -------
    parset_dict : dict
        Dict of parset parameters

    """
    if not os.path.isfile(parset_file):
        log.critical("Missing parset file (%s), I don't know what to do :'(" % (parset_file))
        sys.exit(1)

    log.info("Reading parset file: %s" % (parset_file))
    parset = ConfigParser.RawConfigParser()
    parset.read(parset_file)
    parset_dict = {}

    # Handle global parameters
    parset_dict.update(get_global_options(parset))

    # Handle directions-related parameters
    parset_dict.update(get_directions_options(parset))

    # Handle cluster-specific parameters
    parset_dict.update(get_cluster_options(parset))

    # Set up working directory. All output will be placed in this directory
    if not os.path.isdir(parset_dict['dir_working']):
        os.mkdir(parset_dict['dir_working'])
    try:
        os.chdir(parset_dict['dir_working'])
        for subdir in ['/logs', '/state', '/results', '/datamaps', '/regions']:
            if not os.path.isdir(parset_dict['dir_working']+subdir):
                os.mkdir(parset_dict['dir_working']+subdir)
    except:
        log.critical("Cannot use the working dir %s" % (parset_dict['dir_working']))
        sys.exit(1)
    set_log_file(parset_dict['dir_working']+'/factor.log')
    log.info("=========================================================\n")
    log.info("Working directory is {0}".format(parset_dict['dir_working']))

    # Get all the MS files in the input directory. These are identified by the
    # extensions 'ms', 'MS', 'dpppconcat' or 'dpppcopycol' (both used by the
    # pre-facet pipeline)
    ms_files = []
    for exten in ['MS', 'ms', 'dpppconcat', 'dpppcopycol']:
        ms_files += glob.glob(os.path.join(parset_dict['dir_ms'], '*.{}'.format(exten)))
    parset_dict['mss'] = sorted(ms_files)
    if len(parset_dict['mss']) == 0:
        log.error('No MS files found in {0}!'.format(parset_dict['dir_ms']))
        sys.exit(1)
    log.info("Working on %i band(s)" % (len(parset_dict['mss'])))

    # Handle MS-specific parameters
    parset_dict.update(get_ms_options(parset, parset_dict['mss']))

    # Check for unused sections
    given_sections = parset._sections.keys()
    allowed_sections = ['global', 'directions', 'cluster'] + parset_dict['mss']
    for section in given_sections:
        if section not in allowed_sections:
            log.warning('Section "{}" was given in the parset but is not a valid '
                'section name'.format(option))

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
    parset_dict = parset._sections['global']

    # Paths to the LOFAR software
    if 'lofarroot' not in parset_dict:
        if 'LOFARROOT' in os.environ:
            parset_dict['lofarroot'] = os.environ['LOFARROOT']
        else:
            log.critical("The LOFAR root directory cannot be determined. Please "
                "specify it in the parset as lofarroot")
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
                "Please specify it in the parset as lofarpythonpath")
            sys.exit(1)

    # Parmdb name for dir-indep. selfcal solutions (stored inside the input band
    # measurement sets, so path should be relative to those; default = instrument)
    if 'parmdb_name' not in parset_dict:
        parset_dict['parmdb_name'] = 'instrument'

    # Use interactive mode (default = False). Factor will ask for confirmation of
    # internally derived DDE calibrators and facets
    if 'interactive' in parset_dict:
        parset_dict['interactive'] = parset.getboolean('global', 'interactive')
    else:
        parset_dict['interactive'] = False

    # Make final mosaic (default = True)
    if 'make_mosaic' in parset_dict:
        parset_dict['make_mosaic'] = parset.getboolean('global', 'make_mosaic')
    else:
        parset_dict['make_mosaic'] = True

    # Exit if selfcal fails for any direction (default = True). If False, processing
    # will continue and the failed direction will receive the selfcal solutions of
    # the nearest successful direction unless skip_selfcal_check is True, in which
    # case processing continues as if the selfcal succeeded
    if 'exit_on_selfcal_failure' in parset_dict:
        parset_dict['exit_on_selfcal_failure'] = parset.getboolean('global',
            'exit_on_selfcal_failure')
    else:
        parset_dict['exit_on_selfcal_failure'] = True
    if 'skip_selfcal_check' in parset_dict:
        parset_dict['skip_selfcal_check'] = parset.getboolean('global',
            'skip_selfcal_check')
    else:
        parset_dict['skip_selfcal_check'] = False

    # Max number of bands per WSClean image when wide-band clean is used (default =
    # 5). Smaller values produce better results but require longer run times.
    # Wide-band clean is activated when there are more than 5 bands
    if 'wsclean_nbands' in parset_dict:
        parset_dict['wsclean_nbands'] = parset.getint('global', 'wsclean_nbands')
    else:
        parset_dict['wsclean_nbands'] = 3

    # Use WSClean or CASA for imaging of entire facet (default = wsclean). For large
    # bandwidths, the CASA imager is typically faster
    if 'facet_imager' not in parset_dict:
        parset_dict['facet_imager'] = 'wsclean'

    # Keep calibrated data for each facet to allow re-imaging by hand (default =
    # True for averaged data and False for unaveraged data). If a target is
    # specified (see below), the averaged data for the target is always kept,
    # regardless of this setting
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
        parset_dict['max_selfcal_loops'] = parset.getint('global', 'max_selfcal_loops')
    else:
        parset_dict['max_selfcal_loops'] = 10

    # Use baseline-dependent preaveraging to increase the signal-to-noise of the
    # phase-only solve (default = True). If True, averaging in time is done to
    # exploit the time coherence in the TEC solutions
    if 'preaverage' in parset_dict:
        parset_dict['preaverage'] = parset.getboolean('global', 'preaverage')
    else:
        parset_dict['preaverage'] = True

    # Check for unused options
    given_options = parset.options('global')
    allowed_options = ['dir_working', 'dir_ms', 'lofarroot',
        'lofarpythonpath', 'parmdb_name', 'interactive', 'make_mosaic',
        'exit_on_selfcal_failure', 'skip_selfcal_check', 'wsclean_nbands',
        'facet_imager', 'keep_avg_facet_data', 'keep_unavg_facet_data',
        'max_selfcal_loops', 'preaverage']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [global] section of the '
                'parset but is not a valid global option'.format(option))

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
        parset_dict = {'direction_specific': parset._sections['directions']}
    else:
        parset_dict = {'direction_specific': {}}

    # Check whether any sources from the initial subtract sky model fall on facet
    # edges. If any are found, the facet regions are adjusted to avoid them (default
    # is False)
    if 'check_edges' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['check_edges'] = parset.getboolean('directions',
            'check_edges')
    else:
        parset_dict['direction_specific']['check_edges'] = False

    # Flux and size cuts for selecting directions internally (min flux, max size of
    # a source, and max separation between sources below which they are grouped into
    # one direction; required if no directions_file is given). The number of
    # internally derived directions can be limited to a maximum number of directions
    # if desired with max_num (default = all). These parameters will determine the
    # faceting of the field
    if 'flux_min_jy' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['flux_min_jy'] = parset.getfloat('directions',
            'flux_min_jy')
    if 'max_num' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['max_num'] = parset.getint('directions',
            'max_num')
    else:
        parset_dict['direction_specific']['max_num'] = None
    if 'size_max_arcmin' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['size_max_arcmin'] = parset.getfloat('directions',
            'size_max_arcmin')
    if 'separation_max_arcmin' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['separation_max_arcmin'] = parset.getfloat('directions',
            'separation_max_arcmin')

    # Grouping of directions into groups that are selfcal-ed in parallel, defined as
    # grouping:n_total_per_grouping. For example, groupings = 1:5, 4:0 means two
    # groupings are used, with the first 5 directions put into groups of one (i.e.,
    # each direction processed in series) and the rest of the directions divided
    # into groups of 4 (i.e., 4 directions processed in parallel). Default is one at
    # a time (i.e., groupings = 1:0)
    if 'groupings' in parset_dict['direction_specific']:
        groupings={}
        keys = []
        vals = []
        kvs = parset_dict['direction_specific']['groupings'].split(',')
        for kv in kvs:
            key, val = kv.split(':')
            keys.append(key.strip())
            vals.append(val.strip())
        for key, val in zip(keys, vals):
            groupings[key] = int(val)
        parset_dict['direction_specific']['groupings'] = groupings
    else:
        parset_dict['direction_specific']['groupings'] = {'1':0}
    log.info("Using the following groupings for directions: {0}"
        .format(parset_dict['direction_specific']['groupings']))

    # If groups are used to process more than one direction in parallel, reordering
    # of the directions in the groups can be done to maximize the flux-weighted
    # separation between directions in each group (default = True)
    if 'allow_reordering' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['allow_reordering'] = parset.getboolean('directions',
            'allow_reordering')
    else:
        parset_dict['direction_specific']['allow_reordering'] = True

    # Total number of directions to selfcal (default = all)
    if 'ndir_selfcal' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['ndir_selfcal'] = parset.getint('directions', 'ndir_selfcal')
        log.info("Self calibrating up to %s direction(s)" % (parset_dict['direction_specific']['ndir_selfcal']))
    else:
        parset_dict['direction_specific']['ndir_selfcal'] = -1

    # Total number of directions to process (default = all). If this number is
    # greater than ndir_selfcal, then the remaining directions will not be selfcal-
    # ed but will instead be imaged with the selfcal solutions from the nearest
    # direction for which selfcal succeeded (if a target is specified and
    # target_has_own_facet = True, it will be imaged in this way after ndir_total
    # number of directions are processed)
    if 'ndir_total' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['ndir_total'] = parset.getint('directions', 'ndir_total')
        log.info("Processing up to %s direction(s) in total" % (parset_dict['direction_specific']['ndir_total']))
    else:
        parset_dict['direction_specific']['ndir_total'] = -1

    # A target can be specified to ensure that it falls entirely within a single
    # facet. The values should be those of a circular region that encloses the
    # source and not those of the target itself. Lastly, the target can be placed in
    # a facet of its own. In this case, it will not go through selfcal but will
    # instead use the selfcal solutions of the nearest facet for which selfcal was
    # done
    if 'target_radius_arcmin' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['target_radius_arcmin'] = parset.getfloat('directions',
            'target_radius_arcmin')
    if 'target_has_own_facet' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['target_has_own_facet'] = parset.getboolean('directions',
            'target_has_own_facet')
    else:
        parset_dict['direction_specific']['target_has_own_facet'] = False

    # Radius in degrees within which the direction-dependent solutions will be
    # transferred before starting selfcal (default = 0; i.e., disabled). If a
    # direction is within this distance of a calibrator for which selfcal was
    # successful, the dir-dep selfcal solutions from this calibrator will be used
    # instead of the dir-indep ones
    if 'transfer_radius' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['transfer_radius'] = parset.getfloat('directions',
            'transfer_radius')
    else:
        parset_dict['direction_specific']['transfer_radius'] = 0.0

    # Check for unused options
    given_options = parset.options('directions')
    allowed_options = ['directions_file', 'flux_min_Jy', 'size_max_arcmin',
        'separation_max_arcmin', 'max_num', 'check_edges', 'ndir_total',
        'ndir_selfcal', 'transfer_radius', 'groupings', 'allow_reordering',
        'target_ra', 'target_dec', 'target_radius_arcmin', 'target_has_own_facet']
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
        parset_dict = {'cluster_specific': parset._sections['cluster']}
    else:
        parset_dict = {'cluster_specific': {}}

    # Full path to a local disk on the nodes for I/O-intensive processing. The path
    # must be the same for all nodes. If not given, the default directory in the
    # working directory is used
    if 'ncpu' in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['ncpu'] = parset.getint('cluster', 'ncpu')
    else:
        import multiprocessing
        parset_dict['cluster_specific']['ncpu'] = multiprocessing.cpu_count()
    log.info("Using up to %s CPU(s) per node" % (parset_dict['cluster_specific']['ncpu']))

    # Maximum fraction of the total memory per node that WSClean may use (default =
    # 0.9)
    if 'fmem' in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['fmem'] = parset.getfloat('cluster', 'fmem')
        if parset_dict['cluster_specific']['fmem'] > 1.0:
            parset_dict['cluster_specific']['fmem'] = 1.0
    else:
        parset_dict['cluster_specific']['fmem'] = 0.9
    log.info("Using up to {0}% of the memory per node for WSClean".format(parset_dict['cluster_specific']['fmem']*100.0))

    # Number of directions to process in parallel on each node (default = 1). If
    # directions are split into groups to be processed in parallel (with the
    # groupings parameter), this parameter controls how many directions are run
    # simultaneously on a single node. Note that the number of CPUs (set with the
    # ncpu parameter) will be divided among the directions on each node
    if 'ndir_per_node' in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['ndir_per_node'] = parset.getint('cluster',
            'ndir_per_node')
    else:
        parset_dict['cluster_specific']['ndir_per_node'] = 1
    log.info("Processing up to %s direction(s) in parallel per node" %
        (parset_dict['cluster_specific']['ndir_per_node']))

    # Number of imager jobs to run per node (affects initsubtract and facetimage
    # operationa; default = 1). If your nodes have many CPUs and > 32 GB of memory,
    # it may be advantageous to set this to 2 or more
    if 'nimg_per_node' in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['nimg_per_node'] = parset.getint('cluster',
            'nimg_per_node')
    else:
        parset_dict['cluster_specific']['nimg_per_node'] = 1

    # Full path to cluster description file. Use clusterdesc_file = PBS to use the
    # PBS / torque reserved nodes. If not given, the clusterdesc file for a single
    # (i.e., local) node is used
    if 'clusterdesc_file' not in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['clusterdesc_file'] = parset_dict['lofarroot'] + '/share/local.clusterdesc'
        parset_dict['cluster_specific']['node_list'] = ['localhost']

    # Full path to a local disk on the nodes for I/O-intensive processing. The path
    # must be the same for all nodes. If not given, the default directory in the
    # working directory is used
    if 'dir_local' not in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['dir_local'] = None

    # Check for unused options
    given_options = parset.options('cluster')
    allowed_options = ['ncpu', 'fmem', 'ndir_per_node', 'nimg_per_node',
        'clusterdesc_file', 'dir_local']
    for option in given_options:
        if option not in allowed_options:
            log.warning('Option "{}" was given in the [cluster] section of the '
                'parset but is not a valid cluster option'.format(option))

    return parset_dict


def get_ms_options(parset, ms_files):
    """
    Handle the ms-specific options

    Parameters
    ----------
    parset : RawConfigParser object
        Input parset
    ms_files : list of str
        MS filenames

    Returns
    -------
    parset_dict : dict
        Dictionary with all ms-specific options

    """
    parset_dict = {'ms_specific': {}}
    for ms in ms_files:
        ms = os.path.basename(ms)
        if not ms in parset._sections.keys():
            continue
        parset_dict['ms_specific'][ms] = parset._sections[ms]

    return parset_dict
