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

    parset_dict = parset._sections['global']

    # Get LOFAR and Python paths
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

    # Set-up the working dir (other paths are relative to this)
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

    # Set default dir-indep instrument parmdb name
    if 'parmdb_name' not in parset_dict:
        parset_dict['parmdb_name'] = 'instrument'

    # Some check on types and defaults
    if 'interactive' in parset_dict:
        parset_dict['interactive'] = parset.getboolean('global', 'interactive')
    else:
        parset_dict['interactive'] = False
    if 'make_mosaic' in parset_dict:
        parset_dict['make_mosaic'] = parset.getboolean('global', 'make_mosaic')
    else:
        parset_dict['make_mosaic'] = True
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
    if 'use_chgcentre' in parset_dict:
        parset_dict['use_chgcentre'] = parset.getboolean('global', 'use_chgcentre')
    else:
        parset_dict['use_chgcentre'] = False
    if 'wsclean_nbands' in parset_dict:
        parset_dict['wsclean_nbands'] = parset.getint('global', 'wsclean_nbands')
    else:
        parset_dict['wsclean_nbands'] = 3
    if 'facet_imager' not in parset_dict:
        parset_dict['facet_imager'] = 'wsclean'
    if 'keep_avg_facet_data' in parset_dict:
        parset_dict['keep_avg_facet_data'] = parset.getboolean('global', 'keep_avg_facet_data')
    else:
        parset_dict['keep_avg_facet_data'] = True
    if 'keep_unavg_facet_data' in parset_dict:
        parset_dict['keep_unavg_facet_data'] = parset.getboolean('global', 'keep_unavg_facet_data')
    else:
        parset_dict['keep_unavg_facet_data'] = False
    if 'max_selfcal_loops' in parset_dict:
        parset_dict['max_selfcal_loops'] = parset.getint('global', 'max_selfcal_loops')
    else:
        parset_dict['max_selfcal_loops'] = 10
    if 'preaverage' in parset_dict:
        parset_dict['preaverage'] = parset.getboolean('global', 'preaverage')
    else:
        parset_dict['preaverage'] = True

    # Handle directions-related parameters
    if 'directions' in parset._sections.keys():
        parset_dict['direction_specific'] = parset._sections['directions']
    else:
        parset_dict['direction_specific'] = {}
    if 'check_edges' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['check_edges'] = parset.getboolean('directions',
            'check_edges')
    else:
        parset_dict['direction_specific']['check_edges'] = False
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
    if 'allow_reordering' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['allow_reordering'] = parset.getboolean('directions',
            'allow_reordering')
    else:
        parset_dict['direction_specific']['allow_reordering'] = True
    if 'ndir_selfcal' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['ndir_selfcal'] = parset.getint('directions', 'ndir_selfcal')
        log.info("Self calibrating up to %s direction(s)" % (parset_dict['direction_specific']['ndir_selfcal']))
    else:
        parset_dict['direction_specific']['ndir_selfcal'] = -1
    if 'ndir_total' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['ndir_total'] = parset.getint('directions', 'ndir_total')
        log.info("Processing up to %s direction(s) in total" % (parset_dict['direction_specific']['ndir_total']))
    else:
        parset_dict['direction_specific']['ndir_total'] = -1
    if 'target_radius_arcmin' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['target_radius_arcmin'] = parset.getfloat('directions',
            'target_radius_arcmin')
    if 'target_has_own_facet' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['target_has_own_facet'] = parset.getboolean('directions',
            'target_has_own_facet')
    else:
        parset_dict['direction_specific']['target_has_own_facet'] = False
    if 'transfer_radius' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['transfer_radius'] = parset.getfloat('directions',
            'transfer_radius')
    else:
        parset_dict['direction_specific']['transfer_radius'] = 0.0

    # Handle cluster-related parameters
    if 'cluster' in parset._sections.keys():
        parset_dict['cluster_specific'] = parset._sections['cluster']
    else:
        parset_dict['cluster_specific'] = {}
    if 'ncpu' in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['ncpu'] = parset.getint('cluster', 'ncpu')
    else:
        import multiprocessing
        parset_dict['cluster_specific']['ncpu'] = multiprocessing.cpu_count()
    log.info("Using up to %s CPU(s) per node" % (parset_dict['cluster_specific']['ncpu']))
    if 'fmem' in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['fmem'] = parset.getfloat('cluster', 'fmem')
        if parset_dict['cluster_specific']['fmem'] > 1.0:
            parset_dict['cluster_specific']['fmem'] = 1.0
    else:
        parset_dict['cluster_specific']['fmem'] = 0.9
    log.info("Using up to {0}% of the memory per node for WSClean".format(parset_dict['cluster_specific']['fmem']*100.0))
    if 'ndir_per_node' in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['ndir_per_node'] = parset.getint('cluster',
            'ndir_per_node')
    else:
        parset_dict['cluster_specific']['ndir_per_node'] = 1
    log.info("Processing up to %s direction(s) in parallel per node" %
        (parset_dict['cluster_specific']['ndir_per_node']))
    if 'nimg_per_node' in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['nimg_per_node'] = parset.getint('cluster',
            'nimg_per_node')
    else:
        parset_dict['cluster_specific']['nimg_per_node'] = 1
    if 'clusterdesc_file' not in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['clusterdesc_file'] = parset_dict['lofarroot'] + '/share/local.clusterdesc'
        parset_dict['cluster_specific']['node_list'] = ['localhost']
    if 'dir_local' not in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['dir_local'] = None

    # load MS-specific parameters
    parset_dict['ms_specific'] = {}
    for ms in parset_dict['mss']:
        ms = os.path.basename(ms)
        if not ms in parset._sections.keys():
            continue
        paset_dict_ms = parset._sections[ ms ]
        parset_dict['ms_specific'][ms] = paset_dict_ms

    return parset_dict
