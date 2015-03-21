"""
Module that holds all parset-related functions
"""
import sys
import os
import glob
import logging
import ConfigParser
from factor._logging import set_log_file

log = logging.getLogger('parset')


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
                "specify it in the parset")
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
                "Please specify it in the parset")
            sys.exit(1)

    # Set-up the working dir (other paths are relative to this)
    if not os.path.isdir(parset_dict['dir_working']):
        os.mkdir(parset_dict['dir_working'])
    try:
        os.chdir(parset_dict['dir_working'])
        for subdir in ['/logs', '/images', '/models', '/state', '/parsets',
            '/visdata', '/parmdbs', '/pipeline', '/datamaps', '/regions']:
            if not os.path.isdir(parset_dict['dir_working']+subdir):
                os.mkdir(parset_dict['dir_working']+subdir)
    except:
        log.critical("Cannot use the working dir %s" % (parset_dict['dir_working']))
        sys.exit(1)
    set_log_file(parset_dict['dir_working']+'/factor.log')
    log.info("Working directory is {0}".format(parset_dict['dir_working']))

    # Get all the MS files in the input directory
    parset_dict['mss'] = glob.glob(parset_dict['dir_ms']+'/*[MS|ms]')
    if len(parset_dict['mss']) == 0:
        log.error('No MS files found in {0}!'.format(parset_dict['dir_ms']))
        sys.exit(1)
    log.info("Working on %i band(s)" % (len(parset_dict['mss'])))

    # Get dir-indep instrument parmdbs
    if 'parmdb_name' not in parset_dict:
        parset_dict['parmdb_name'] = 'instrument'

    # Some check on types and defaults
    if 'interactive' in parset_dict:
        parset_dict['interactive'] = parset.getboolean('global', 'interactive')
        log.debug("Using interactive mode")
    else:
        parset_dict['interactive'] = False
    if 'make_mosaic' in parset_dict:
        parset_dict['make_mosaic'] = parset.getboolean('global', 'make_mosaic')
        log.debug("Making final mosaic")
    else:
        parset_dict['make_mosaic'] = False
    if 'use_ftw' in parset_dict:
        parset_dict['use_ftw'] = parset.getboolean('global', 'use_ftw')
        log.debug("Using FT / FTW")
    else:
        parset_dict['use_ftw'] = True
    if 'imager' not in parset_dict:
        parset_dict['imager'] = 'awimager'
    if parset_dict['imager'].lower() not in ['awimager', 'casapy', 'wsclean']:
        log.error('Imager "{0}" not understood'.format(parset_dict['imager']))
        sys.exit(1)
    if 'imagerroot' not in parset_dict:
        parset_dict['imagerroot'] = parset_dict['lofarroot']

    # Handle directions-related parameters
    if 'directions' in parset._sections.keys():
        parset_dict['direction_specific'] = parset._sections['directions']
    else:
        parset_dict['direction_specific'] = {}
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
        parset_dict['direction_specific']['one_at_a_time'] = False
        log.debug("Using the following groupings for directions: {0}".format(groupings))
    else:
        parset_dict['direction_specific']['one_at_a_time'] = True
    if 'ndir' in parset_dict['direction_specific']:
        parset_dict['direction_specific']['ndir'] = parset.getint('directions', 'ndir')
        log.debug("Processing up to %s directions in total" % (parset_dict['direction_specific']['ndir']))
    else:
        parset_dict['direction_specific']['groupings'] = {'1': 0}
        parset_dict['direction_specific']['ndir'] = -1

    # Handle cluster-related parameters
    if 'cluster' in parset._sections.keys():
        parset_dict['cluster_specific'] = parset._sections['cluster']
    else:
        parset_dict['cluster_specific'] = {}
    if 'ncpu' in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['ncpu'] = parset.getint('cluster', 'ncpu')
        log.debug("Using up to %s CPUs per node" % (parset_dict['cluster_specific']['ncpu']))
    else:
        parset_dict['cluster_specific']['ncpu'] = 1
    if 'clusterdesc_file' not in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['clusterdesc_file'] = parset_dict['lofarroot'] + '/share/local.clusterdesc'
        parset_dict['cluster_specific']['node_list'] = ['localhost']
    if 'distribute' in parset_dict['cluster_specific']:
        parset_dict['cluster_specific']['distribute'] = parset.getboolean('cluster', 'distribute')
    else:
        parset_dict['cluster_specific']['distribute'] = False

    # load MS-specific parameters
    parset_dict['ms_specific'] = {}
    for ms in parset_dict['mss']:
        ms = os.path.basename(ms)
        if not ms in parset._sections.keys():
            continue
        paset_dict_ms = parset._sections[ ms ]
        parset_dict['ms_specific'][ms] = paset_dict_ms

    return parset_dict
