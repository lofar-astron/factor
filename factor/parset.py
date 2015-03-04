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

    # set-up the working dir (other paths are relative to this)
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

    # get all the MS in the directory
    parset_dict['mss'] = glob.glob(parset_dict['dir_ms']+'/*[MS|ms]')
    log.info("Working on %i bands" % (len(parset_dict['mss'])))

    # some check on types
    if 'ndir' in parset_dict:
        parset_dict['ndir'] = parset.getint('global', 'ndir')
        log.debug("Processing up to %s directions in total" % (parset_dict['ndir']))
    else:
        parset_dict['ndir'] = -1
    if 'ndir_parallel' in parset_dict:
        parset_dict['ndir_parallel'] = parset.getint('global', 'ndir_parallel')
        log.debug("Processing up to {0} directions in parallel".format(
            parset_dict['ndir_parallel']))
    else:
        parset_dict['ndir_parallel'] = 1
    if 'ncpu' in parset_dict:
        parset_dict['ncpu'] = parset.getint('global', 'ncpu')
        log.debug("Using up to %s CPUs per node" % (parset_dict['ncpu']))
    else:
        parset_dict['ncpu'] = 1
    if 'interactive' in parset_dict:
        parset_dict['interactive'] = parset.getboolean('global', 'interactive')
        log.debug("Using interactive mode")
    else:
        parset_dict['interactive'] = False
    if 'use_ftw' in parset_dict:
        parset_dict['use_ftw'] = parset.getboolean('global', 'use_ftw')
        log.debug("Using ftw to FFT model image")
    else:
        parset_dict['use_ftw'] = False
    if 'directions_flux_min_jy' in parset_dict:
        parset_dict['directions_flux_min_jy'] = parset.getfloat('global',
            'directions_flux_min_jy')
    if 'directions_max_num' in parset_dict:
        parset_dict['directions_max_num'] = parset.getint('global',
            'directions_max_num')
    else:
        parset_dict['directions_max_num'] = None
    if 'directions_size_max_arcmin' in parset_dict:
        parset_dict['directions_size_max_arcmin'] = parset.getfloat('global',
            'directions_size_max_arcmin')
    if 'directions_separation_max_arcmin' in parset_dict:
        parset_dict['directions_separation_max_arcmin'] = parset.getfloat('global',
            'directions_separation_max_arcmin')

    # load MS-specific parameters
    parset_dict['ms_specific'] = {}
    for ms in parset_dict['mss']:
        ms = os.path.basename(ms)
        if not ms in parset._sections.keys():
            continue
        paset_dict_ms = parset._sections[ ms ]
        parset_dict['ms_specific'][ms] = paset_dict_ms

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


    return parset_dict
