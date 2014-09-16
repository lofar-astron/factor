# This module handle the reading of the parset

import sys, os, glob
import logging
import ConfigParser

log = logging.getLogger('parset')

def parset_read(parset_file):

    if not os.path.isfile(parset_file):
        log.critical("Missing parset file (%s), I don't know what to do :'(" % (parset_file))
        sys.exit(1)

    log.info("Reading parset file: %s." % (parset_file))

    parset = ConfigParser.RawConfigParser()
    parset.read(parset_file)

    parset_dict = parset._sections['global']

    # set-up the working dir (first thing to do, other path may be relative to this)
    try:
        os.chdir(parset_dict['dir_working'])
        if not os.path.isdir(parset_dict['dir_working']+'/log'): os.mkdir(parset_dict['dir_working']+'/log')
        if not os.path.isdir(parset_dict['dir_working']+'/img'): os.mkdir(parset_dict['dir_working']+'/img')
    except:
        log.critical("Cannot use the working dir %s." % (parset_dict['dir_working']))
        sys.exit(1)

    # get all the MS in the directory
    parset_dict['mss'] = glob.glob(parset_dict['dir_ms']+'/*[MS|ms]')
    log.info("Working on %i MSs" % (len(parset_dict['mss'])))

    # some check on types
    if 'ncpu' in parset_dict: parset_dict['ncpu'] = parset.getint('global', 'ncpu')
    log.debug("Using %i processors for multi-thread." % (parset_dict['ncpu']))

    # load MS-specific parameters
    parset_dict['ms_specific'] = {}
    for ms in parset_dict['mss']:
        ms = os.path.basename(ms)
        if not ms in parset._sections.keys(): continue
        paset_dict_ms = parset._sections[ ms ]
        parset_dict['ms_specific'][ms] = paset_dict_ms

    return parset_dict
