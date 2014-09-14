# This module handle the reading of the parset

import sys, os, glob
import logging
import ConfigParser

def parset_read(parset_file):

    if not os.path.isfile(parset_file):
        logging.critical("Missing parset file (%s), I don't know what to do :'(" % (parset_file))
        sys.exit(1)

    logging.info("Reading parset file: %s" % (parset_file))

    parset = ConfigParser.RawConfigParser()
    parset.read(parset_file)

    parset_dict = parset._sections['global']

    # get all the MS in the directory
    parset_dict['mss'] = glob.glob(parset_ditc['ms_dir']+'/*MS')

    return parset_dict
