# This module handle the reading of the parset

import logging
import ConfigParser

def parset_read(parset_file):

    if not os.path.isfile(parset_file):
        logging.critical("Missing parset file (%s), I don't know what to do :'(" % (parset_file))
        sys.exit(1)

    logging.info("Reading parset file: %s" % (parset_file))

    config = ConfigParser.RawConfigParser()
    config.read(parset_file)

    return config._sections['global']
