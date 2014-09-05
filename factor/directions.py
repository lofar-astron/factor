# This module handle the reading of directions

import logging
import os

def directions_read(directions_file):

    if not os.path.isfile(directions_file):
        logging.critical("Missing directions file (%s)." % (directions_file))
        sys.exit(1)

    logging.info("Reading directions file: %s" % (directions_file))

    return ['dir', 'dir']
