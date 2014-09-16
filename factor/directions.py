# This module handle the reading of directions

import os
import logging 
from factor.lib.direction import direction as d

log = logging.getLogger('parset')

def directions_read(directions_file):

    if not os.path.isfile(directions_file):
        log.critical("Missing directions file (%s)." % (directions_file))
        sys.exit(1)

    log.info("Reading directions file: %s" % (directions_file))

    return [d('aaa'), d('bbb')]
