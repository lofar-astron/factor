# This module handle the reading of directions

import os
import numpy as np
import logging 
from factor.lib.direction import direction as d

log = logging.getLogger('parset')

def directions_read(directions_file):

    if not os.path.isfile(directions_file):
        log.critical("Missing directions file (%s)." % (directions_file))
        sys.exit(1)

    log.info("Reading directions file: %s" % (directions_file))
    types = np.dtype({'names':['name','ra','dec','reg','multiscale','solint_a','solint_p','hdr'], \
                    'formats':['S100',np.float,np.float,'S100',np.bool,np.int,np.int,np.bool]})
    directions = np.genfromtxt(directions_file, comments='#', delimiter=',', unpack=False,
                      converters={0: lambda x: x.strip(), 3: lambda x: x.strip()}, dtype=types)

    data = []
    for direction in directions:
        data.append( d(*direction) )
    
    return data
