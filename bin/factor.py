#!/usr/bin/env python
# encoding: utf-8
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Authors:
# Francesco de Gasperin
# Jose Sabater Montes
# Wendy Williams
# Martin Hardcastle
#
# The code is based on an original idea of Reinout van Weeren

_author = "Francesco de Gasperin (fdg@hs.uni-hamburg.de)"

import sys, os, time
import logging
import factor

if __name__=='__main__':

    import optparse
    opt = optparse.OptionParser(usage='%prog [-v|-q] parset [default: facotr.parset] \n'
            +_author, version='%prog '+_version.__version__)
    opt.add_option('-q', help='Quiet', action='store_true', default=False)
    opt.add_option('-v', help='Verbose', action='store_true', default=False)
    (options, args) = opt.parse_args()

    if options.q:
        _logging.set_level('warning')
    if options.v:
        _logging.set_level('debug')

    # Check options
    if len(args) not in [0, 1]:
        opt.print_help()
        sys.exit()

    try: parset_file = args[1]
    except: parset_file = 'factor.parset'

    if not os.path.isfile(parset_file) and not options.i:
        logging.critical("Missing parset file, I don't know what to do :'(")
        sys.exit(1)

    parset = factor.parset.parset_read( parset_file )

    for step in steps:
        operation = parset.getString( '.'.join( [ "LoSoTo.Steps", step, "Operation" ] ) )
        if not operation in operations:
            logging.error('Unkown operation: '+operation)
            continue
        logging.info("--> Starting \'" + step + "\' step (operation: " + operation + ").")
        start = time.clock()
        returncode = operations[ operation ].run( step, parset, H )
        if returncode != 0:
           logging.error("Step \'" + step + "\' incomplete. Try to continue anyway.")
        else:
           logging.info("Step \'" + step + "\' completed successfully.")
        elapsed = (time.clock() - start)
        logging.debug("Time for this step: "+str(elapsed)+" s.")

    del H
    logging.info("Done.")

