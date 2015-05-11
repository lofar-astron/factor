#!/usr/bin/env python
"""
Script to finalize output for BBS run
"""
import os
import sys

if __name__ == '__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog <ms_filename> [<parmdb_filename>]')
    (options, args) = opt.parse_args()

    if len(args) not in [1, 2]:
        opt.print_help()
        sys.exit(1)

    ms_filename = args[0]
    if ms_filename[-1] == '/':
        ms_filename = ms_filename[:-1]

    if len(args) > 1:
        # Copy 'ms_filename/instrument' to output parmdb if needed
        parmdb_filename = args[1]
        if parmdb_filename is not None:
            if not os.path.exists(parmdb_filename):
                instrument_filename = os.path.join(ms_filename, 'instrument')
                os.system('cp -r {0} {1}'.format(instrument_filename, parmdb_filename))
