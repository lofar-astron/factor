#!/usr/bin/env python
"""
Script to finalize output for LoSoTo
"""
import os
import sys

if __name__ == '__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog <parmdb_filename> <ms_filename> <solset>')
    (options, args) = opt.parse_args()

    if len(args) != 3:
        opt.print_help()
        sys.exit(1)

    parmdb_filename = args[0]
    ms_filename = args[1]
    solset = args[2]

    if os.path.exists('{0}/instrument'.format(ms_filename)):
        os.system('rm -rf {0}/instrument'.format(ms_filename))
    os.system('cp -r {0} {1}/instrument'.format(parmdb_filename, ms_filename))

    if os.path.exists(parmdb_filename):
        os.system('rm -rf {0}'.format(parmdb_filename))
    os.system('cp -r {0}/{1}_instrument {2}'.format(ms_filename, solset, parmdb_filename))
