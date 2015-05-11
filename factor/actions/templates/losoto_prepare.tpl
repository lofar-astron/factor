#!/usr/bin/env python
"""
Script to prepare input for LoSoTo
"""
import os
import sys

if __name__ == '__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog <parmdb_filename> <ms_filename>')
    (options, args) = opt.parse_args()

    if len(args) != 2:
        opt.print_help()
        sys.exit(1)

    parmdb_filename = args[0]
    ms_filename = args[1]

    if os.path.exists('{0}/instrument'.format(ms_filename)):
        os.system('rm -rf {0}/instrument'.format(ms_filename))
    os.system('cp -r {0} {1}/instrument'.format(parmdb_filename, ms_filename))
