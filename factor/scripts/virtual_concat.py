#!/usr/bin/env python
"""
Script to perform a virtual concatenation
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import os


def main(ms_files, outfile, clobber=True):
    """
    Performs a virtual concatenation

    Parameters
    ----------
    ms_files : list
        List of files to merge
    outfile : str
        Output file
    clobber : bool, optional
        If True, existing files are overwritten

    """
    if type(ms_files) is str:
        ms_files = [f.strip(' \'\"') for f in ms_files.strip('[]').split(',')]
    if type(clobber) is str:
        if clobber.lower() == 'true':
            clobber = True
        else:
            clobber = False
    if os.path.exists(outfile):
        if clobber:
            os.system('rm -rf {0}'.format(outfile))
        else:
            return

    # TODO: order ms_files by time?
    pt.msutil.msconcat(ms_files, outfile)


if __name__ == '__main__':
    descriptiontext = "Perform virtual concatenation.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms_files', nargs='+', help='list of ms files to concatenate')
    parser.add_argument('outfile', help='output filename')
    parser.add_argument('-c', '--clobber', help='overwrite existing outfile?', type=bool, default=True)

    args = parser.parse_args()
    print "ms_files:",args.ms_files
    print "outfile:",args.outfile
    main(args.ms_files, args.outfile, clobber=args.clobber)
