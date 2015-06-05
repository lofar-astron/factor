#!/usr/bin/env python
"""
Script to merge parmdbs in time
"""
import argparse
from argparse import RawTextHelpFormatter
import os
import lofar.parmdb
import sys


def main(inparmdbs, outparmdb, clobber=True):
    """
    Merges parmdbs in time into a single parmdb

    Parameters
    ----------
    inparmdbs : list
        List of input parmdb file names
    outparmdb : str
        Name of output merged parmdb
    clobber : bool, optional
        If True, overwrite existing output file

    """
    if type(inparmdbs) is str:
        inparmdbs = inparmdbs.strip('[]').split(',')
        inparmdbs = [f.strip() for f in inparmdbs]

    if os.path.exists(outparmdb):
        if clobber:
            os.system('rm -rf {0}'.format(outparmdb))
        else:
            return outparmdb
    os.system('cp -r {0} {1}'.format(inparmdbs[0], outparmdb))

    if len(inparmdbs) > 1:
        pdb_concat = lofar.parmdb.parmdb(outparmdb)
        for inparmdb in inparmdbs[1:]:
            pdb = lofar.parmdb.parmdb(inparmdb)
            for parmname in pdb.getNames():
                v = pdb.getValuesGrid(parmname)
                pdb_concat.addValues(v.copy())
        pdb_concat.flush()


if __name__ == '__main__':
    descriptiontext = "Merge parmdbs in time.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('inparmdbs', help='list of parmdbs to merge')
    parser.add_argument('outparmdb', help='output parmdb')
    parser.add_argument('-c', '--clobber', help='clobber existing file?', type=bool, default=True)

    args = parser.parse_args()
    main(args.inparmdbs, args.outparmdb, clobber=args.clobber)
