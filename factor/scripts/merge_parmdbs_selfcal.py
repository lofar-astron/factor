#!/usr/bin/env python
"""
Script to merge parmdbs in time
"""
import argparse
from argparse import RawTextHelpFormatter
import os
import lofar.parmdb
import sys


def main(parmdb_p, parmdb_a, parmdb_out, clobber=True):
    """
    Merges facet selfcal parmdbs into a parmdb for a single band

    Parameters
    ----------
    parmdb_p : str
        Filename of CommonScalarPhase and TEC parmdb
    parmdb_a : str
        Filename of Gain parmdb. The nearset match in frequency to that of the
        input band will be used
    parmdb_out : str
        Filename of output file
    clobber : bool, optional
        If True, overwrite existing output file

    """
    if type(clobber) is str:
        if clobber.lower() == 'true':
            clobber = True
        else:
            clobber = False

    if os.path.exists(parmdb_out):
        if clobber:
            os.system('rm -rf {0}'.format(parmdb_out))
        else:
            return parmdb_out
    pdb_out = lofar.parmdb.parmdb(parmdb_out, create=True)

    # Copy over the CommonScalar phases and TEC
    pdb_p = lofar.parmdb.parmdb(parmdb_p)
    for parmname in pdb_p.getNames():
        parms = pdb_p.getValuesGrid(parmname)
        pdb_out.addValues(parms)

    # Copy over the Gains
    pdb_a = lofar.parmdb.parmdb(parmdb_a)
    for parmname in pdb_a.getNames():
        parms = pdb_a.getValuesGrid(parmname)
        pdb_out.addValues(parms)

    # Write values
    pdb_out.flush()


if __name__ == '__main__':
    descriptiontext = "Merge parmdbs in time.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('parmdb_p', help='phase parmdb')
    parser.add_argument('parmdb_a', help='gain parmdb')
    parser.add_argument('parmdb_out', help='output parmdb')
    parser.add_argument('-c', '--clobber', help='clobber existing file?', type=bool, default=True)

    args = parser.parse_args()
    main(args.parmdb_p, args.parmdb_a, args.parmdb_out, clobber=args.clobber)
