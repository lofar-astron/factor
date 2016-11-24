#!/usr/bin/env python
"""
Script to merge selfcal parmdbs
"""
import argparse
from argparse import RawTextHelpFormatter
import os
import lofar.parmdb as pdb
import casacore.tables as pt
import shutil
import numpy as np


def main(parmdb_p, parmdb_a, parmdb_out, clobber=True):
    """
    Merges facet selfcal parmdbs into a single parmdb

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
            shutil.rmtree(parmdb_out)
        else:
            return

    os.system('cp -r {0} {1}'.format(parmdb_p,parmdb_out))

    ## Copy over the Gains
    pdb_out = pdb.parmdb(parmdb_out)
    pdb_a = pdb.parmdb(parmdb_a)
    parms = pdb_a.getValuesGrid('*')
    for parmname in pdb_a.getNames():
        # Set flagged solutions to NaN
        flagged = np.where(np.logical_or(parms[parmname]['values'] == 0.0,
            np.isnan(parms[parmname]['values'])))
        parms[parmname]['values'][flagged] = np.nan
        ValueHolder = pdb_out.makeValue(values=parms[parmname]['values'],
                                        sfreq=parms[parmname]['freqs'],
                                        efreq=parms[parmname]['freqwidths'],
                                        stime=parms[parmname]['times'],
                                        etime=parms[parmname]['timewidths'],
                                        asStartEnd=False)

        pdb_out.addValues(parmname,ValueHolder)
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
