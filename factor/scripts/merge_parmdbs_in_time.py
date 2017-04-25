#!/usr/bin/env python
"""
Script to merge parmdbs in time
"""
import argparse
from argparse import RawTextHelpFormatter
import os
import casacore.tables as pt
import shutil
import lofar.parmdb as pdb
import sys
import numpy as np


def main(input_mslist, parmdb_name, outparmdb, clobber=True, scratch_dir=None):
    """
    Merges parmdbs in time into a single parmdb

    The parmdbs are assumed to be located in the input MS with name
    parmdb_name

    Parameters
    ----------
    input_mslist : list
        List of input MS file names
    parmdb_name : str
        Name of parmdb (relative to MS files)
    outparmdb : str
        Name of output merged parmdb
    clobber : bool, optional
        If True, overwrite existing output file
    scratch_dir : str, optional
        Scratch directory for temp storage

    """
    if type(input_mslist) is str:
        input_mslist = input_mslist.strip('[]').split(',')
        input_mslist = [f.strip() for f in input_mslist]
    inparmdbs = [os.path.join(ms, parmdb_name) for ms in input_mslist]

    if type(clobber) is str:
        if clobber.lower() == 'true':
            clobber = True
        else:
            clobber = False

    if os.path.exists(outparmdb):
        if clobber:
            os.system('rm -rf {0}'.format(outparmdb))
        else:
            return

    # Copy to scratch directory if specified
    if scratch_dir is not None:
        inparmdbs_orig = inparmdbs
        inparmdbs = [os.path.join(scratch_dir, os.path.basename(inp)+'_{}'.format(i))
            for i, inp in enumerate(inparmdbs_orig)]
        outparmdb_orig = outparmdb
        outparmdb = os.path.join(scratch_dir, os.path.basename(outparmdb_orig))
        for inp_orig, inp in zip(inparmdbs_orig, inparmdbs):
            shutil.copytree(inp_orig, inp)

    pdb_concat = pdb.parmdb(outparmdb, create=True)

    for i, inparmdb in enumerate(inparmdbs):
        pdb_add = pdb.parmdb(inparmdb)
        parms = pdb_add.getValuesGrid('*')
        for parmname in pdb_add.getNames():
            flagged = np.where(np.logical_or(parms[parmname]['values'] == 0.0,
                np.isnan(parms[parmname]['values'])))
            parms[parmname]['values'][flagged] = np.nan
            ValueHolder = pdb_concat.makeValue(values=parms[parmname]['values'],
                                               sfreq=parms[parmname]['freqs'],
                                               efreq=parms[parmname]['freqwidths'],
                                               stime=parms[parmname]['times'],
                                               etime=parms[parmname]['timewidths'],
                                               asStartEnd=False)
            pdb_concat.addValues(parmname, ValueHolder)
        pdb_concat.flush()
        pdb_add = False
    pdb_concat = False

    # Copy output to original path and delete copies if scratch directory is specified
    if scratch_dir is not None:
        shutil.copytree(outparmdb, outparmdb_orig)
        shutil.rmtree(outparmdb)
        for inp in inparmdbs:
            shutil.rmtree(inp)


if __name__ == '__main__':
    descriptiontext = "Merge parmdbs in time.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('mslist', help='list of ms with parmdbs to merge')
    parser.add_argument('parmdb_name', help='name of parmdbs to merge')
    parser.add_argument('outparmdb', help='output parmdb')
    parser.add_argument('-c', '--clobber', help='clobber existing file?', type=bool, default=True)

    args = parser.parse_args()
    main(args.mslist, args.parmdb_name, args.outparmdb, clobber=args.clobber)
