#!/usr/bin/env python
"""
Script to merge parmdbs in time
"""
import argparse
from argparse import RawTextHelpFormatter
import os
import pyrap.tables as pt
import lofar.parmdb as pdb
import sys


def main(input_mslist, parmdb_name, outparmdb, clobber=True):
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
    os.system('cp -r {0} {1}'.format(inparmdbs[0], outparmdb))

### old copy, which puts all values into the same domain (i.e. ignoring the timestamps)
#    if len(inparmdbs) > 1:
#        pdb_concat = lofar.parmdb.parmdb(outparmdb)
#        for inparmdb in inparmdbs[1:]:
#            pdb = lofar.parmdb.parmdb(inparmdb)
#            for parmname in pdb.getNames():
#                v = pdb.getValuesGrid(parmname)
#                pdb_concat.addValues(v.copy())
#        pdb_concat.flush()

### version with TAQL, which avoids the previous problem, but may introduce new ones
#    if len(inparmdbs) > 1:
#        for inparmdb in inparmdbs[1:]:
#            pt.taql('insert into '+outparmdb+' select from '+inparmdb)
#    outdb = pt.table(outparmdb,readonly=False,ack=False)
#    errdesc = outdb.getcoldesc('ERRORS')
#    errdesc['name'] = 'ERRORS'
#    outdb.removecols('ERRORS')
#    outdb.addcols(errdesc)
#    outdb.close()

### version which adds new data by generating a new data-holder first
    if len(inparmdbs) > 1:
        pdb_concat = pdb.parmdb(outparmdb)
        for inparmdb in inparmdbs[1:]:
            pdb_add = pdb.parmdb(inparmdb)
            for parmname in pdb_add.getNames():
                parms = pdb_add.getValuesGrid(parmname)
                ValueHolder = pdb_concat.makeValue(values=parms[parmname]['values'],
                                                   sfreq=parms[parmname]['freqs'],
                                                   efreq=parms[parmname]['freqwidths'],
                                                   stime=parms[parmname]['times'],
                                                   etime=parms[parmname]['timewidths'],
                                                   asStartEnd=False)
                pdb_concat.addValues(parmname,ValueHolder)       
            pdb_concat.flush()
            pdb_add = False 
        pdb_concat = False


if __name__ == '__main__':
    descriptiontext = "Merge parmdbs in time.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('mslist', help='list of ms with parmdbs to merge')
    parser.add_argument('parmdb_name', help='name of parmdbs to merge')
    parser.add_argument('outparmdb', help='output parmdb')
    parser.add_argument('-c', '--clobber', help='clobber existing file?', type=bool, default=True)

    args = parser.parse_args()
    main(args.mslist, args.parmdb_name, args.outparmdb, clobber=args.clobber)
