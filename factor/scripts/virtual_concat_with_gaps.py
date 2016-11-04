#!/usr/bin/env python
"""
Script to perform a virtual concatenation with frequecy gaps
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables as pt
import numpy as np
import os
import sys
import uuid


def main(ms_files, outfile, clobber=True):
    """
    Performs a virtual concatenation with possible frequency gaps

    Parameters
    ----------
    ms_files : list
        List of files to merge, ordered by frequency. Files that do not exist
        are identified as gaps and are filled with flagged dummy data
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
            pt.tabledelete(outfile)
        else:
            return

    # Find at least one existing ms
    ms_exists = None
    for ms in ms_files:
        if os.path.exists(ms):
            ms_exists = ms
            sw = pt.table('{}::SPECTRAL_WINDOW'.format(ms))
            ms_exists_ref_freq = sw.getcol('REF_FREQUENCY')[0]
            sw.close()
            break
    if ms_exists is None:
        print('ERROR: no files exist')
        sys.exit(1)

    # Identify gaps
    ms_files_to_concat = []
    for i, ms in enumerate(ms_files):
        if not os.path.exists(ms):
            # Missing file means gap, so create an appropriate dummy dataset with
            # a random name
            ms_new = '{0}_{1}.ms'.format(os.path.splitext(ms)[0], uuid.uuid4().urn.split('-')[-1])
            pt.tableutil.tablecopy(ms_exists, ms_new)

            # Alter SPECTRAL_WINDOW subtable as appropriate to fill gap
            sw = pt.table('{}::SPECTRAL_WINDOW'.format(ms_new), readonly=False)
            tot_bandwidth = sw.getcol('TOTAL_BANDWIDTH')[0]
            if i > 0:
                sw_low = pt.table('{}::SPECTRAL_WINDOW'.format(ms_files[i-1]))
                ref_freq = sw_low.getcol('REF_FREQUENCY') + tot_bandwidth
                sw_low.close()
            else:
                for j in range(1, len(ms_files)-1):
                    if os.path.exists(ms_files[j]):
                        sw_high = pt.table('{}::SPECTRAL_WINDOW'.format(ms_files[j]))
                        ref_freq = sw_high.getcol('REF_FREQUENCY') - tot_bandwidth * j
                        sw_high.close()
                        break
            chan_freq = sw.getcol('CHAN_FREQ') - ms_exists_ref_freq + ref_freq
            sw.putcol('REF_FREQUENCY', ref_freq)
            sw.putcol('CHAN_FREQ', chan_freq)
            sw.close()

            # Flag all data
            t = pt.table(ms_new, readonly=False)
            t.putcol('FLAG_ROW', np.ones(len(t), dtype=bool))
            t.close()

            ms_files_to_concat.append(ms_new)
        else:
            ms_files_to_concat.append(ms)

    # Concat
    pt.msutil.msconcat(ms_files_to_concat, outfile)


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
