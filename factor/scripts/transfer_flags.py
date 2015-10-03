#! /usr/bin/env python
"""
Script to tranfer flags from one MS file to another
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import sys
import numpy as np


def main(ms1, ms2):
    """
    Transfer flags

    Parameters
    ----------
    ms1 : str
        Name of MS file from which the flags will be taken
    ms2 : str or list
        Name of MS file or list of files to which the flags will be transferred

    """
    if type(ms2) is str:
        if '[' in ms2:
            ms2 = ms2.strip('[]').split(',')
            ms2 = [m.strip() for m in ms2]
        else:
            ms2 = [ms2]

    t1 = pt.table(ms1, readonly=True, ack=False)
    flags1 = t1.getcol('FLAG')

    numberofchans1 = numpy.int(numpy.shape(flags1)[1])
    chanperms = numpy.int(len(mslist))/numberofchans1

    for ms_id, ms in enumerate(ms2):
        if os.path.isdir(ms):
            flagsin = flags1.getcolslice('FLAG', [chanperms*ms_id, 0], [(chanperms*(ms_id+1))-1, 3])

            t2 = pt.table(ms, readonly=False, ack=False)
            flags2 = t2.getcol('FLAG')

            # Expand flags to match output MS
            numberofchans2 = numpy.int(numpy.shape(flags2)[1])
            flagsout = np.repeat(flagsin, numberofchans2, axis=0)

            # Perform logical OR to pick up ms flags
            flagsout = np.logical_or(flagsout, flags2)

            # Write data
            t2.putcol('FLAG', flagsout)
            t2.flush()
            t2.close()
    t1.close()


if __name__ == '__main__':
    descriptiontext = "Transfer flags from ms1 to ms2.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms1', help='name of MS file 1')
    parser.add_argument('ms2', help='name of MS file 2')
    args = parser.parse_args()

    main(args.ms1, args.ms2)





