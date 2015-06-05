#! /usr/bin/env python
"""
Script to split a dataset into chunks of equal time
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import numpy as np
import sys
import os


def main(dataset, blockl, clobber=True):
    """
    Split dataset into time chunks of length blockl time slots

    Parameters
    ----------
    dataset : str
        Name of MS file to split
    blockl : int
        Number of time slots per chunk
    clobber : bool, optional
        If True, existing files are overwritten

    """
    if type(clobber) is str:
        if clobber.lower() == 'true':
            clobber = True
        else:
            clobber = False

    blockl = int(blockl)
    if blockl < 1:
        blockl = 1

    # Get time per sample and number of samples
    t = pt.table(dataset, readonly=True, ack=False)
    for t2 in t.iter(["ANTENNA1","ANTENNA2"]):
        if (t2.getcell('ANTENNA1',0)) < (t2.getcell('ANTENNA2',0)):
            timepersample = t2[1]['TIME']-t2[0]['TIME'] # sec
            nsamples = t2.nrows()
            break
    t.close()

    nchunks = int(np.ceil((np.float(nsamples) / np.float(blockl))))
    tlen = timepersample * np.float(blockl) / 3600.0 # length of block in hours
    tobs = timepersample * nsamples / 3600.0 # length of obs in hours

    files = []
    for c in range(nchunks):
        chunk_file = '{0}_chunk{1}.ms'.format(os.path.splitext(dataset)[0], c)
        files.append(chunk_file)
        t0 = tlen * np.float(c) # hours
        t1 = t0 + tlen # hours
        if c == 0:
            t0 = -0.1 # make sure first chunk gets first slot
        if c == nchunks-1 and t1 < tobs:
            t1 = tobs + 0.1 # make sure last chunk gets all that remains
        split_ms(dataset, chunk_file, t0, t1, clobber=clobber)

    return {'files': files}


def split_ms(msin, msout, start_out, end_out, clobber=True):
    """
    Splits an MS between start and end times in hours relative to first time

    Parameters
    ----------
    msin : str
        Name of MS file to split
    msout : str
        Name of output MS file
    start_out : float
        Start time in hours relative to first time
    end_out : float
        End time in hours relative to first time
    clobber : bool, optional
        If True, existing files are overwritten

    """
    if os.path.exists(msout):
        if clobber:
            os.system('rm -rf {0}'.format(msout))
        else:
            return

    t = pt.table(msin, ack=False)
    starttime = t[0]['TIME']

    t1 = t.query('TIME >= ' + str(starttime+start_out*3600.0) + ' && '
      'TIME < ' + str(starttime+end_out*3600.0), sortlist='TIME,ANTENNA1,ANTENNA2')

    t1.copy(msout, True)
    t1.close()
    t.close()


if __name__ == '__main__':
    descriptiontext = "Chunk a dataset in time.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms_filename', help='Dataset name')
    parser.add_argument('-w', '--width', help='width of chunks in number of samples', type=int, default=10)

    args = parser.parse_args()
    main(args.ms_filename, blockl=args.width)
