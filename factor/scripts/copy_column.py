#! /usr/bin/env python
"""
Script to copy a column between MS files
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import numpy
import sys


def copy_column_to_ms(ms, inputcol, outputcol, ms_from=None):
    """
    Copies one column to another, within an MS file or between two MS files

    Parameters
    ----------
    ms : str
        MS file receiving copy
    inputcol : str
        Column name to copy from
    outputcol : str
        Column name to copy to
    ms_from : str, optional
        MS file to copy from. If None, the column is copied internally

    """
    t = pt.table(ms, readonly=False, ack=False)
    if ms_from is not None:
        tf = pt.table(ms_from, readonly=False, ack=False)
        data = tf.getcol(inputcol)
        desc = tf.getcoldesc(inputcol)
    else:
        data = t.getcol(inputcol)
        desc = t.getcoldesc(inputcol)

    # Add the output column if needed
    if outputcol not in t.colnames():
        desc['name'] = outputcol
        cd = pt.tableutil.makecoldesc(desc['name'], desc)
        tdesc = pt.tableutil.maketabdesc(cd)
        t._addcols(tdesc, {}, True)
        t._makerow()

    t.putcol(outputcol, data)
    t.flush()
    t.close()


def copy_column_to_bands(mslist, ms_from, inputcol, outputcol):
    """
    Copies one column from an MS file to multiple MS files (bands)

    Parameters
    ----------
    mslist : list
        MS files receiving copy
    ms_from : str
        MS file to copy from.
    inputcol : str
        Column name to copy from
    outputcol : str
        Column name to copy to

    """
    datain = pt.table(ms_from)
    data = datain.getcol(inputcol, nrow=1)
    numberofchans = numpy.int(numpy.shape(data)[1])
    chanperms = numberofchans/numpy.int(len(mslist))

    for ms_id, ms in enumerate(mslist):
        if os.path.isdir(ms):
            data = datain.getcolslice(inputcol, [chanperms*ms_id,0], [(chanperms*(ms_id+1))-1,3])
            dataout = pt.table(ms, readonly=False)
            dataout.putcol(outputcol, data)
            dataout.flush()
            dataout.close()


def main(ms_from, ms_to, column_from, column_to):
    """
    Copy a column between MS files

    Parameters
    ----------
    ms_from : str
        Name of MS file to copy from
    ms_to : str or list
        Name of MS file or list of MS files to copy to. May be given as a list
        or as a string (e.g., '[ms1, ms2]'
    column_from : str
        Name of column to copy from
    column_to : str
        Name of column to copy to

    """
    if type(ms_to) is str:
        if '[' in ms_to:
            ms_to = ms_to.strip('[]').split(',')
            ms_to = [m.strip() for m in ms_to]

    if type(ms_to) is list:
        # List means call copy_column_to_bands()
        copy_column_to_bands(ms_to, ms_from, column_from, column_to)
    else:
        if ms_to == ms_from:
            ms_from = None
        copy_column_to_ms(ms_to, column_from, column_to, ms_from)


if __name__ == '__main__':
    descriptiontext = "Copy a column between MS files.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms_from', help='name of the MS file from copy from')
    parser.add_argument('ms_to', help='name of the MS file from copy to')
    parser.add_argument('column_from', help='name of the column to copy from')
    parser.add_argument('column_to', help='name of the column to copy to')
    args = parser.parse_args()

    main(args.ms_from, args.ms_to, args.column_from, args.column_to)





