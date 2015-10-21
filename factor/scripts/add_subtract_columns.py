#! /usr/bin/env python
"""
Script to add or subtract two columns between one or two MS files
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import sys


def main(ms1, ms2, column1, column2, column_out, op='add'):
    """
    Add/subtract columns (column_out = column1 +/- column2)

    Note: we could also use TaQL to do this. E.g.:

        taql 'update ms1, ms2 t2 set column_out = column1 + t2.column2'

    but we have to create the output column first.

    Parameters
    ----------
    ms1 : str
        Name of MS file from which column 1 will be taken. This MS file will
        also receive the output column
    ms2 : str or list
        Name of MS file from which column 2 will be taken. Can be the same as
        ms1
    column1 : str
        Name of column 1
    column2 : str
        Name of column 2
    column_out : str
        Name of output column (written to ms1)
    op : str
        Operation to perform: 'add' or 'subtract'

    """
    if ms1 == ms2:
        ms2 = None

    # Read in the data
    t1 = pt.table(ms1, readonly=False, ack=False)
    data1 = t1.getcol(column1)
    if ms2 is not None:
        t2 = pt.table(ms2, readonly=False, ack=False)
        data2 = t2.getcol(column2)
    else:
        data2 = t1.getcol(column2)

    # Add the output column if needed
    if column_out not in t1.colnames():
        desc = t1.getcoldesc(column1)
        desc['name'] = column_out
        cd = pt.tableutil.makecoldesc(desc['name'], desc)
        tdesc = pt.tableutil.maketabdesc(cd)
        t1._addcols(tdesc, {}, True)
        t1._makerow()

    # Add or subtract columns
    if op.lower() == 'add':
        t1.putcol(column_out, data1 + data2)
    elif op.lower() == 'subtract':
        t1.putcol(column_out, data1 - data2)
    else:
        print('Operation not understood. Must be either "add" or "subtract"')
        sys.exit(1)
    t1.flush()
    t1.close()


if __name__ == '__main__':
    descriptiontext = "Add/subtract columns (column_out = column1 +/- column2).\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms1', help='name of MS file 1')
    parser.add_argument('ms2', help='name of MS file 2')
    parser.add_argument('column1', help='name of column 1')
    parser.add_argument('column2', help='name of column 2')
    parser.add_argument('column_out', help='name of the output column (written to ms1)')
    parser.add_argument('op', help='operation: "add" or "subtract"')
    args = parser.parse_args()

    main(args.ms1, args.ms2, args.column1, args.column2, args.column_out, args.op)
