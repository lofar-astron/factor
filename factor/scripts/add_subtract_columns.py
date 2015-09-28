#! /usr/bin/env python
"""
Script to add two columns between one or two MS files
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables as pt
import sys


def main(ms1, ms2, column1, column2, column_out, op='add'):
    """
    Add/subtract columns (column_out = column1 +/- column2)

    Parameters
    ----------
    ms1 : str
        Name of MS file to copy from
    ms2 : str or list
        Name of MS file or list of MS files to copy to. May be given as a list
        or as a string (e.g., '[ms1, ms2]'
    column1 : str
        Name of column 1
    column2 : str
        Name of column 2
    column_out : str
        Name of output column
    op : str
        Operation to perform: 'add' or 'subtract'

    """
    if ms1 == ms2:
        ms2 = None

    t1 = pt.table(ms1, readonly=False, ack=False)
    data1 = t1.getcol(column1)
    cd = t1.getcoldesc(column1)
    cd['name'] = column_out
    try:
        t.addcols(cd)
    except:
        pass

    if ms2 is not None:
        t2 = pt.table(ms2, readonly=False, ack=False)
        data2 = t2.getcol(column2)
    else:
        data2 = t1.getcol(column2)

    if op.lower() == 'add':
        data_out = data1 + data2
    elif op.lower() == 'subtract':
        data_out = data1 - data2
    else:
        print('Operation not understood. Must be either "add" or "subtract"')
        sys.exit(1)

    t.putcol(column_out, data_out)
    t.flush()
    t.close()


if __name__ == '__main__':
    descriptiontext = "Add/subtract columns (column_out = column1 +/- column2).\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms1', help='name of MS file 1')
    parser.add_argument('ms2', help='name of MS file 2')
    parser.add_argument('column1', help='name of column 1')
    parser.add_argument('column2', help='name of column 2')
    parser.add_argument('column_out', help='name of the output column')
    parser.add_argument('op', help='operation: "add" or "subtract"')
    args = parser.parse_args()

    main(args.ms1, args.ms2, args.column1, args.column2, args.column_out, args.op)





