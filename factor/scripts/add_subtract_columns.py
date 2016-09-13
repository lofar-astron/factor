#! /usr/bin/env python
"""
Script to add or subtract two columns between one or two MS files
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables as pt
import sys
import os


def main(ms1, ms2, column1, column2, column_out, op='add', in_memory=False):
    """
    Add/subtract columns (column_out = column1 +/- column2)

    Parameters
    ----------
    ms1 : str
        Name of MS file from which column 1 will be taken. This MS file will
        also receive the output column
    ms2 : str or list
        Name of MS file from which column 2 will be taken.
    column1 : str
        Name of column 1
    column2 : str
        Name of column 2
    column_out : str
        Name of output column (written to ms1)
    op : str, optional
        Operation to perform: 'add' or 'subtract'
    in_memory : bool, optional
        If True, do the subtraction in memory rather than with taql

    """
    if type(in_memory) is str:
        if in_memory.lower() == 'true':
            in_memory = True
        else:
            in_memory = False

    # Add the output column to ms1 if needed
    t1 = pt.table(ms1, readonly=False, ack=False)
    if column_out not in t1.colnames():
        desc = t1.getcoldesc(column1)
        desc['name'] = column_out
        t1.addcols(desc)
    t1.close()

    if in_memory:
        # Add or subtract columns in memory
        t1 = pt.table(ms1, readonly=False, ack=False)
        data1 = t1.getcol(column1)
        t2 = pt.table(ms2, ack=False)
        data2 = t2.getcol(column2)
        t2.close()
        if op.lower() == 'add':
            t1.putcol(column_out, data1 + data2)
        elif op.lower() == 'subtract':
            t1.putcol(column_out, data1 - data2)
        else:
            print('Operation not understood. Must be either "add" or "subtract"')
            sys.exit(1)
        t1.flush()
        t1.close()
    else:
        # Add or subtract columns with TaQL
        if op.lower() == 'add':
            op_sym = '+'
        elif op.lower() == 'subtract':
            op_sym = '-'
        else:
            print('Operation not understood. Must be either "add" or "subtract"')
            sys.exit(1)
        os.system("taql 'update {0}, {1} t2 set {2}={3}{4}t2.{5}'".format(
            ms1, ms2, column_out, column1, op_sym, column2))


if __name__ == '__main__':
    descriptiontext = "Add/subtract columns (column_out = column1 +/- column2).\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms1', help='name of MS file 1')
    parser.add_argument('ms2', help='name of MS file 2')
    parser.add_argument('column1', help='name of column 1')
    parser.add_argument('column2', help='name of column 2')
    parser.add_argument('column_out', help='name of the output column (written to ms1)')
    parser.add_argument('op', help='operation: "add" or "subtract"')
    parser.add_argument('in_memory', help='do operation in memory')
    args = parser.parse_args()

    main(args.ms1, args.ms2, args.column1, args.column2, args.column_out, args.op, args.in_memory)
