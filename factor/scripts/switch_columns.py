#! /usr/bin/env python
"""
Script to switch the names of normal and baseline-averaged columns
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables as pt
import numpy
import sys


def main(ms_file):
    """
    Switch the names of normal and baseline-averaged columns

    Parameters
    ----------
    ms_file : str
        Name of MS file

    """
    t = pt.table(ms_file, readonly=False, ack=False)

    for (column1, column2) in zip(['DATA', 'WEIGHT_SPECTRUM'], ['BLAVG_DATA', 'BLAVG_WEIGHT_SPECTRUM']):
        if column1 not in t.colnames() and column2 not in t.colnames():
            print('Both columns must be present in MS')
            sys.exit(1)

        # Rename column1 to temp col
        t.renamecol(column1, column1+'_TEMP')

        # Rename column2 to column1
        t.renamecol(column2, column1)

        # Rename temp col to column2
        t.renamecol(column1+'_TEMP', column2)
    t.flush()
    t.close()


if __name__ == '__main__':
    descriptiontext = "Switch the names of two columns.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms_file', help='name of the MS file')
    parser.add_argument('column1', help='name of column 1')
    parser.add_argument('column2', help='name of column 2')
    args = parser.parse_args()

    main(args.ms_file, args.column1, args.column2)
