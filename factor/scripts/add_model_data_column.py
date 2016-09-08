#! /usr/bin/env python
"""
Script to add MODEL_DATA column to MS
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables as pt
import sys
import os


def main(ms):
    """
    Add MODEL_DATA column

    Parameters
    ----------
    ms : str
        Name of input MS file

    """
    # Add the MODEL_DATA column if needed
    t1 = pt.table(ms, readonly=False, ack=False)
    if 'MODEL_DATA' not in t1.colnames():
        desc = t1.getcoldesc('DATA')
        desc['name'] = 'MODEL_DATA'
        t1.addcols(desc)
    t1.close()


if __name__ == '__main__':
    descriptiontext = "Add MODEL_DATA column.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('ms', help='name of MS file')
    args = parser.parse_args()

    main(args.ms_list)
