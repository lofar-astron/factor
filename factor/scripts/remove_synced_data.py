#! /usr/bin/env python
"""
Script to delete synced data file and parent directory if empty
"""
import argparse
from argparse import RawTextHelpFormatter
import sys
import glob
import shutil
import os


def main(filename):
    """
    Delete synced data file and parent directory if empty

    Parameters
    ----------
    filename : str
        Filename of file to delete

    """
    # Delete file
    if os.path.exists(filename):
        if os.path.isdir(filename):
            shutil.rmtree(filename)
        else:
            os.remove(filename)

    # Next delete the parent directory if empty
    if os.path.exists(os.path.dirname(filename)):
        try:
            os.rmdir(os.path.dirname(filename))
        except OSError:
            pass


if __name__ == '__main__':
    descriptiontext = "Delete synced data file.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('filename', help='filename of file to delete')
    args = parser.parse_args()

    main(args.filename)

