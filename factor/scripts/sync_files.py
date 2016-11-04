#! /usr/bin/env python
"""
Script to sync files
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.tables as pt
import numpy as np
import sys
import os


def main(file_from, file_to):
    """
    Sync a file

    Parameters
    ----------
    file_from : str
        Name of file to copy from. Can be list of files such as '[file1,file2]'
    file_to : str
        Name of file to copy to

    """
    if file_from.startswith('[') and file_from.endswith(']'):
        # Assume both inputs are lists
        file_from = file_from.strip('[]').split(',')
        files_from = [f.strip() for f in file_from]
        file_to = file_to.strip('[]').split(',')
        files_to = [f.strip() for f in file_to]
    else:
        files_from = [file_from]
        files_to = [file_to]

    for file_from, file_to in zip(files_from, files_to):
        if file_from.endswith('/'):
            file_from = file_from.strip('/')
        if file_to.endswith('/'):
            file_to = file_to.strip('/')

        destination_dir = os.path.dirname(file_to)
        if len(destination_dir) > 0:
            os.system('/bin/mkdir -p {0}'.format(destination_dir))

        # Copy files
        try:
            os.system('/bin/cp -r {0} {1}'.format(file_from, file_to))
        except:
            print('ERROR: could not sync file to {}, possibly due to lack of space.'.format(destination_dir))
            if os.path.exists(file_to):
                os.system('rm -rf {}'.format(file_to))
            sys.exit(1)


if __name__ == '__main__':
    descriptiontext = "Sync a file.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('file_from', help='name of the image to copy')
    parser.add_argument('file_to', help='loop counter')
    args = parser.parse_args()

    main(args.file_from, args.file_to)
