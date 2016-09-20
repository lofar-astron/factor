#! /usr/bin/env python
"""
Script to sync files
"""
import argparse
from argparse import RawTextHelpFormatter
import os


def main(file_from, file_to):
    """
    Sync a file

    Parameters
    ----------
    file_from : str
        Name of file to from
    file_to : str
        Name of file to copy to

    """
    destination_dir = os.path.dirname(file_to)
    if not os.path.exists(file_to):
        # Copy file
        os.system('/bin/cp -rT {0} {1}'.format(file_from, file_to))
    else:
        # Sync file
        os.system('/usr/bin/rsync -a {0} {1}'.format(file_from, destination_dir))


if __name__ == '__main__':
    descriptiontext = "Sync a file.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('file_from', help='name of the image to copy')
    parser.add_argument('file_to', help='loop counter')
    args = parser.parse_args()

    main(args.file_from, args.file_to)





