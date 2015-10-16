#! /usr/bin/env python
"""
Script to delete images for selfcal looping
"""
import argparse
from argparse import RawTextHelpFormatter
import sys
import glob
import shutil
import os


def main(imageroot, counter):
    """
    Delete images with given root name

    Deletion is only done when counter > 0 (i.e., after the first selfcal pass).

    Parameters
    ----------
    imageroot : str
        Root name of images to delete
    counter : int
        Index of loop

    """
    counter = int(counter)

    if counter > 0:
        images1  = glob.glob(imageroot + '.*')
        for image in images1:
            if os.path.exists(image):
                shutil.rmtree(image)
        images2  = glob.glob(imageroot.replace('image41', 'image42') + '.*')
        for image in images2:
            if os.path.exists(image):
                shutil.rmtree(image)


if __name__ == '__main__':
    descriptiontext = "Delete images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('imageroot', help='root name of the images to delete')
    parser.add_argument('counter', help='loop index')
    args = parser.parse_args()

    main(args.image, args.counter)

