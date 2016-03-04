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


def main(imageroot, counter, indx):
    """
    Delete images with given root name

    The root name is assumed to contain the string "image?2." (e.g., "image42.").

    Deletion is only done when counter > 0 (i.e., after the first selfcal pass).

    Parameters
    ----------
    imageroot : str
        Root name of images to delete
    counter : int
        Index of loop
    indx : int
        Index of imaging step (e.g., 4)

    """
    counter = int(counter)
    indx = int(indx)

    # Find the index of the imaging step in imageroot. This may not necessarily
    # match the current index (given by indx) due to the need to run this
    # script with an existing mapfile from a previous imaging step. Therefore,
    # we find the root index in imageroot and replace it with the current one
    for i in range(5):
        if imageroot.find('image{}2.'.format(i)) >= 0:
            root_indx = i
            break
    imageroot.replace('image{}2.'.format(root_indx), 'image{}2.'.format(indx))

    if counter > 0:
        # First delete the "image?2" images
        images1  = glob.glob(imageroot + '.*')
        for image in images1:
            if os.path.exists(image):
                if os.path.isdir(image):
                    shutil.rmtree(image)

        # Next delete the "image?1" images
        images2  = glob.glob(imageroot.replace('image{}2.'.format(indx), 'image{}1.'.format(indx)) + '.*')
        for image in images2:
            if os.path.exists(image):
                if os.path.isdir(image):
                    shutil.rmtree(image)


if __name__ == '__main__':
    descriptiontext = "Delete images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('imageroot', help='root name of the images to delete')
    parser.add_argument('counter', help='loop index')
    args = parser.parse_args()

    main(args.image, args.counter)

