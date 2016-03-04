#! /usr/bin/env python
"""
Script to copy an image for selfcal looping
"""
import argparse
from argparse import RawTextHelpFormatter
import sys
import shutil
import os


def main(image, counter, indx):
    """
    Copy an image

    Parameters
    ----------
    image : str
        Name of image to copy
    counter : int
        Index of loop
    indx : int
        Index of imaging step (e.g., 4)

    Returns
    -------
    result : dict
        Dict with image name from previous loop

    """
    counter = int(counter)
    indx = int(indx)

    image_copy = image.replace('image{0}2'.format(indx), 'image{0}2_iter{1}'.format(indx, counter))
    if os.path.exists(image_copy):
        shutil.rmtree(image_copy)
    shutil.copytree(image, image_copy)

    if counter > 0:
        # Use image from previous iteration of the current imaging step
        image_prev = image.replace('image{0}2'.format(indx), 'image{0}2_iter{1}'.format(indx, counter-1))
    else:
        # Use image from previous imaging step
        image_prev = image.replace('image{0}2'.format(indx), 'image{0}2'.format(indx-1))

    return {'previous_image': image_prev}


if __name__ == '__main__':
    descriptiontext = "Copy an image.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('image', help='name of the image to copy')
    parser.add_argument('counter', help='loop index')
    args = parser.parse_args()

    main(args.image, args.counter)





