#! /usr/bin/env python
"""
Script to copy an image for selfcal looping
"""
import argparse
from argparse import RawTextHelpFormatter
import sys


def main(image, counter):
    """
    Copy an image

    Parameters
    ----------
    image : str
        Name of image to copy
    counter : int
        Index of loop

    Returns
    -------
    result : dict
        Dict with image name from previous loop

    """
    counter = int(counter)

    image_copy = '{0}_{1}'.format(image, counter)

    if counter > 0:
        image_prev = '{0}_{1}'.format(image, counter-1)
    else:
        image_prev = image_copy

    return {'previous_image': image_prev}


if __name__ == '__main__':
    descriptiontext = "Copy an image.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('image', help='name of the image to copy')
    parser.add_argument('counter', help='loop index')
    args = parser.parse_args()

    main(args.image, args.counter)





