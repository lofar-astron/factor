#!/usr/bin/env python
"""
Script to merge two sky models
"""
import argparse
from argparse import RawTextHelpFormatter
import numpy
import pyrap.images as pim
import os


def main(image_pre, image_post, res_val):
    """
    Check quantities in residual images
    """
    imgpre = pim.image(image_pre)
    pixelspre = numpy.copy(imgpre.getdata())
    maxvalpre = numpy.copy(numpy.max(pixelspre))

    img = pim.image(image_post)
    pixels = numpy.copy(img.getdata())
    maxval = numpy.copy(numpy.max(pixels))

    if (maxval > res_val): # or ((maxval*0.95) > maxvalpre) :
        logging.info('WARNING RESIDUAL TOO LARGE')
        logging.info('Max = {0}'.format(maxval))
        logging.info('Previous max = {0}'.format(maxvalpre))
        return {'break': False, 'maxval': maxval, 'maxvalpre': maxvalpre}
    else:
        return {'break': True, 'maxval': maxval, 'maxvalpre': maxvalpre}


if __name__ == '__main__':
    descriptiontext = "Check residual images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('image_pre', help='image before selfcal')
    parser.add_argument('image_post', help='image after selfcal')
    parser.add_argument('res_val', help='maximum acceptible residual in Jy')

    args = parser.parse_args()
    main(args.image_pre, args.image_post, args.res_val)
