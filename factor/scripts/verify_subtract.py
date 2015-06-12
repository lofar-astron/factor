#!/usr/bin/env python
"""
Script to verify selfcal subtract
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
    maxvalpre = numpy.max(imgpre.getdata())
    imgpost = pim.image(image_post)
    maxvalpost = numpy.max(imgpost.getdata())

    if (maxvalpost > res_val): # or ((maxvalpost*0.95) > maxvalpre) :
        return {'break': False, 'maxval': maxvalpost, 'maxvalpre': maxvalpre}
    else:
        return {'break': True, 'maxval': maxvalpost, 'maxvalpre': maxvalpre}


if __name__ == '__main__':
    descriptiontext = "Check residual images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('image_pre', help='image before selfcal')
    parser.add_argument('image_post', help='image after selfcal')
    parser.add_argument('res_val', help='maximum acceptible residual in Jy')

    args = parser.parse_args()
    main(args.image_pre, args.image_post, args.res_val)
