#! /usr/bin/env python
"""
Script to make a clean mask from an image
"""
import argparse
from argparse import RawTextHelpFormatter
from lofar import bdsm
import numpy
import sys
import os


def main(image_name, mask_name, atrous_do=False, threshisl=0.0, threshpix=0.0, rmsbox=None,
         iterate_threshold=False, adaptive_rmsbox=False, img_format='fits',
         threshold_format='float'):
    """
    Run PyBDSM to make an island clean mask

    Parameters
    ----------
    TODO

    Returns
    -------
    result : dict
        Dict with 5-sigma threshold

    """
    if atrous_do:
        threshisl = 4.0

    if rmsbox is not None:
        rmsbox = eval(rmsbox)

    if type(atrous_do) is str:
        if atrous_do.lower() == 'true':
            atrous_do = True
        else:
            atrous_do = False

    if type(iterate_threshold) is str:
        if iterate_threshold.lower() == 'true':
            iterate_threshold = True
        else:
            iterate_threshold = False

    if type(adaptive_rmsbox) is str:
        if adaptive_rmsbox.lower() == 'true':
            adaptive_rmsbox = True
        else:
            adaptive_rmsbox = False

    if iterate_threshold:
        # Start with high threshold and lower it until we get at least one island
        threshpix_orig = threshpix
        threshisl_orig = threshisl
        nisl = 0
        threshpix = 25
        threshisl = 15
        while nisl == 0:
            img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                                     thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
                                     atrous_do=atrous_do, ini_method='curvature',
                                     adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=20, quiet=True)
            nisl = img.nisl
            threshpix /= 1.2
            threshisl /= 1.2
        threshpix = threshpix_orig
        threshisl = threshisl_orig
    else:
        img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                                 thresh_pix=numpy.float(threshpix), thresh_isl=numpy.float(threshisl),
                                 atrous_do=atrous_do, ini_method='curvature',
                                 adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=20, quiet=True)

    img.export_image(img_type='island_mask', mask_dilation=0, outfile=mask_name,
                     img_format=img_format, clobber=True)

    if threshold_format == 'float':
        return {'threshold_5sig': 5.0 * img.clipped_rms}
    elif threshold_format == 'str_with_units':
        return {'threshold_5sig': "'{0}Jy'".format(5.0 * img.clipped_rms)}


if __name__ == '__main__':
    descriptiontext = "Make a clean mask.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('image_name', help='Image name')
    parser.add_argument('mask_name', help='Mask name')
    parser.add_argument('-a', '--atrous_do', help='use wavelet fitting', type=bool, default=False)
    parser.add_argument('-i', '--threshisl', help='', type=float, default=3.0)
    parser.add_argument('-p', '--threshpix', help='', type=float, default=5.0)
    parser.add_argument('-r', '--rmsbox', help='rms box width and step (e.g., "(60, 20)")',
        type=str, default='(60, 20)')
    parser.add_argument('-t', '--iterate_threshold', help='iteratively decrease threshold until at least '
        'one island is found', type=bool, default=False)
    parser.add_argument('-o', '--adaptive_rmsbox', help='use an adaptive rms box', type=bool, default=False)
    parser.add_argument('-f', '--img_format', help='format of output mask', type=str, default='casa')
    parser.add_argument('-d', '--threshold_format', help='format of return value', type=str, default='float')

    args = parser.parse_args()
    main(args.image_name, args.mask_name, atrous_do=args.atrous_do,
         threshisl=args.threshisl, threshpix=args.threshpix, rmsbox=args.rmsbox, iterate_threshold=args.iterate_threshold,
         adaptive_rmsbox=args.adaptive_rmsbox, img_format=args.img_format,
         threshold_format=args.threshold_format)
