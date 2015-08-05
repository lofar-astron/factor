#! /usr/bin/env python
"""
Script to copy beam info to WSClean MFS FITS images
"""
import argparse
from argparse import RawTextHelpFormatter
from astropy.io import fits as pyfits
import os


def main(image):
    """
    Copy the beam info to a WSClean MFS FITS image

    This script assumes standard WSClean naming conventions

    Parameters
    ----------
    image : str
        Filename of image in which the beam info will be entered

    """
    if os.path.isdir(image):
        # Image is a casa image, so skip it
        return

    templateim = image.replace('-MFS-', '-0000-')
    hduimtemplate = pyfits.open(templateim)
    hduim = pyfits.open(image, mode='update')
    headertemplate = hduimtemplate[0].header
    header = hduim[0].header
    bmaj = headertemplate['BMAJ']
    bmin = headertemplate['BMIN']
    bpa = headertemplate['BPA']

    header.update('BMAJ', bmaj, "")
    header.update('BMIN', bmin, "")
    header.update('BPA', bpa, "")

    hduim.flush()
    hduim.close()
    hduimtemplate.close()


if __name__ == '__main__':
    descriptiontext = "Copy beam info.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('image', help='filename of the WSClean MFS FITS image to copy to')
    args = parser.parse_args()

    main(args.image, args.templateim)
