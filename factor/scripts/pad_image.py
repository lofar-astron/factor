#! /usr/bin/env python
"""
Script to pad WSClean model images
"""
import argparse
from argparse import RawTextHelpFormatter
from astropy.io import fits as pyfits
import numpy as np
import glob
import sys
import os


def main(root, scalefactor=1.5):
    model_images = glob.glob(root + '-model.fits') + glob.glob(root + '-0*-model.fits') + glob.glob(root + '-MFS-model.fits')

    scalefactor = float(scalefactor)

    for infile in model_images:
        hdu = pyfits.open(infile)
        imdata = hdu[0].data[0, 0]

        (xsize, ysize) = imdata.shape
        assert(xsize == ysize)
        print 'size is', xsize

        padsize = int(xsize * scalefactor)
        offset = (padsize - xsize) / 2
        print 'padding to', padsize
        print 'offset is', offset

        newdata=np.zeros((1, 1, padsize, padsize))

        newdata[0, 0, offset:offset+xsize, offset:offset+xsize] = imdata
        hdu[0].data = newdata
        hdu[0].header['CRPIX1'] += offset
        hdu[0].header['CRPIX2'] += offset
        hdu.writeto(infile, clobber=True)

    return {'padsize': '{0} {0}'.format(padsize)}


if __name__ == '__main__':
    descriptiontext = "Pad WSClean model images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('root', help='root name of input images')
    parser.add_argument('scalefactor', help='padding factor')
    args = parser.parse_args()

    main(args.root, args.scalefactor)
