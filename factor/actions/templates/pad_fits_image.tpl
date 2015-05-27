#!/usr/bin/python

# pad a FITS file

import pyfits
import sys
import os
import numpy as np

def padfits(infile, outfile, scalefactor=1.2):

    hdu = pyfits.open(infile)
    imdata = hdu[0].data[0,0]

    (xsize,ysize) = imdata.shape
    assert(xsize == ysize)

    padsize = int(xsize * 1.2)
    offset = (padsize - xsize) / 2

    newdata = np.zeros((1, 1, padsize, padsize))

    newdata[0,0, offset:offset+xsize, offset:offset+xsize] = imdata
    hdu[0].data = newdata
    hdu[0].header['CRPIX1'] += offset
    hdu[0].header['CRPIX2'] += offset
    hdu.writeto(outfile, clobber=True)
    return padsize

if __name__=='__main__':
    filename = sys.argv[1]
    split = filename.split('.fits')
    outfile = split[0] + '_padded.fits'
    padfits(filename, outfile)
