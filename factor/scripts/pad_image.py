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


def get_optimum_size(size):
    """
    Gets the nearest optimum image size

    Taken from the casa source code (cleanhelper.py)

    Parameters
    ----------
    size : int
        Target image size in pixels

    Returns
    -------
    optimum_size : int
        Optimum image size nearest to target size

    """
    def prime_factors(n, douniq=True):
        """ Return the prime factors of the given number. """
        factors = []
        lastresult = n
        sqlast=int(np.sqrt(n))+1
        if n == 1:
            return [1]
        c=2
        while 1:
             if (lastresult == 1) or (c > sqlast):
                 break
             sqlast=int(np.sqrt(lastresult))+1
             while 1:
                 if(c > sqlast):
                     c=lastresult
                     break
                 if lastresult % c == 0:
                     break
                 c += 1

             factors.append(c)
             lastresult /= c

        if (factors==[]): factors=[n]
        return np.unique(factors).tolist() if douniq else factors

    n = int(size)
    if (n%2 != 0):
        n+=1
    fac=prime_factors(n, False)
    for k in range(len(fac)):
        if (fac[k] > 7):
            val=fac[k]
            while (np.max(prime_factors(val)) > 7):
                val +=1
            fac[k]=val
    newlarge=np.product(fac)
    for k in range(n, newlarge, 2):
        if ((np.max(prime_factors(k)) < 8)):
            return k
    return newlarge


def main(root, scalefactor=1.5):
    model_images = glob.glob(root + '-model.fits') + glob.glob(root + '-0*-model.fits') + glob.glob(root + '-MFS-model.fits')

    scalefactor = float(scalefactor)

    for infile in model_images:
        hdu = pyfits.open(infile)
        imdata = hdu[0].data[0, 0]

        (xsize, ysize) = imdata.shape
        assert(xsize == ysize)
        print 'size is', xsize

        padsize = get_optimum_size(int(xsize * scalefactor))
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
