#!/usr/bin/python
"""
Script to apply a primary-beam correction to a mosaic image
"""
import argparse
from argparse import RawTextHelpFormatter
from astropy.io import fits as pf
import astropy.wcs as pywcs
import os
import sys
import numpy as np
import scipy.ndimage
from scipy import interpolate


def main(mosaicfits, pbfits, outroot):
    """
    Corrects mosaic image with primary beam

    Parameters
    ----------
    mosaicfits : str
        Filename of mosaic image
    pbfits : str
        Filename of zeroed avgpb image
    outroot : str
        Filename root of output files (outroot.pbcor.fits and outroot.pbcut.fits)

    """
    pb_rescaled_fits = mosaicfits.replace('.fits', '') + '.pb.fits'

    if not os.path.exists(mosaicfits):
        raise Exception, "missing file: {m}".format(m=mosaicfits)
    if not os.path.exists(pbfits):
        raise Exception, "missing file: {m}".format(m=pbfits)

    mosaicpbfits = '{0}.pbcor.fits'.format(outroot)
    if os.path.isfile(mosaicpbfits):
        print "warning: overwriting {m}".format(m=mosaicpbfits)
        os.system("rm -rf %s" %(mosaicpbfits))
    mosaicpbcutfits = '{0}.pbcut.fits'.format(outroot)
    if os.path.isfile(mosaicpbcutfits):
        print "warning: overwriting {m}".format(m=mosaicpbcutfits)
        os.system("rm -rf %s" %(mosaicpbcutfits))

    mosaichead = pf.getheader(mosaicfits)
    mosaicdat = pf.getdata(mosaicfits)
    hdulistout = pf.open(mosaicfits)
    wcsout = pywcs.WCS(mosaichead)
    ra0out = mosaichead.get('OBSRA')
    dec0out = mosaichead.get('OBSDEC')

    S = mosaicdat.shape
    if len(S) == 4:
        print "dim 4"
        nc, nf, ny, nx = S
        dim = 4
    elif len(S) == 2:
        print "dim 2"
        ny, nx = S
        dim = 2
    else:
        raise Exception, "I don't know how to handle an image with this shape: "+str(S)

    pbcorout = mosaicdat.copy()
    mosaiccut = mosaicdat.copy()
    pbcorout = np.zeros_like(mosaicdat)

    pbhead = pf.getheader(pbfits)
    if pbhead['CDELT4'] == 0.0:
        # Causes WCS init problems if zero
        pbhead['CDELT4'] = -8.236827542606E+07
    pbcordat = pf.getdata(pbfits)
    pbhdulist = pf.open(pbfits)
    pbwcs = pywcs.WCS(pbhead)

    pbS = pbcordat.shape
    if len(pbS) == 4:
        print "dim 4"
        pbnc, pbnf, pbny, pbnx = pbS
        ndim = 4
    elif len(pbS) == 2:
        print "dim 2"
        pbny, pbnx = pbS
        ndim = 2
    else:
        raise Exception, "I don't know how to handle an image with this shape: "+str(pbS)

    rescale = False
    if pbnx != nx:
        rescale = True
    if pbny != ny:
        rescale = True

    rescale = True

    if rescale:
        # resample the pb image on the input image grid
        pbcordat_resampled = np.nan*np.ones_like(mosaicdat)

        if ndim == 4:
            print "resampling pb image... this may take a while"
            # loop by column
            prog = 0
            for xi in range(nx):
                if 100*xi/nx > prog:
                    prog = 100*xi/nx
                    if prog%10 == 0:
                        #print prog ,
                        sys.stdout.write(str(prog))
                        sys.stdout.flush()
                    else:
                        #print '.',
                        sys.stdout.write('.')
                        sys.stdout.flush()
                c = np.zeros(ny, dtype=int)
                f = np.zeros(ny, dtype=int)
                x = xi*np.ones(ny, dtype=int)
                y = np.arange(ny, dtype=int)
                pixcrd = np.array([x, y, c, f]).transpose()
                ra, dec, c, f = wcsout.wcs_pix2world(pixcrd, 0).transpose()
                worldcrd = np.array([ra, dec, c*0, f*0]).transpose()
                pbx, pby, pbc, pbf = pbwcs.wcs_world2pix(worldcrd, 0).transpose()
                pbx = np.array(pbx, dtype=int)
                pby = np.array(pby, dtype=int)

                for i in range(ny):
                    # outside the pb coverage
                    if (pbx[i] < 0) or (pby[i] < 0 ) or (pbx[i] >= pbnx ) or (pby[i] >= pbny):
                        pass
                    else:
                        pbcordat_resampled[0,0,y[i],x[i]] = pbcordat[0,0,pby[i],pbx[i]]

        elif ndim == 2:
            print "resampling pb image... this may take a while"
            # loop by column
            prog = 0
            for xi in range(nx):
                if 100*xi/nx > prog:
                    prog = 100*xi/nx
                    if prog%10 == 0:
                        sys.stdout.write(str(prog))
                        sys.stdout.flush()
                    else:
                        sys.stdout.write('.')
                        sys.stdout.flush()
                x = xi*np.ones(ny, dtype=int)
                y = np.arange(ny, dtype=int)
                pixcrd = np.array([x, y]).transpose()
                ra, dec = wcsout.wcs_pix2world(pixcrd, 0).transpose()
                worldcrd = np.array([ra, dec]).transpose()
                pbx, pby = pbwcs.wcs_world2pix(worldcrd, 0).transpose()
                pbx = np.array(pbx, dtype=int)
                pby = np.array(pby, dtype=int)

                for i in range(ny):
                    # outside the pb coverage
                    if (pbx[i] < 0) or (pby[i] < 0 ) or (pbx[i] >= pbnx ) or (pby[i] >= pbny):
                        pass
                    else:
                        pbcordat_resampled[y[i],x[i]] = pbcordat[pbx[i],pby[i]]

    Pcut = 0.4  # cut at Pcut power point of PB
    if rescale:
        pbcordat_resampled[pbcordat_resampled**0.5<Pcut] = np.nan  # set nan beyond
        mosaiccut[pbcordat_resampled**0.5<Pcut] = np.nan  # set nan beyond
        mosaiccut[np.isnan(pbcordat_resampled)] = np.nan  # set nan beyond
        mosaiccor = mosaicdat/(pbcordat_resampled**0.5)   #assuming awimager output is avgpb

        if os.path.isfile(pb_rescaled_fits):
            print "warning: overwriting {m}".format(m=pb_rescaled_fits)
            os.system("rm -rf %s" %(pb_rescaled_fits))
        pf.writeto(pb_rescaled_fits, pbcordat_resampled**0.5, header=mosaichead)
    else:
        pbcordat[pbcordat**0.5<Pcut] = np.nan
        mosaiccut[pbcordat**0.5<Pcut] = np.nan
        mosaiccut[np.isnan(pbcordat)] = np.nan

        mosaiccor = mosaicdat/(pbcordat**0.5)   #assuming awimager output is avgpb

    pf.writeto(mosaicpbfits, mosaiccor, header=mosaichead)
    pf.writeto(mosaicpbcutfits, mosaiccut, header=mosaichead)


if __name__ == '__main__':
    descriptiontext = "Apply a primary-beam correction to a mosaic image.\n"
    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('mosaicfits', help='filenames of mosaic image')
    parser.add_argument('pbfits', help='filenames of pbcor image')
    parser.add_argument('outroot', help='Output root name of corrected mosaic fits file')

    args = parser.parse_args()
    main(args.mosaicfits, args.pbfits, args.outroot)
