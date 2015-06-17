#! /usr/bin/env python
"""
Script to create a mosaic from facet images
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.tables
import pyrap.images as pim
from pyrap import quanta
import numpy as np
import pyfits
import os
from factor.directions import mask_vertices


def main(images, vertices, outfits, maxwidth=0):
    """
    Creates mosaic

    Parameters
    ----------
    images : str or list of str
        List of filenames of facet images. May be given as a list or as a string
        (e.g., '[image1, image2]'
    vertices : str or list of str
        List of filenames of facet vertices files. May be given as
        a list or as a string (e.g., '[vert1, vert2]'
    outfits : str
        Filename of output FITS mosaic
    maxwidth : int, optional
        Maximum number of pixels to consider for the width of the mosaic
        [default 0 = unlimited] This can be helpful at high declination.

    """
    if type(images) is str:
        images = images.strip('[]').split(',')
        images = [im.strip() for im in images]
    if type(vertices) is str:
        vertices = vertices.strip('[]').split(',')
        vertices = [v.strip() for v in vertices]

    formstr = '{0:45s}  {1:45s} {2:s}  {3:s} {4:s} {5:s}'
    print formstr.format("-----","--------","------------","-------","-------","------")
    print formstr.format("Image", "FC reg","Norm. weight", "Maj(ac)", "Min(ac)","PA(deg)")
    print formstr.format("-----","--------","------------","-------","-------","------")

    psf_fwhm = [] # resolution
    frequency = [] # frequency of images (should be equal?)
    for i in range(len(images)):
        this_pim = pim.image(images[i])
        info_dict = this_pim.info()['imageinfo']['restoringbeam']
        bpar_ma = quanta.quantity(info_dict['major']).get_value('deg')
        bpar_mi = quanta.quantity(info_dict['minor']).get_value('deg')
        bpar_pa = quanta.quantity(info_dict['positionangle']).get_value('deg')
        psf_fwhm.append([bpar_ma, bpar_mi, bpar_pa])
        frequency.append(this_pim.info()['coordinates']['spectral2']['restfreq'])
        print '{0:45.45s}  {1:45.45s} {2:0.2f}          {3:0.2f}    {4:0.2f}    {5:0.2f}'.format(images[i], vertices[i], 0, bpar_ma*60, bpar_mi*60,bpar_pa)

    psf_fwhm = np.array(psf_fwhm)
    frequency = np.array(frequency)
    mean_psf_fwhm = np.mean(psf_fwhm, axis=0)
    mean_frequency = np.mean(frequency)

    # Initialize some vectors
    declims = [] # store the limits of the declination axes
    raleft = []
    raright = []
    rainc = [] # store the r.a. increments in case they differ
    decinc = [] # store the dec increments in case they differ
    pims = [] # stores the pyrap images of the data

    # Get image frames for input images
    for im in images:
        image = pim.image(im)
        sptcoords = image.coordinates().get_coordinate('spectral')
        nc = sptcoords.get_axis_size()

        # Get Stokes axis. Ensure we are working with the Stokes parameter requested.
        stkcoords = image.coordinates().get_coordinate('stokes')
        if stkcoords.get_axis_size() == 1:
            assert(stkcoords.get_stokes()[0] == 'I')
        else:
            stks = stkcoords.get_stokes().index('I')
            image = image.subimage(blc=(0, stks), trc=(nc-1, stks), dropdegenerate=False)
        ns = 1

        dircoords = image.coordinates().get_coordinate('direction')
        nx = dircoords.get_axis_size(axis=1)
        ny = dircoords.get_axis_size(axis=0)
        inc = dircoords.get_increment()
        ref = dircoords.get_referencepixel()
        val = dircoords.get_referencevalue()
        # wsclean image header is weird
        if val[1]<0:
            val[1]+=2*np.pi
        ra_axis = (range(nx)-ref[1])*inc[1]+val[1]
        dec_axis = (range(ny)-ref[0])*inc[0]+val[0]
        rainc.append(inc[1])
        decinc.append(inc[0])
        declims.append(min(dec_axis))
        declims.append(max(dec_axis))
        mean_ra = np.mean(ra_axis)
        raleft.append((ra_axis[0]-mean_ra)*np.cos(val[0])+mean_ra)
        raright.append((ra_axis[-1]-mean_ra)*np.cos(val[0])+mean_ra)
        pims.append(image)

    # Generate the mosaic coordinate frame
    master_dec = np.arange(min(declims),max(declims),min(decinc))
    if max(raleft)-min(raright) > 5.*np.pi/3.: # crossed RA=0
        for i in range(len(raright)):
            raright[i] = raright[i]-2.*np.pi
    master_ra = np.arange(max(raleft),min(raright),max(rainc))
    lmra = len(master_ra)
    if args.maxwidth != 0:
        if lmra > args.maxwidth:
            xboundary = (lmra-args.maxwidth)/2
            master_ra = master_ra[xboundary:-xboundary]
    print "Found ra,dec pixel increments (arcsec):"
    print np.array(rainc)*206265.,np.array(decinc)*206265.
    ma = pims[-1].coordinates()
    ma['direction'].set_referencepixel([len(master_dec)/2,len(master_ra)/2])
    ma['direction'].set_increment([decinc[np.argmin(np.abs(decinc))],rainc[np.argmin(np.abs(rainc))]])
    ma['direction'].set_referencevalue([master_dec[len(master_dec)/2],master_ra[len(master_ra)/2]])

    # Initialize the arrays for the output image, sensitivity, and weights
    master_im = np.zeros((len(master_dec),len(master_ra)))

    # Reproject the images onto the master grid, weight and normalize
    for i, im in enumerate(pims):
        print 'doing image',i
        im, mask = mask_vertices(im, vertices[i])
        im = im.regrid([2,3],ma,outshape=(nc,ns,len(master_dec),len(master_ra)))
        mask = mask.regrid([2,3],ma,outshape=(nc,ns,len(master_dec),len(master_ra)))
        master_im += np.squeeze(im.getdata())
        master_mask += np.squeeze(mask.getdata())

    blank=np.ones_like(im)*np.nan
    master_im=np.where(master_mask,master_im,blank)

    # Write fits files
    arrax = np.zeros( (1,1, len(master_im[:,0]), len(master_im[0,:])) )
    arrax[0,0,:,:] = master_im

    # Open new casa image for mosaic
    new_pim = pim.image('',shape=(1,1, len(master_dec),len(master_ra)), coordsys=ma)
    new_pim.putdata(arrax)
    # Write fits
    new_pim.tofits(args.outfits, overwrite=True)

    # need to add new beam info (not sure if this is possible with pyrap)
    hdu = pyfits.open(args.outfits,mode='update')
    header = hdu[0].header
    header.update('BMAJ',mean_psf_fwhm[0])
    header.update('BMIN',mean_psf_fwhm[1])
    header.update('BPA',mean_psf_fwhm[2])
    header.update('BUNIT',pims[-1].info()['unit'])
    header.update('RESTFRQ',mean_frequency)
    header.update('RESTFREQ',mean_frequency)
    newhdu = pyfits.PrimaryHDU(data=hdu[0].data, header=header)
    newhdu.writeto(args.outfits,clobber=True)


if __name__ == '__main__':
    descriptiontext = "Create a mosaic from facet images.\n"
    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i','--images', help='List of filenames of facet images.')
    parser.add_argument('-v','--vertices', help='List of filenames of facet vertices.')
    parser.add_argument('-o','--outfits', help='Output name of mosaic fits file '
        '[default mosaic.fits]')
    parser.add_argument('-m','--maxwidth', help='Maximum number of pixels to '
        'consider for the width of the mosaic [default 0 = unlimited] This can '
        'be helpful at high declination.', default=0, type=int)

    args = parser.parse_args()
    main(args.images, args.vertices, args.outfits, args.maxwidth)
