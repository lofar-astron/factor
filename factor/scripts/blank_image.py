#! /usr/bin/env python
"""
Script to blank regions (with zeros or NaNs) in a fits or casa image
"""
import argparse
from argparse import RawTextHelpFormatter
import casacore.images as pim
from casacore import quanta
import numpy as np
import sys
import os
import pickle
import glob
from factor.lib.polygon import Polygon
from astropy.io import fits as pyfits


def read_vertices(filename):
    """
    Returns facet vertices stored in input file
    """
    with open(filename, 'r') as f:
        direction_dict = pickle.load(f)
    return direction_dict['vertices']


def main(input_image_file, vertices_file, output_image_file, blank_value='zero',
    img_format='fits', image_is_casa_model=False, image_is_wsclean_model=False,
    nterms=1):
    """
    Blank a region in an image

    Parameters
    ----------
    input_image_file : str
        Filename of input image from which mask will be made. If the image does
        not exist, a template image with center at (reference_ra_deg,
        reference_dec_deg) will be made internally
    vertices_file : str, optional
        Filename of file with vertices (must be a pickle file containing
        a dictionary with the vertices in the 'vertices' entry)
    output_image_file : str
        Filename of output image
    blank_value : str, optional
        Value for blanks (one of 'zero' or 'nan')
    img_format : str, optional
        Format of output mask image (one of 'fits' or 'casa')
    image_is_casa_model : bool, optional
        If True, the input and output image files are treated as the root name
        of a casa model image (or images)
    image_is_wsclean_model : bool, optional
        If True, the input and output image files are treated as the root name
        of a WSClean model image (or images)
    nterms : int, optional
        If image_is_casa_model is True, this argument sets the number of nterms
        for the model

    """
    if type(image_is_casa_model) is str:
        if image_is_casa_model.lower() == 'true':
            image_is_casa_model = True
        else:
            image_is_casa_model = False
    if type(image_is_wsclean_model) is str:
        if image_is_wsclean_model.lower() == 'true':
            image_is_wsclean_model = True
        else:
            image_is_wsclean_model = False

    if type(nterms) is str:
        nterms = int(nterms)

    if image_is_casa_model:
        if nterms == 1:
            input_image_files = [input_image_file+'.model']
            output_image_files = [output_image_file+'.model']
        if nterms == 2:
            input_image_files = [input_image_file+'.model.tt0',
                input_image_file+'.model.tt1']
            output_image_files = [output_image_file+'.model.tt0',
                output_image_file+'.model.tt1']
        if nterms == 3:
            input_image_files = [input_image_file+'.model.tt0',
                input_image_file+'.model.tt1', input_image_file+'.model.tt2']
            output_image_files = [output_image_file+'.model.tt0',
                output_image_file+'.model.tt1', output_image_file+'.model.tt2']
    elif image_is_wsclean_model:
        input_image_files = glob.glob(input_image_file+'*-model.fits')
        output_image_files = [f.replace(input_image_file, output_image_file) for f in input_image_files]
    else:
        input_image_files = [input_image_file]
        output_image_files = [output_image_file]

    if blank_value == 'zero':
        blank_val = 0.0
    elif blank_value == 'nan':
        blank_val = np.nan
    else:
        print('Blank value type "{}" not understood.'.format(blank_with))
        sys.exit(1)

    for input_image, output_image in zip(input_image_files, output_image_files):
        im = pim.image(input_image)
        data = im.getdata()
        coordsys = im.coordinates()
        imshape = im.shape()
        new_im = pim.image('', shape=imshape, coordsys=coordsys)

        # Construct polygon
        vertices = read_vertices(vertices_file)
        RAverts = vertices[0]
        Decverts = vertices[1]
        xvert = []
        yvert = []
        for RAvert, Decvert in zip(RAverts, Decverts):
            pixels = new_im.topixel([1, 1, Decvert*np.pi/180.0,
                RAvert*np.pi/180.0])
            xvert.append(pixels[2]) # x -> Dec
            yvert.append(pixels[3]) # y -> RA
        poly = Polygon(xvert, yvert)

        # Find size of facet and blank regions outside
        xmin = int(np.min(xvert)) - 2
        xmax = int(np.max(xvert)) + 2
        ymin = int(np.min(yvert)) - 2
        ymax = int(np.max(yvert)) + 2
        data[0, 0, :, :ymin] = blank_val
        data[0, 0, :, ymax:] = blank_val
        data[0, 0, :xmin, :] = blank_val
        data[0, 0, xmax:, :] = blank_val

        # Find distance to nearest poly edge and blank those that
        # are outside the facet (dist < 0)
        pix_ind = np.indices((xmax-xmin, ymax-ymin))
        pix_ind[0] += xmin
        pix_ind[1] += ymin
        dist = poly.is_inside(pix_ind[0], pix_ind[1])
        outside_ind = np.where(dist < 0.0)
        if len(outside_ind[0]) > 0:
            data[0, 0, pix_ind[0][outside_ind], pix_ind[1][outside_ind]] = blank_val

        # Save changes
        new_im.putdata(data)
        if img_format == 'fits':
            new_im.tofits(output_image, overwrite=True)
            info_dict = im.info()['imageinfo']['restoringbeam']
            hdu = pyfits.open(output_image, mode='update')
            header = hdu[0].header
            header['BMAJ'] = quanta.quantity(info_dict['major']).get_value('deg')
            header['BMIN'] = quanta.quantity(info_dict['minor']).get_value('deg')
            header['BPA'] = quanta.quantity(info_dict['positionangle']).get_value('deg')
            header['BUNIT'] = im.info()['unit']
            header['RESTFRQ'] = im.info()['coordinates']['spectral2']['restfreq']
            header['RESTFREQ'] = im.info()['coordinates']['spectral2']['restfreq']
            newhdu = pyfits.PrimaryHDU(data=hdu[0].data, header=header)
            newhdu.writeto(output_image,clobber=True)

        elif img_format == 'casa':
            new_im.saveas(output_image, overwrite=True)
        else:
            print('Output image format "{}" not understood.'.format(img_format))
            sys.exit(1)


if __name__ == '__main__':
    descriptiontext = "Blank regions of an image.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('input_image_file', help='Filename of input image')
    parser.add_argument('vertices_file', help='Filename of vertices file')
    parser.add_argument('output_image_file', help='Filename of output image')
    parser.add_argument('-b', '--blank_value', help='value for blank pixesl', type=str, default='zeros')
    parser.add_argument('-f', '--img_format', help='format of output mask', type=str, default='casa')
    args = parser.parse_args()
    main(args.input_image_file, args.vertices_file, args.output_image_file,
         blank_value=args.blank_value, img_format=args.img_format)
