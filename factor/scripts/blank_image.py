#! /usr/bin/env python
"""
Script to blank regions (with zeros or NaNs) in a fits image
"""
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import sys
import os
import pickle
import glob
from factor.lib.polygon import Polygon
from astropy.io import fits as pyfits
from astropy import wcs


def read_vertices(filename):
    """
    Returns facet vertices stored in input file
    """
    with open(filename, 'r') as f:
        direction_dict = pickle.load(f)
    return direction_dict['vertices']


def main(input_image_file, vertices_file, output_image_file, blank_value='zero',
    image_is_wsclean_model=False):
    """
    Blank a region in an image

    Parameters
    ----------
    input_image_file : str
        Filename of input image to blank
    vertices_file : str, optional
        Filename of file with vertices (must be a pickle file containing
        a dictionary with the vertices in the 'vertices' entry)
    output_image_file : str
        Filename of output image
    blank_value : str, optional
        Value for blanks (one of 'zero' or 'nan')
    image_is_wsclean_model : bool, optional
        If True, the input and output image files are treated as the root name
        of a WSClean model image (or images)

    """
    if type(image_is_wsclean_model) is str:
        if image_is_wsclean_model.lower() == 'true':
            image_is_wsclean_model = True
        else:
            image_is_wsclean_model = False

    if image_is_wsclean_model:
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

    # Construct polygon of facet region
    header = pyfits.getheader(input_image_files[0], 0)
    w = wcs.WCS(header)
    RAind = w.axis_type_names.index('RA')
    Decind = w.axis_type_names.index('DEC')
    vertices = read_vertices(vertices_file)
    RAverts = vertices[0]
    Decverts = vertices[1]
    xvert = []
    yvert = []
    for RAvert, Decvert in zip(RAverts, Decverts):
        ra_dec = np.array([[0.0, 0.0, 0.0, 0.0]])
        ra_dec[0][RAind] = RAvert
        ra_dec[0][Decind] = Decvert
        xvert.append(w.wcs_world2pix(ra_dec, 0)[0][Decind])
        yvert.append(w.wcs_world2pix(ra_dec, 0)[0][RAind])
    poly = Polygon(xvert, yvert)

    for input_image, output_image in zip(input_image_files, output_image_files):
        hdu = pyfits.open(input_image)
        data = hdu[0].data

        # Find limits of facet poly and blank pixels outside them
        xmin = max(int(np.min(xvert)) - 2, 0)
        xmax = min(int(np.max(xvert)) + 2, data.shape[2])
        ymin = max(int(np.min(yvert)) - 2, 0)
        ymax = min(int(np.max(yvert)) + 2, data.shape[3])
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

        hdu[0].data = data
        hdu.writeto(output_image, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Blank regions of an image.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('input_image_file', help='Filename of input image')
    parser.add_argument('vertices_file', help='Filename of vertices file')
    parser.add_argument('output_image_file', help='Filename of output image')
    parser.add_argument('-b', '--blank_value', help='value for blank pixesl', type=str, default='zeros')
    args = parser.parse_args()
    main(args.input_image_file, args.vertices_file, args.output_image_file,
         blank_value=args.blank_value)
