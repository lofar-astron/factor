#! /usr/bin/env python
"""
Script to make a clean mask from an image
"""
import argparse
from argparse import RawTextHelpFormatter
from lofar import bdsm
import pyrap.images as pim
import pickle
import numpy as np
import sys
import os
from factor.directions import Polygon


def read_vertices(filename):
    """
    Returns facet vertices
    """
    with open(filename, 'r') as f:
        direction_dict = pickle.load(f)
    return direction_dict['vertices']


def main(image_name, mask_name, atrous_do=False, threshisl=0.0, threshpix=0.0, rmsbox=None,
         iterate_threshold=False, adaptive_rmsbox=False, img_format='fits',
         threshold_format='float', trim_by=25, vertices_file=None, atrous_jmax=6):
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

    trim_by = int(trim_by)
    atrous_jmax = int(atrous_jmax)
    threshpix = float(threshpix)
    threshisl = float(threshisl)

    if iterate_threshold:
        # Start with high threshold and lower it until we get at least one island
        nisl = 0
        while nisl == 0:
            img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                                     thresh_pix=threshpix, thresh_isl=threshisl,
                                     atrous_do=atrous_do, ini_method='curvature', thresh='hard',
                                     adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=150,
                                     rms_box_bright=(35,7), rms_map=True, quiet=True,
                                     atrous_jmax=atrous_jmax)
            nisl = img.nisl
            threshpix /= 1.2
            threshisl /= 1.2
            if threshpix < 5.0:
                break
    else:
        img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox,
                                 thresh_pix=threshpix, thresh_isl=threshisl,
                                 atrous_do=atrous_do, ini_method='curvature', thresh='hard',
                                 adaptive_rms_box=adaptive_rmsbox, adaptive_thresh=150,
                                 rms_box_bright=(35,7), rms_map=True, quiet=True,
                                 atrous_jmax=atrous_jmax)

    if img.nisl == 0:
        print('No islands found. Clean mask cannot be made.')
        sys.exit(1)

    img.export_image(img_type='island_mask', mask_dilation=0, outfile=mask_name,
                     img_format=img_format, clobber=True)

    if vertices_file is not None or trim_by > 0:
        # Modify the clean mask to exclude regions outside of the polygon and
        # trim edges
        mask_tmp_name = mask_name + '.tmp1'
        if os.path.exists(mask_tmp_name):
            os.system('rm -rf {0}'.format(mask_tmp_name))
        os.system('cp -r {0} {1}'.format(mask_name, mask_tmp_name))

        mask_im = pim.image(mask_tmp_name)
        img_type = mask_im.imagetype()
        if img_type == 'FITSImage':
            mask_im.saveas(mask_name+'.tmp2')
            mask_im = pim.image(mask_name+'.tmp2')
        data = mask_im.getdata()

        if vertices_file is not None:
            vertices = read_vertices(vertices_file)
            RAverts = vertices[0]
            Decverts = vertices[1]
            xvert = []
            yvert = []
            for RAvert, Decvert in zip(RAverts, Decverts):
                pixels = mask_im.topixel([0, 1, Decvert*np.pi/180.0,
                    RAvert*np.pi/180.0])
                xvert.append(pixels[2]) # x -> Dec
                yvert.append(pixels[3]) # y -> RA
            poly = Polygon(xvert, yvert)

            # Find masked regions
            masked_ind = np.where(data[0, 0])

            # Find distance to nearest poly edge and unmask those that
            # are outside the facet (dist < 0)
            dist = poly.is_inside(masked_ind[0], masked_ind[1])
            outside_ind = np.where(dist < 0.0)
            if len(outside_ind[0]) > 0:
                data[0, 0, masked_ind[0][outside_ind], [masked_ind[1][outside_ind]]] = 0

        if trim_by > 0:
            sh = np.shape(data)
            data[0, 0, 0:sh[2], 0:trim_by] = 0
            data[0, 0, 0:trim_by, 0:sh[3]] = 0
            data[0, 0, 0:sh[2], sh[3]-trim_by:sh[3]] = 0
            data[0, 0, sh[2]-trim_by:sh[2], 0:sh[3]] = 0

        # Save changes
        mask_im.putdata(data)
        if img_format == 'fits':
            mask_im.tofits(mask_name, overwrite=True)
        else:
            mask_im.saveas(mask_name, overwrite=True)

    if threshold_format == 'float':
        return {'threshold_5sig': 5.0 * img.clipped_rms}
    elif threshold_format == 'str_with_units':
        # This is done to get around the need for quotes around strings in casapy scripts
        # 'casastr/' is removed by the generic pipeline
        return {'threshold_5sig': 'casastr/{0}Jy'.format(5.0 * img.clipped_rms)}


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
    parser.add_argument('-b', '--trim_by', help='Trim masked region by this number of pixels', type=int, default=25)
    parser.add_argument('-v', '--vertices_file', help='file containing facet polygon vertices', type=str, default=None)
    parser.add_argument('-j', '--atrous_jmax', help='Max wavelet scale', type=int, default=3)

    args = parser.parse_args()
    main(args.image_name, args.mask_name, atrous_do=args.atrous_do,
         threshisl=args.threshisl, threshpix=args.threshpix, rmsbox=args.rmsbox,
         iterate_threshold=args.iterate_threshold,
         adaptive_rmsbox=args.adaptive_rmsbox, img_format=args.img_format,
         threshold_format=args.threshold_format, trim_by=args.trim_by,
         vertices_file=args.vertices_file, atrous_jmax=args.atrous_jmax)
