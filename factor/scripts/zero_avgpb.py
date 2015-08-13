#!/usr/bin/env python
"""
Script to zero the corners of avgpb images
"""
import argparse
from argparse import RawTextHelpFormatter
import pyrap.images as pim
import numpy as np


def main(image, radius=0.5, output=None):
    """
    Zero corners of avgpb images

    Parameters
    ----------
    image : str
        avgpb image
    radius : float, optional
        Radius beyond which to zero avgpb values (expressed as fraction of
        image width)
    output : str, optional
        Output image name. If None, add a 'z' to the end of the input filename

    """
    if type(radius) is str:
        radius = float(radius)

    pb = pim.image(image)
    pbdata = pb.getdata()
    (nx, ny) = pbdata.shape[2:]
    pbrad = radius * nx
    cy = ny / 2
    cx = nx / 2
    pbcounter = 0
    for j in range(nx, cx, -1):
        k = ny
        while np.sqrt((j-cx)**2+(k-cy)**2) > pbrad and k > cy:
            pbcounter += 1
            pbdata[:, :, k-1, j-1] = 0.
            pbdata[:, :, k-1, nx-j] = 0.
            pbdata[:, :, ny-k, j-1] = 0.
            pbdata[:, :, ny-k, nx-j] = 0.
            k -= 1
        if k == ny:
            break

    print pbcounter,'x 4 =',pbcounter*4,'zeros replaced'

    if output is None:
        while image[-1] == '/':
            image = image[: -1]
        outim = image + 'z'
    else:
        outim = output
    pout = pim.image(outim, values=pbdata, coordsys=pb.coordinates())


if __name__ == '__main__':
    descriptiontext = "Zero corners of avgpb images.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('image', help='Image to adjust')
    parser.add_argument('-o', '--output', help='output filename', default=None)
    parser.add_argument('-r', '--radius',help='Radius beyond which to zero avgpb values', default=0.5, type=float)

    args = parser.parse_args()
    main(args.image, radius=args.radius, output=args.output)
