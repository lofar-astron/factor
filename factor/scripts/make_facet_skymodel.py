#! /usr/bin/env python
"""
Script to make a sky model for a facet
"""
import argparse
from argparse import RawTextHelpFormatter
import lsmtool
from lsmtool.operations_lib import radec2xy
import matplotlib.path as mplPath
from numpy import array, zeros
import sys
import os
import pickle


def read_vertices(filename):
    """
    Returns facet vertices
    """
    with open(filename, 'r') as f:
        direction_dict = pickle.load(f)
    return direction_dict['vertices']


def main(fullskymodel, outmodel, vertices_file, cal_only=False, facet_ra=0.0,
    facet_dec=0.0, cal_radius_deg=0.0):

    s = lsmtool.load(fullskymodel)
    vertices = read_vertices(vertices_file)

    # Get all facet sources
    x, y, midRA, midDec = s._getXY()
    xv, yv = radec2xy(vertices[0], vertices[1], midRA, midDec)
    xyvertices = array([[xp, yp] for xp, yp in zip(xv, yv)])
    bbPath = mplPath.Path(xyvertices)
    inside = zeros(len(s), dtype=bool)
    for i in range(len(s)):
        inside[i] = bbPath.contains_point((x[i], y[i]))
    s.select(inside, force=True)

    if len(s) == 0:
        print('No sources found for this facet')
        return

    if cal_only:
        # Get calibrator sources
        dist = s.getDistance(facet_ra, facet_dec)
        print(dist)
        print(facet_ra)
        print(facet_dec)
        print(midRA)
        print(midDec)
        print(cal_radius_deg)
        print(dist < cal_radius_deg)
        s.select(dist < cal_radius_deg)

    if len(s) == 0:
        print('No sources found for this facet')
    else:
        s.write(outmodel, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Make a sky model for a facet.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fullskymodel', help='name of the full skymodel')
    parser.add_argument('outmodel', help='name for the output')
    parser.add_argument('vertices_file', help='file containing facet vertices')
    parser.add_argument('-c', '--cal_only', help='return calibrator model only', type=bool, default=False)
    parser.add_argument('-r', '--facet_ra', help='RA of facet center in degrees', type=float, default=0.0)
    parser.add_argument('-d', '--facet_dec', help='Dec of facet center in degrees', type=float, default=0.0)
    parser.add_argument('-g', '--cal_radius_deg', help='Radius of calibrator in degrees', type=float, default=0.0)

    args = parser.parse_args()

    main(args.fullskymodel, args.outmodel, args.vertices_file, cal_only=args.cal_only,
        facet_ra=args.facet_ra, facet_dec=args.facet_dec, cal_radius_deg=args.cal_radius_deg)
