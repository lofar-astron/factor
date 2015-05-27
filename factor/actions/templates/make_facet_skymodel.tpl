#!/usr/bin/env python
"""
Script to make a sky model for a facet
"""
import lsmtool
from lsmtool.operations_lib import radec2xy
import matplotlib.path as mplPath
from numpy import array, zeros
import sys
import os


def read_vertices(filename):
    """
    Returns facet vertices
    """
    pass


if __name__ == '__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog <full_skymodel_filename> <out_skymodel_filename> <vert_filename> <ra> <dec>')
    opt.add_option('--cal_radius', help='radius in deg of calibrator', default=None)

    (options, args) = opt.parse_args()

    if len(args) != 3:
        opt.print_help()
        sys.exit(1)

    fullskymodel = args[0]
    outmodel = args[1]
    vertices = read_vertices(args[2])
    facet_ra = args[3]
    facet_dec = args[4]
    cal_radius_deg = options.cal_radius

    s = lsmtool.load(fullskymodel)

    if cal_radius_deg is not None:
        # Get calibrator model
        dist = s.getDistance(facet_ra, facet_dec)
        s.select(dist < cal_radius_deg)
    else:
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
    else:
        s.write(outmodel, clobber=True)
