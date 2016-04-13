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


def read_vertices(filename, cal_only=False):
    """
    Returns facet vertices

    Parameters
    ----------
    filename : str
        Filename of pickled file with direction vertices

    """
    with open(filename, 'r') as f:
        direction_dict = pickle.load(f)
    if cal_only:
        return direction_dict['vertices_cal']
    else:
        return direction_dict['vertices']


def main(fullskymodel, outmodel, vertices_file, cal_only=False, remove_cal=False):
    """
    Makes a makesourcedb sky model for components inside input polygon

    Parameters
    ----------
    fullskymodel : str
        Filename of makesourcedb sky model file containing the full-field model
    outmodel : str
        Filename of output sky model
    vertices_file : str
        Filename of pickled file with direction vertices that define polygon
    cal_only : bool, optional
        If True, only components within the calibrator region will be selected
    remove_cal : bool, optional
        If True, remove components from within the calibrator region from the
        full facet. Note: this option can only be activated if cal_only is
        False

    """
    if type(cal_only) is str:
        if cal_only.lower() == 'true':
            cal_only = True
        else:
            cal_only = False
    if type(remove_cal) is str:
        if remove_cal.lower() == 'true':
            remove_cal = True
        else:
            remove_cal = False
    if cal_only and remove_cal:
        print('cal_only and remove_cal cannot both be True')
        sys.exit(1)

    s = lsmtool.load(fullskymodel)
    vertices = read_vertices(vertices_file, cal_only=cal_only)

    # Select sources inside poly defined by vertices
    x, y, midRA, midDec = s._getXY()
    xv, yv = radec2xy(vertices[0], vertices[1], midRA, midDec)
    xyvertices = array([[xp, yp] for xp, yp in zip(xv, yv)])
    bbPath = mplPath.Path(xyvertices)
    inside = zeros(len(s), dtype=bool)
    for i in range(len(s)):
        inside[i] = bbPath.contains_point((x[i], y[i]))
    s.select(inside, force=True)

    if remove_cal and len(s) > 0:
        vertices = read_vertices(vertices_file, cal_only=True)

        # Select sources inside poly defined by vertices
        x, y, midRA, midDec = s._getXY()
        xv, yv = radec2xy(vertices[0], vertices[1], midRA, midDec)
        xyvertices = array([[xp, yp] for xp, yp in zip(xv, yv)])
        bbPath = mplPath.Path(xyvertices)
        inside = zeros(len(s), dtype=bool)
        for i in range(len(s)):
            inside[i] = bbPath.contains_point((x[i], y[i]))
        s.remove(inside, force=True)

    if len(s) == 0:
        print('No sources found for this facet')
        os.system('touch {0}'.format(outmodel))
    else:
        s.write(outmodel, clobber=True)


if __name__ == '__main__':
    descriptiontext = "Make a sky model for a facet.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('fullskymodel', help='name of the full skymodel')
    parser.add_argument('outmodel', help='name for the output')
    parser.add_argument('vertices_file', help='file containing facet vertices')
    parser.add_argument('-c', '--cal_only', help='return calibrator model only', type=bool, default=False)
    args = parser.parse_args()

    main(args.fullskymodel, args.outmodel, args.vertices_file, cal_only=args.cal_only)
