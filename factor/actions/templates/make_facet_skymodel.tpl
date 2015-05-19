import lsmtool
from lsmtool.operations_lib import radec2xy
import matplotlib.path as mplPath
from numpy import array, zeros
import sys
import os


fullskymodel = sys.argv[1]
outmodel = sys.argv[2]
cal_only = bool(int(sys.argv[3]))
vertices = {{ vertices }}
facet_ra = {{ ra }}
facet_dec = {{ dec }}
cal_radius_deg = {{ cal_radius }}

s = lsmtool.load(fullskymodel)

if cal_only:
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
