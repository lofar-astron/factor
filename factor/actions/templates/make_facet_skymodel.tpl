import lsmtool
from lsmtool.operations_lib import radec2xy
import matplotlib.path as mplPath
from numpy import array, zeros
import sys
import os


fullskymodel = sys.argv[1]
outmodel = sys.argv[2]
cal_only = bool(sys.argv[3])
vertices = {{ vertices }}

s = lsmtool.load(fullskymodel)

if cal_only:
    # Get calibrator model
    dist = s.getDistance(self.facet_ra, self.facet_dec)
    s.select(dist < self.cal_radius_deg)
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
