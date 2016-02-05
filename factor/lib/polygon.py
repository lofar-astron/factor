"""
General polygon library

Contains the master class for all facet polygons.

Note: code based on:
http://code.activestate.com/recipes/578381-a-point-in-polygon-program-sw-sloan-algorithm/

"""
import numpy as np


class Polygon:
    """
    Generic polygon class

    Polygons are used to define the facet boundaries.

    Parameters
    ----------
    x : array
        A sequence of nodal x-coords.
    y : array
        A sequence of nodal y-coords.

    """
    def __init__(self, x, y):
        if len(x) != len(y):
            raise IndexError('x and y must be equally sized.')
        self.x = np.asfarray(x)
        self.y = np.asfarray(y)

        # Closes the polygon if were open
        x1, y1 = x[0], y[0]
        xn, yn = x[-1], y[-1]
        if x1 != xn or y1 != yn:
            self.x = np.concatenate((self.x, [x1]))
            self.y = np.concatenate((self.y, [y1]))

        # Anti-clockwise coordinates
        if _det(self.x, self.y) < 0:
            self.x = self.x[::-1]
            self.y = self.y[::-1]


    def is_inside(self, xpoint, ypoint, smalld=1e-12):
        """
        Check if point is inside a general polygon.

        An improved version of the algorithm of Nordbeck and Rydstedt.

        REF: SLOAN, S.W. (1985): A point-in-polygon program. Adv. Eng.
             Software, Vol 7, No. 1, pp 45-47.

        Parameters
        ----------
        xpoint : array or float
            The x-coord of the point to be tested.
        ypoint : array or float
            The y-coords of the point to be tested.
        smalld : float
            Tolerance within which point is considered to be on a side.

        Returns
        -------
        mindst : array or float
            The distance from the point to the nearest point of the polygon:
                If mindst < 0 then point is outside the polygon.
                If mindst = 0 then point in on a side of the polygon.
                If mindst > 0 then point is inside the polygon.

        """
        xpoint = np.asfarray(xpoint)
        ypoint = np.asfarray(ypoint)

        # Scalar to array
        if xpoint.shape is tuple():
            xpoint = np.array([xpoint], dtype=float)
            ypoint = np.array([ypoint], dtype=float)
            scalar = True
        else:
            scalar = False
        # Check consistency
        if xpoint.shape != ypoint.shape:
            raise IndexError('x and y  must be equally sized.')

        # If snear = True: Dist to nearest side < nearest vertex
        # If snear = False: Dist to nearest vertex < nearest side
        snear = np.ma.masked_all(xpoint.shape, dtype=bool)

        # Initialize arrays
        mindst = np.ones_like(xpoint, dtype=float) * np.inf
        j = np.ma.masked_all(xpoint.shape, dtype=int)
        x = self.x
        y = self.y
        n = len(x) - 1  # Number of sides/vertices defining the polygon

        # Loop over each side defining polygon
        for i in range(n):
            d = np.ones_like(xpoint, dtype=float) * np.inf

            # Start of side has coords (x1, y1)
            # End of side has coords (x2, y2)
            # Point has coords (xpoint, ypoint)
            x1 = x[i]
            y1 = y[i]
            x21 = x[i + 1] - x1
            y21 = y[i + 1] - y1
            x1p = x1 - xpoint
            y1p = y1 - ypoint

            # Points on infinite line defined by
            #     x = x1 + t * (x1 - x2)
            #     y = y1 + t * (y1 - y2)
            # where
            #     t = 0    at (x1, y1)
            #     t = 1    at (x2, y2)
            # Find where normal passing through (xpoint, ypoint) intersects
            # infinite line
            t = -(x1p * x21 + y1p * y21) / (x21 ** 2 + y21 ** 2)
            tlt0 = t < 0
            tle1 = (0 <= t) & (t <= 1)

            # Normal intersects side
            d[tle1] = ((x1p[tle1] + t[tle1] * x21) ** 2 +
                       (y1p[tle1] + t[tle1] * y21) ** 2)

            # Normal does not intersects side
            # Point is closest to vertex (x1, y1)
            # Compute square of distance to this vertex
            d[tlt0] = x1p[tlt0] ** 2 + y1p[tlt0] ** 2

            # Store distances
            mask = d < mindst
            mindst[mask] = d[mask]
            j[mask] = i

            # Point is closer to (x1, y1) than any other vertex or side
            snear[mask & tlt0] = False

            # Point is closer to this side than to any other side or vertex
            snear[mask & tle1] = True

        if np.ma.count(snear) != snear.size:
            raise IndexError('Error computing distances')
        mindst **= 0.5

        # Point is closer to its nearest vertex than its nearest side, check if
        # nearest vertex is concave.
        # If the nearest vertex is concave then point is inside the polygon,
        # else the point is outside the polygon.
        jo = j.copy()
        jo[j == 0] -= 1
        area = _det([x[j + 1], x[j], x[jo - 1]], [y[j + 1], y[j], y[jo - 1]])
        mindst[~snear] = np.copysign(mindst, area)[~snear]

        # Point is closer to its nearest side than to its nearest vertex, check
        # if point is to left or right of this side.
        # If point is to left of side it is inside polygon, else point is
        # outside polygon.
        area = _det([x[j], x[j + 1], xpoint], [y[j], y[j + 1], ypoint])
        mindst[snear] = np.copysign(mindst, area)[snear]

        # Point is on side of polygon
        mindst[np.fabs(mindst) < smalld] = 0

        # If input values were scalar then the output should be too
        if scalar:
            mindst = float(mindst)
        return mindst


def _det(xvert, yvert):
    """
    Compute twice the area of the triangle defined by points using the
    determinant formula.

    Parameters
    ----------
    xvert : array
        A vector of nodal x-coords.
    yvert : array
        A vector of nodal y-coords.

    Returns
    -------
    area : float
        Twice the area of the triangle defined by the points:
            area is positive if points define polygon in anticlockwise order.
            area is negative if points define polygon in clockwise order.
            area is zero if at least two of the points are concident or if
            all points are collinear.

    """
    xvert = np.asfarray(xvert)
    yvert = np.asfarray(yvert)
    x_prev = np.concatenate(([xvert[-1]], xvert[:-1]))
    y_prev = np.concatenate(([yvert[-1]], yvert[:-1]))
    return np.sum(yvert * x_prev - xvert * y_prev, axis=0)
